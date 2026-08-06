#include "run_logger.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdexcept>

#include <sys/stat.h>
#include <unistd.h>

namespace vcftools_ng {
namespace {

double timeval_seconds(const timeval& value) {
    return static_cast<double>(value.tv_sec) +
           static_cast<double>(value.tv_usec) / 1'000'000.0;
}

std::string shell_quote(const std::string& value) {
    if (value.empty()) {
        return "''";
    }
    bool safe = true;
    for (const unsigned char character : value) {
        if (!(std::isalnum(character) || character == '_' ||
              character == '-' || character == '.' ||
              character == '/' || character == ':' ||
              character == '=' || character == ',')) {
            safe = false;
            break;
        }
    }
    if (safe) {
        return value;
    }
    std::string result = "'";
    for (const char character : value) {
        if (character == '\'') {
            result += "'\\''";
        } else {
            result += character;
        }
    }
    result += '\'';
    return result;
}

std::string command_line(int argc, char** argv) {
    std::ostringstream command;
    for (int index = 0; index < argc; ++index) {
        if (index != 0) {
            command << ' ';
        }
        command << shell_quote(argv[index]);
    }
    return command.str();
}

}  // namespace

class RunLogger::TeeBuffer final : public std::streambuf {
public:
    TeeBuffer(std::streambuf* terminal, std::streambuf* file)
        : terminal_(terminal), file_(file) {}

protected:
    int overflow(int character) override {
        if (traits_type::eq_int_type(
                character, traits_type::eof())) {
            return sync() == 0
                       ? traits_type::not_eof(character)
                       : traits_type::eof();
        }
        std::lock_guard lock(mutex_);
        const char value = traits_type::to_char_type(character);
        const bool terminal_ok =
            !traits_type::eq_int_type(
                terminal_->sputc(value), traits_type::eof());
        if (!file_failed_ && file_ != nullptr &&
            traits_type::eq_int_type(
                file_->sputc(value), traits_type::eof())) {
            disable_file();
        }
        return terminal_ok ? character : traits_type::eof();
    }

    std::streamsize xsputn(
        const char* data, std::streamsize size) override {
        std::lock_guard lock(mutex_);
        const std::streamsize terminal_written =
            terminal_->sputn(data, size);
        if (!file_failed_ && file_ != nullptr &&
            file_->sputn(data, size) != size) {
            disable_file();
        }
        return terminal_written;
    }

    int sync() override {
        std::lock_guard lock(mutex_);
        const int terminal_status = terminal_->pubsync();
        if (!file_failed_ && file_ != nullptr &&
            file_->pubsync() != 0) {
            disable_file();
        }
        return terminal_status;
    }

public:
    [[nodiscard]] bool file_complete() const noexcept {
        return !file_failed_;
    }

private:
    void disable_file() noexcept {
        if (file_failed_) {
            return;
        }
        file_failed_ = true;
        static constexpr char warning[] =
            "\nWarning: log file write failed; continuing with stderr "
            "only. The log file is incomplete.\n";
        (void)terminal_->sputn(
            warning, static_cast<std::streamsize>(sizeof(warning) - 1));
        (void)terminal_->pubsync();
    }

    std::streambuf* terminal_;
    std::streambuf* file_;
    std::mutex mutex_;
    bool file_failed_ = false;
};

std::string format_timestamp(
    std::chrono::system_clock::time_point value) {
    const std::time_t timestamp =
        std::chrono::system_clock::to_time_t(value);
    std::tm local {};
    localtime_r(&timestamp, &local);
    char date[40] {};
    if (std::strftime(
            date, sizeof(date), "%Y-%m-%dT%H:%M:%S%z",
            &local) == 0) {
        return "unknown";
    }
    std::string result(date);
    if (result.size() >= 5) {
        result.insert(result.size() - 2, ":");
    }
    return result;
}

RunLogger::RunLogger(
    std::optional<std::string> file_path,
    const std::string& version, int argc, char** argv,
    std::chrono::steady_clock::time_point steady_start,
    std::chrono::system_clock::time_point wall_start)
    : steady_start_(steady_start) {
    (void)getrusage(RUSAGE_SELF, &usage_start_);
    (void)getrusage(RUSAGE_CHILDREN, &child_usage_start_);
    original_stderr_ = std::cerr.rdbuf();
    if (file_path.has_value()) {
        file_path_ = *file_path;
        struct stat stderr_status {};
        struct stat path_status {};
        file_via_stderr_ =
            fstat(STDERR_FILENO, &stderr_status) == 0 &&
            stat(file_path_.c_str(), &path_status) == 0 &&
            S_ISREG(stderr_status.st_mode) &&
            stderr_status.st_dev == path_status.st_dev &&
            stderr_status.st_ino == path_status.st_ino;
        if (!file_via_stderr_) {
            file_.emplace(
                file_path_,
                std::ios::out | std::ios::trunc);
            if (!*file_) {
                file_.reset();
                throw std::runtime_error(
                    "Could not create log file: " + file_path_);
            }
            tee_ = std::make_unique<TeeBuffer>(
                original_stderr_, file_->rdbuf());
            std::cerr.rdbuf(tee_.get());
        }
    }

    std::cerr
        << "vcftools-ng " << version << "\n"
        << "Log format: vcftools-ng-text-v1\n"
        << "Start time: " << format_timestamp(wall_start) << "\n"
        << "Command: " << command_line(argc, argv) << "\n";
    std::error_code directory_error;
    const auto directory =
        std::filesystem::current_path(directory_error);
    std::cerr
        << "Working directory: "
        << (directory_error ? "unknown" : directory.string())
        << "\n"
        << "Log file: "
        << (file_path.has_value()
                ? file_path_
                : "disabled (--no-log-file)")
        << "\n\n";
}

RunLogger::~RunLogger() {
    if (!finished_) {
        try {
            finish("failed");
        } catch (...) {
        }
    }
    if (original_stderr_ != nullptr) {
        std::cerr.rdbuf(original_stderr_);
    }
}

void RunLogger::finish(const std::string& status) {
    if (finished_) {
        return;
    }
    struct rusage usage_end {};
    struct rusage child_usage_end {};
    (void)getrusage(RUSAGE_SELF, &usage_end);
    (void)getrusage(RUSAGE_CHILDREN, &child_usage_end);
    const double wall_seconds =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - steady_start_)
            .count();
    const double user_seconds =
        timeval_seconds(usage_end.ru_utime) -
            timeval_seconds(usage_start_.ru_utime) +
        timeval_seconds(child_usage_end.ru_utime) -
            timeval_seconds(child_usage_start_.ru_utime);
    const double system_seconds =
        timeval_seconds(usage_end.ru_stime) -
            timeval_seconds(usage_start_.ru_stime) +
        timeval_seconds(child_usage_end.ru_stime) -
            timeval_seconds(child_usage_start_.ru_stime);
    const long peak_rss =
        std::max(usage_end.ru_maxrss, child_usage_end.ru_maxrss);
    const double cpu_percent =
        wall_seconds > 0.0
            ? (user_seconds + system_seconds) /
                  wall_seconds * 100.0
            : 0.0;

    std::cerr
        << "\nWall time: " << std::fixed << std::setprecision(3)
        << wall_seconds << " seconds\n"
        << "User CPU time: " << user_seconds << " seconds\n"
        << "System CPU time: " << system_seconds << " seconds\n"
        << "Average CPU: " << cpu_percent << "%\n"
        << "Peak RSS: " << peak_rss << " KiB\n"
        << "Exit status: " << status << "\n"
        << "End time: "
        << format_timestamp(std::chrono::system_clock::now())
        << "\n";
    std::cerr.flush();
    if (tee_ != nullptr && !tee_->file_complete()) {
        std::cerr << "Log status: incomplete (file write failed)\n";
        std::cerr.flush();
    }
    finished_ = true;
}

bool RunLogger::file_enabled() const noexcept {
    return file_.has_value() || file_via_stderr_;
}

const std::string& RunLogger::file_path() const noexcept {
    return file_path_;
}

}  // namespace vcftools_ng
