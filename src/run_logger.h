#pragma once

#include <chrono>
#include <fstream>
#include <memory>
#include <optional>
#include <streambuf>
#include <string>

#include <sys/resource.h>

namespace vcftools_ng {

class RunLogger {
public:
    RunLogger(
        std::optional<std::string> file_path,
        const std::string& version, int argc, char** argv,
        std::chrono::steady_clock::time_point steady_start,
        std::chrono::system_clock::time_point wall_start);
    ~RunLogger();

    RunLogger(const RunLogger&) = delete;
    RunLogger& operator=(const RunLogger&) = delete;

    void finish(const std::string& status);
    [[nodiscard]] bool file_enabled() const noexcept;
    [[nodiscard]] const std::string& file_path() const noexcept;

private:
    class TeeBuffer;

    std::optional<std::ofstream> file_;
    std::unique_ptr<TeeBuffer> tee_;
    std::streambuf* original_stderr_ = nullptr;
    std::string file_path_;
    std::chrono::steady_clock::time_point steady_start_;
    struct rusage usage_start_ {};
    struct rusage child_usage_start_ {};
    bool file_via_stderr_ = false;
    bool finished_ = false;
};

std::string format_timestamp(
    std::chrono::system_clock::time_point value);

}  // namespace vcftools_ng
