#include "output_transaction.h"

#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <system_error>
#include <vector>

#include <unistd.h>

namespace vcftools_ng::output {
namespace {

struct Entry {
    std::string final_path;
    std::string staged_path;
    std::string backup_path;
    bool backup_created = false;
    bool published = false;
};

std::mutex registry_mutex;
std::vector<Entry> registry;
std::atomic<unsigned long> next_identifier{0};

std::string private_path(
    const std::string& path, const char* role) {
    return path + ".vcftools-ng." + role + "." +
           std::to_string(static_cast<long long>(::getpid())) + "." +
           std::to_string(next_identifier.fetch_add(1));
}

std::string unused_private_path(
    const std::string& path, const char* role) {
    for (unsigned attempt = 0; attempt < 1024; ++attempt) {
        const std::string candidate = private_path(path, role);
        std::error_code error;
        const auto status =
            std::filesystem::symlink_status(candidate, error);
        if ((!error && !std::filesystem::exists(status)) ||
            error == std::errc::no_such_file_or_directory) {
            return candidate;
        }
        if (error) {
            throw std::runtime_error(
                "Could not inspect private output path " + candidate +
                ": " + error.message());
        }
    }
    throw std::runtime_error(
        "Could not allocate a private output path beside " + path);
}

bool remove_if_present(const std::string& path) noexcept {
    std::error_code error;
    std::filesystem::remove(path, error);
    if (error) {
        std::cerr
            << "Warning: could not remove private output path "
            << path << ": " << error.message() << "\n";
        return false;
    }
    return true;
}

void restore_backups(std::vector<Entry>& entries) noexcept {
    const bool inject_restore_failure = [] {
        const char* value =
            std::getenv("VCFTOOLS_NG_TEST_FAIL_OUTPUT_RESTORE");
        return value != nullptr && *value != '\0' &&
               std::string(value) != "0";
    }();
    bool injected = false;
    for (auto iterator = entries.rbegin(); iterator != entries.rend();
         ++iterator) {
        if (iterator->published) {
            if (!remove_if_present(iterator->final_path)) {
                std::cerr
                    << "Error: output recovery could not remove the newly "
                       "published destination; preserved backup, if any: "
                    << iterator->backup_path << "\n";
                continue;
            }
            iterator->published = false;
        }
        if (!iterator->backup_created) {
            continue;
        }
        std::error_code error;
        if (inject_restore_failure && !injected) {
            injected = true;
            error = std::make_error_code(std::errc::io_error);
        } else {
            std::filesystem::rename(
                iterator->backup_path, iterator->final_path, error);
        }
        if (error) {
            std::cerr
                << "Error: output recovery is incomplete; restore "
                << iterator->backup_path << " to "
                << iterator->final_path << " manually: "
                << error.message() << "\n";
            continue;
        }
        iterator->backup_created = false;
    }
}

}  // namespace

std::string stage_path(const std::string& final_path) {
    std::lock_guard lock(registry_mutex);
    for (const auto& entry : registry) {
        if (entry.final_path == final_path) {
            throw std::runtime_error(
                "Output destination requested more than once: " +
                final_path);
        }
    }
    Entry entry;
    entry.final_path = final_path;
    entry.staged_path = unused_private_path(final_path, "tmp");
    registry.push_back(entry);
    return entry.staged_path;
}

std::ofstream open_stream(
    const std::string& final_path, std::ios::openmode mode) {
    const std::string staged = stage_path(final_path);
    std::ofstream stream;
    stream.exceptions(std::ios::badbit | std::ios::failbit);
    try {
        stream.open(staged, mode);
    } catch (const std::ios_base::failure& error) {
        throw std::runtime_error(
            "Could not open staged output for " + final_path + ": " +
            error.what());
    }
    return stream;
}

void finish_stream(
    std::ofstream& stream, const std::string& final_path) {
    if (!stream.is_open()) {
        return;
    }
    try {
        stream.flush();
        stream.close();
    } catch (const std::ios_base::failure& error) {
        throw std::runtime_error(
            "Could not finish output " + final_path + ": " +
            error.what());
    }
    if (!stream) {
        throw std::runtime_error(
            "Could not finish output " + final_path);
    }
}

void commit_all() {
    std::lock_guard lock(registry_mutex);
    if (const char* injected =
            std::getenv("VCFTOOLS_NG_TEST_FAIL_OUTPUT_COMMIT");
        injected != nullptr && *injected != '\0' &&
        std::string(injected) != "0") {
        throw std::runtime_error(
            "Injected scientific-output commit failure");
    }
    std::size_t injected_failure_after = 0;
    if (const char* injected =
            std::getenv("VCFTOOLS_NG_TEST_FAIL_OUTPUT_COMMIT_AFTER");
        injected != nullptr && *injected != '\0') {
        char* end = nullptr;
        const unsigned long long parsed =
            std::strtoull(injected, &end, 10);
        if (end != injected && *end == '\0' && parsed > 0 &&
            parsed <= std::numeric_limits<std::size_t>::max()) {
            injected_failure_after = static_cast<std::size_t>(parsed);
        }
    }

    try {
        for (auto& entry : registry) {
            std::error_code status_error;
            const auto status = std::filesystem::symlink_status(
                entry.final_path, status_error);
            if (status_error &&
                status_error != std::errc::no_such_file_or_directory) {
                throw std::runtime_error(
                    "Could not inspect output destination " +
                    entry.final_path + ": " + status_error.message());
            }
            if (!status_error && std::filesystem::exists(status)) {
                if (std::filesystem::is_directory(status)) {
                    throw std::runtime_error(
                        "Output destination is a directory: " +
                        entry.final_path);
                }
                entry.backup_path = unused_private_path(
                    entry.final_path, "backup");
                std::filesystem::rename(
                    entry.final_path, entry.backup_path);
                entry.backup_created = true;
            }
        }
        std::size_t published_count = 0;
        for (auto& entry : registry) {
            std::filesystem::rename(
                entry.staged_path, entry.final_path);
            entry.published = true;
            ++published_count;
            if (injected_failure_after == published_count) {
                throw std::runtime_error(
                    "Injected scientific-output partial commit failure");
            }
        }
        for (auto& entry : registry) {
            if (entry.backup_created) {
                if (remove_if_present(entry.backup_path)) {
                    entry.backup_created = false;
                }
            }
        }
        registry.clear();
    } catch (...) {
        restore_backups(registry);
        for (const auto& entry : registry) {
            remove_if_present(entry.staged_path);
        }
        registry.clear();
        throw;
    }
}

void abandon_all() noexcept {
    std::lock_guard lock(registry_mutex);
    restore_backups(registry);
    for (const auto& entry : registry) {
        remove_if_present(entry.staged_path);
    }
    registry.clear();
}

}  // namespace vcftools_ng::output
