#pragma once

#include <fstream>
#include <ios>
#include <string>

namespace vcftools_ng::output {

// Stage every scientific artifact beside its final destination.  The main
// entry point publishes all staged files only after every writer has closed
// successfully.
std::string stage_path(const std::string& final_path);

std::ofstream open_stream(
    const std::string& final_path,
    std::ios::openmode mode = std::ios::out | std::ios::trunc);

void finish_stream(
    std::ofstream& stream, const std::string& final_path);

void commit_all();
void abandon_all() noexcept;

}  // namespace vcftools_ng::output
