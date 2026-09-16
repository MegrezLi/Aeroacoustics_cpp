#include "turbine/batch.hpp"
#include <algorithm>
#include <atomic>
#include <cctype>
#include <cwctype>
#include <thread>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
#ifdef AERO_USE_MKL
#include <mkl.h>
#endif
namespace turbine {
namespace {
std::filesystem::path normalized(const std::filesystem::path &path) {
#ifdef _WIN32
    // Resolve the existing ancestor, including junctions, before appending new folders.
    // MinGW's canonical() can report ENOENT for existing Windows paths.
    auto ancestor = std::filesystem::absolute(path).lexically_normal();
    std::filesystem::path suffix;
    while (!std::filesystem::exists(ancestor)) {
        if (ancestor == ancestor.root_path())
            throw std::invalid_argument("Batch path has no existing root");
        suffix = ancestor.filename() / suffix;
        ancestor = ancestor.parent_path();
    }
    HANDLE handle = CreateFileW(ancestor.c_str(), 0, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
                                nullptr, OPEN_EXISTING, FILE_FLAG_BACKUP_SEMANTICS, nullptr);
    if (handle == INVALID_HANDLE_VALUE)
        throw std::system_error(GetLastError(), std::system_category(), "Cannot resolve batch path");
    struct CloseHandleOnExit {
        HANDLE handle;
        ~CloseHandleOnExit() { CloseHandle(handle); }
    } close{handle};
    const DWORD size = GetFinalPathNameByHandleW(handle, nullptr, 0, FILE_NAME_NORMALIZED);
    if (!size)
        throw std::system_error(GetLastError(), std::system_category(), "Cannot resolve batch path");
    std::wstring resolved(size, L'\0');
    const DWORD used = GetFinalPathNameByHandleW(handle, resolved.data(), size, FILE_NAME_NORMALIZED);
    if (!used || used >= size)
        throw std::runtime_error("Batch path resolution failed");
    resolved.resize(used);
    if (resolved.rfind(L"\\\\?\\UNC\\", 0) == 0)
        resolved = L"\\\\" + resolved.substr(8);
    else if (resolved.rfind(L"\\\\?\\", 0) == 0)
        resolved.erase(0, 4);
    auto result = (std::filesystem::path(resolved) / suffix).lexically_normal();
    auto text = result.native();
    std::transform(text.begin(), text.end(), text.begin(), [](wchar_t c) { return std::towlower(c); });
    result = text;
#else
    auto result = std::filesystem::weakly_canonical(std::filesystem::absolute(path));
#endif
    while (result != result.root_path() && result.filename().empty())
        result = result.parent_path();
    return result;
}
bool contains(const std::filesystem::path &parent, const std::filesystem::path &child) {
    auto a = parent.begin(), b = child.begin();
    for (; a != parent.end() && b != child.end() && *a == *b; ++a, ++b) {
    }
    return a == parent.end();
}
struct LocalMathThreads {
#ifdef AERO_USE_MKL
    int previous = mkl_set_num_threads_local(1);
    ~LocalMathThreads() { mkl_set_num_threads_local(previous); }
#endif
};
} // namespace
std::vector<CaseOutcome> run_cases(const std::vector<CaseJob> &jobs, std::size_t workers) {
    if (!workers)
        throw std::invalid_argument("Batch worker count must be positive");
    std::vector<std::filesystem::path> outputs;
    for (const auto &job : jobs) {
        if (job.input.empty() || job.output.empty())
            throw std::invalid_argument("Batch input and output paths must not be empty");
        auto output = normalized(job.output);
        for (const auto &other : outputs)
            if (contains(other, output) || contains(output, other))
                throw std::invalid_argument("Batch output directories overlap: " + job.output.string());
        outputs.push_back(std::move(output));
    }
    // Do not allow an output directory to contain another job's primary input.
    for (const auto &job : jobs)
        for (const auto &output : outputs)
            if (contains(output, normalized(job.input)))
                throw std::invalid_argument("Batch output contains a case input");
    std::vector<CaseOutcome> results(jobs.size());
    std::atomic<std::size_t> next{0};
    const auto execute = [&] {
        LocalMathThreads math_threads;
        while (true) {
            const auto i = next.fetch_add(1, std::memory_order_relaxed);
            if (i >= jobs.size())
                return;
            try {
                results[i].summary = run_case(jobs[i].input, jobs[i].output, jobs[i].options);
            } catch (const std::exception &e) {
                results[i].error = e.what();
            } catch (...) {
                results[i].error = "Unknown case failure";
            }
        }
    };
    workers = std::min(workers, jobs.size());
    std::vector<std::thread> pool;
    pool.reserve(workers > 0 ? workers - 1 : 0);
    try {
        for (std::size_t i = 1; i < workers; ++i)
            pool.emplace_back(execute);
    } catch (...) {
        for (auto &thread : pool)
            thread.join();
        throw;
    }
    execute();
    for (auto &thread : pool)
        thread.join();
    return results;
}
} // namespace turbine
