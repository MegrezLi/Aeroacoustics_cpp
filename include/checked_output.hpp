#pragma once
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace diagnostics {
inline void flush_checked(std::ostream &stream, const std::string &path) {
    try {
        stream.flush();
        if (!stream)
            throw std::runtime_error("stream failed");
    } catch (const std::exception &e) {
        throw std::runtime_error("Cannot flush output " + path + ": " + e.what());
    }
}
class CheckedOutput {
    std::ofstream stream_;
    std::string path_;

  public:
    void open(const std::filesystem::path &path) {
        path_ = path.string();
        try {
            stream_.exceptions(std::ios::badbit | std::ios::failbit);
            stream_.open(path);
        } catch (const std::exception &e) {
            throw std::runtime_error("Cannot open output " + path_ + ": " + e.what());
        }
    }
    template <class T> CheckedOutput &operator<<(const T &value) {
        try {
            stream_ << value;
        } catch (const std::exception &e) {
            throw std::runtime_error("Cannot write output " + path_ + ": " + e.what());
        }
        return *this;
    }
    void finish() {
        flush_checked(stream_, path_);
        try {
            stream_.close();
        } catch (const std::exception &e) {
            throw std::runtime_error("Cannot close output " + path_ + ": " + e.what());
        }
    }
};
} // namespace diagnostics
