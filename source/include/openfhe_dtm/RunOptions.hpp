#pragma once

#include <stdexcept>
#include <string>

namespace openfhe_dtm {
namespace detail {

inline bool endsWithTdbExtension(std::string path) {
    for (char& c : path) {
        if (c >= 'A' && c <= 'Z') {
            c = static_cast<char>(c - 'A' + 'a');
        }
    }
    return path.size() >= 4 && path.substr(path.size() - 4) == ".tdb";
}

} // namespace detail

struct RunOptions {
    bool synthetic = false;
    std::string tdb_path;
    std::string key_dir;
    std::string block_dir;
    std::string preset;
    std::string device;
    int level = 0;
};

inline bool originalDataRequested(const RunOptions& options) {
    return !options.tdb_path.empty();
}

inline void validateRunOptions(const RunOptions& options) {
    if (options.synthetic && originalDataRequested(options)) {
        throw std::invalid_argument("--synthetic cannot be combined with --tdb");
    }
    if (originalDataRequested(options) && options.tdb_path.size() < 5) {
        throw std::invalid_argument("--tdb must point to an original-compatible .tdb file");
    }
    if (originalDataRequested(options) && !detail::endsWithTdbExtension(options.tdb_path)) {
        throw std::invalid_argument("--tdb must point to an original-compatible .tdb file");
    }
}

inline RunOptions parseRunOptions(int argc, char** argv) {
    RunOptions options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto needValue = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                throw std::invalid_argument(std::string(name) + " requires a value");
            }
            return argv[++i];
        };

        if (arg == "--synthetic") {
            options.synthetic = true;
        } else if (arg == "--tdb") {
            options.tdb_path = needValue("--tdb");
        } else if (arg == "--key-dir") {
            options.key_dir = needValue("--key-dir");
        } else if (arg == "--block-dir") {
            options.block_dir = needValue("--block-dir");
        } else if (arg == "--preset") {
            options.preset = needValue("--preset");
        } else if (arg == "--device") {
            options.device = needValue("--device");
        } else if (arg == "--level") {
            options.level = std::stoi(needValue("--level"));
        } else if (arg == "--config") {
            (void)needValue("--config");
        } else {
            throw std::invalid_argument("unknown option: " + arg);
        }
    }
    validateRunOptions(options);
    return options;
}

} // namespace openfhe_dtm
