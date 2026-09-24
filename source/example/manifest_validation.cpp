#include "openfhe_dtm/OpenFheBackend.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    using namespace openfhe_dtm;
    if (argc != 2) {
        std::cerr << "usage: manifest_validation <key-directory>\n";
        return 2;
    }
    try {
        const std::filesystem::path key_directory(argv[1]);
        OpenFheBackend created(8, key_directory, false);
        if (!std::filesystem::is_regular_file(key_directory / "manifest.json"))
            throw std::runtime_error("manifest.json was not created");

        OpenFheBackend restored(8, key_directory, false);
        bool mismatch_rejected = false;
        try {
            OpenFheBackend incompatible(16, key_directory, false);
        } catch (const std::exception&) {
            mismatch_rejected = true;
        }
        if (!mismatch_rejected)
            throw std::runtime_error("incompatible key manifest was accepted");

        std::cout << "Key manifest: created, restored, mismatch rejected\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "manifest validation failed: " << e.what() << "\n";
        return 1;
    }
}
