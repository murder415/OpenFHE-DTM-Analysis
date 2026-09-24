#include "openfhe_dtm/TdbDataset.hpp"

#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: make_validation_tdb <output.tdb>\n";
        return 2;
    }

    try {
        openfhe_dtm::writeValidationTdbDataset(argv[1]);
        std::cout << "wrote validation .tdb: " << argv[1] << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "failed to write validation .tdb: " << e.what() << "\n";
        return 1;
    }
}
