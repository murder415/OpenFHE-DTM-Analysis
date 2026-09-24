#include "openfhe_dtm/OpenFheBackend.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
    using namespace openfhe_dtm;
    try {
        OpenFheBackend he(8, {}, true);
        const std::vector<double> input{0.125, -0.25, 0.5, 0.75, -1.0, 1.25, 1.5, -1.75};
        const auto encrypted = he.encrypt(he.encode(input));
        const auto refreshed = he.bootstrap(encrypted);
        const auto output = he.decrypt(refreshed).values;

        double max_abs_error = 0.0;
        for (std::size_t i = 0; i < input.size(); ++i)
            max_abs_error = std::max(max_abs_error, std::abs(output[i] - input[i]));

        std::cout << "CKKS bootstrap executed: slots=" << he.slotCount()
                  << " max_abs_error=" << max_abs_error << "\n";
        if (max_abs_error > 0.01)
            throw std::runtime_error("bootstrap error exceeded 0.01");
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "bootstrap validation failed: " << e.what() << "\n";
        return 1;
    }
}
