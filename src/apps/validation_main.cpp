#include "turbine/validation.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
int main(int argc, char **argv) {
    try {
        const std::string command = argc > 1 ? argv[1] : "";
        auto factor = [&]() {
            std::size_t n = 0;
            const std::string s = argv[5];
            double k = std::stod(s, &n);
            if (n != s.size() || !std::isfinite(k) || k <= 0)
                throw std::invalid_argument("Invalid coverage factor");
            return k;
        };
        using namespace turbine::validation;
        if (command == "compare" && argc == 6)
            compare_files(argv[2], argv[3], argv[4], factor());
        else if (command == "budget" && argc == 6)
            budget_files(argv[2], argv[3], argv[4], factor());
        else if (command == "ensemble" && argc == 4)
            ensemble_file(argv[2], argv[3]);
        else if (command == "import-map" && argc == 6)
            import_receiver_map(argv[2], argv[3], argv[4], argv[5]);
        else
            throw std::invalid_argument(
                "Usage:\n  aeroacoustics_validate compare OBS.csv PRED.csv NEW_OUTPUT_DIR K\n  "
                "aeroacoustics_validate budget PERTURBATIONS.csv CORRELATIONS.csv NEW_OUTPUT_DIR K\n  "
                "aeroacoustics_validate ensemble SAMPLES.csv NEW_OUTPUT_DIR\n  aeroacoustics_validate "
                "import-map RECEIVER_MAP.csv CONTEXT.csv NEW_PRED.csv PROVENANCE");
        std::cout << "Validation analysis completed (no automatic field-accuracy or compliance claim).\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
