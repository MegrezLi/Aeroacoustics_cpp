#include "turbine/farm.hpp"
#include <iostream>
int main(int argc, char **argv) {
    try {
        if (argc != 3)
            throw std::invalid_argument("Usage: aeroacoustics_farm FARM.dat NEW_OUTPUT_DIRECTORY");
        turbine::run_farm(turbine::FarmOptions::read(argv[1]), argv[2]);
        std::cout << "C++ engineering farm completed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
