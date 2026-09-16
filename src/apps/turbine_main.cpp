#include "turbine/simulation.hpp"
#include <iostream>
int main(int argc, char **argv) {
    try {
        if (argc < 3)
            throw std::runtime_error("Usage: aeroacoustics_turbine CASE.fst OUTPUT_DIRECTORY [duration] "
                                     "[--lookup-policy=clamp|error] [--solver=reference|scaled] "
                                     "[--observer-block-size=N]");
        turbine::RunOptions options;
        bool has_duration = false, has_policy = false, has_solver = false, has_block = false;
        for (int i = 3; i < argc; ++i) {
            const std::string option = argv[i];
            if (option.rfind("--lookup-policy=", 0) == 0 && !has_policy) {
                const auto policy = option.substr(16);
                if (policy != "clamp" && policy != "error")
                    throw std::invalid_argument("Lookup policy must be clamp or error");
                options.lookup_policy =
                    policy == "error" ? diagnostics::LookupPolicy::error : diagnostics::LookupPolicy::clamp;
                has_policy = true;
            } else if (option.rfind("--solver=", 0) == 0 && !has_solver) {
                const auto mode = option.substr(9);
                if (mode != "reference" && mode != "scaled")
                    throw std::invalid_argument("Solver must be reference or scaled");
                options.solver.mode =
                    mode == "scaled" ? turbine::SolverMode::scaled : turbine::SolverMode::reference;
                has_solver = true;
            } else if (option.rfind("--observer-block-size=", 0) == 0 && !has_block) {
                const auto count = option.substr(22);
                if (count.empty() || count.find_first_not_of("0123456789") != std::string::npos)
                    throw std::invalid_argument("Observer block size must be a positive integer");
                options.observer_block_size = std::stoull(count);
                has_block = true;
            } else if (!has_duration && option.rfind("--", 0) != 0) {
                std::size_t used = 0;
                options.duration = std::stod(option, &used);
                if (used != option.size())
                    throw std::invalid_argument("Invalid duration");
                has_duration = true;
            } else
                throw std::invalid_argument("Unknown or duplicate argument: " + option);
        }

        const auto result = turbine::run_case(argv[1], argv[2], options);
        if (result.lookup.calls)
            std::cerr << "Warning: " << result.lookup.calls
                      << " out-of-range lookup calls; endpoint values used. See lookup_diagnostics.csv "
                         "(includes BEM trial and UA internal queries).\n";
        std::cout << "C++ turbine completed " << result.duration << " s, " << result.steps << " steps, "
                  << result.acoustic_samples << " acoustic samples\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
