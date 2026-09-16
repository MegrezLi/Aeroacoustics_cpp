#include "turbine/batch.hpp"
#include <iostream>
int main(int argc, char **argv) {
    try {
        if (argc < 4 || argc % 2 != 0)
            throw std::invalid_argument(
                "Usage: aeroacoustics_batch WORKERS CASE.fst OUTPUT [CASE.fst OUTPUT ...]");
        const std::string count = argv[1];
        if (count.empty() || count.find_first_not_of("0123456789") != std::string::npos)
            throw std::invalid_argument("WORKERS must be a positive integer");
        std::vector<turbine::CaseJob> jobs;
        for (int i = 2; i < argc; i += 2)
            jobs.push_back({argv[i], argv[i + 1], {}});
        const auto results = turbine::run_cases(jobs, std::stoull(count));
        bool failed = false;
        for (std::size_t i = 0; i < results.size(); ++i) {
            std::cout << jobs[i].input.string() << ": ";
            if (results[i].summary)
                std::cout << "completed " << results[i].summary->steps << " steps\n";
            else {
                std::cout << results[i].error << '\n';
                failed = true;
            }
        }
        return failed ? 1 : 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
