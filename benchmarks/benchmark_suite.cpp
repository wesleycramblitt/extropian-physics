// Benchmark demo runner: `benchmark_suite [--case NAME] [--list] [--grid N]
// [--full]`.  Default: run every case in smoke tier (seconds each).
//
// The full sweeps described in docs/benchmark_plan.md are enabled with
// `--full` and are intended for a capable machine; smoke is CI-sized.

#include "cases.hpp"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    const auto cases = bench::all_cases();

    std::string want;
    bool full = false;
    int grid = 0;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--case") == 0 && i + 1 < argc) want = argv[++i];
        else if (std::strcmp(argv[i], "--list") == 0) {
            for (const auto& c : cases)
                std::cout << c.name << " — " << c.family << " (" << c.anchor << ")\n";
            return 0;
        } else if (std::strcmp(argv[i], "--full") == 0) full = true;
        else if (std::strcmp(argv[i], "--grid") == 0 && i + 1 < argc)
            grid = std::atoi(argv[++i]);
        else {
            std::cerr << "usage: benchmark_suite [--case NAME] [--full] [--grid N] [--list]\n";
            return 2;
        }
    }

    bench::RunSpec spec;
    spec.full = full;
    spec.grid = grid;

    int ran = 0;
    for (const auto& c : cases) {
        if (!want.empty() && want != c.name) continue;
        std::cout << "──── " << c.name << " (" << c.family << ") — anchor: "
                  << c.anchor << (full ? " [FULL]" : " [smoke]") << " ────\n";
        c.run(spec);
        ++ran;
    }
    std::cout << "ran " << ran << " case(s); full=" << (full ? "yes" : "no")
              << " grid=" << grid << "\n";
    return 0;
}
