# Benchmark demos (W18 plan — see docs/benchmark_plan.md).
# Runnable smoke tiers by default; full sweeps via `benchmark_suite --full`.
add_executable(benchmark_suite
    benchmark_suite.cpp
    cases_fdm3.cpp
    cases_modules.cpp
)
target_link_libraries(benchmark_suite PRIVATE exd::physics)
