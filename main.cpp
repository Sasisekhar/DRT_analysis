// Without arguments: checks the implementation against the examples in the paper.
// With --json [horizon]: prints the example task set as JSON for draw_drt.py, e.g.
//     ./drt --json 60 | python3 draw_drt.py
#include "drt.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

using namespace drt;

#define CHECK(cond)                                                          \
    do {                                                                     \
        if (!(cond)) {                                                       \
            std::cerr << "CHECK failed: " #cond " (line " << __LINE__ << ")\n"; \
            std::exit(1);                                                    \
        }                                                                    \
    } while (0)

// Figure 2: frame-separated running example.
static Task figure2() {
    Task T("Figure 2");
    const auto v1 = T.add_vertex("v1", 2, 5);
    const auto v2 = T.add_vertex("v2", 1, 8);
    const auto v3 = T.add_vertex("v3", 3, 8);
    const auto v4 = T.add_vertex("v4", 5, 10);
    const auto v5 = T.add_vertex("v5", 1, 5);
    T.add_edge(v1, v2, 10);
    T.add_edge(v2, v3, 15);
    T.add_edge(v3, v1, 11);
    T.add_edge(v2, v4, 20);
    T.add_edge(v4, v2, 20);
    T.add_edge(v1, v5, 20);
    T.add_edge(v5, v4, 10);
    return T;
}

// Figure 8: arbitrary deadlines.
static Task figure8() {
    Task T("Figure 8");
    const auto v1 = T.add_vertex("v1", 5, 7);
    const auto v2 = T.add_vertex("v2", 1, 10);
    const auto v3 = T.add_vertex("v3", 2, 5);
    T.add_edge(v1, v2, 8);
    T.add_edge(v2, v3, 3);
    T.add_edge(v3, v1, 1000);
    return T;
}

static FeasibilityResult report(const char* label, const std::vector<Task>& tasks) {
    const auto r = check_feasibility(tasks);
    std::cout << label << ": U = " << r.utilization << ", " << status_name(r.status);
    if (r.bound > 0) std::cout << " (checked t <= " << r.bound << ")";
    if (r.witness) std::cout << ", witness t_f = " << *r.witness;
    std::cout << "\n  " << r.reason << "\n";
    return r;
}

int main(int argc, char** argv) {
    if (argc > 1 && std::string(argv[1]) == "--json") {
        const Time horizon = argc > 2 ? std::stoll(argv[2]) : 60;
        write_json(std::cout, {figure2(), figure8()}, horizon);
        return 0;
    }

    const Task t2 = figure2();
    const auto u2 = utilization(t2);
    std::cout << "Figure 2: U(T) = " << u2.e << "/" << u2.p << "\n";
    CHECK(u2.e == 1 && u2.p == 6);  // Example V.6: 3/18

    const auto dbf2 = demand_bound_function(t2, 60);
    for (Time t : {5, 10, 18, 43, 45, 51, 60}) std::cout << "  dbf(" << t << ") = " << dbf2(t) << "\n";
    CHECK(dbf2(43) == 9 && dbf2(45) == 9);  // Examples III.3 / Fig. 3
    CHECK(dbf2(51) >= 11);                  // Example IV.1: demand pair <11, 51>

    const Task t8 = figure8();
    CHECK(!t8.frame_separated());
    const auto dbf8 = demand_bound_function(t8, 20);
    std::cout << "Figure 8 (arbitrary deadlines):\n";
    for (Time t : {5, 7, 10, 16, 18}) std::cout << "  dbf(" << t << ") = " << dbf8(t) << "\n";
    CHECK(dbf8(16) == 7 && dbf8(18) == 8);  // Example VI.1: <7, 16> and <8, 18>

    std::cout << "\n";
    CHECK(report("{Fig. 2, Fig. 8}", {t2, t8}).status == FeasibilityResult::Status::Feasible);

    Task heavy("heavy");
    const auto h = heavy.add_vertex("h", 4, 5);
    heavy.add_edge(h, h, 20);
    const auto r = report("{Fig. 2, heavy}", {t2, heavy});
    CHECK(r.status == FeasibilityResult::Status::Infeasible && r.witness == 5);  // 2 + 4 > 5

    Task full("full");
    const auto f = full.add_vertex("f", 3, 3);
    full.add_edge(f, f, 3);
    CHECK(report("{Fig. 2, full}", {t2, full}).status == FeasibilityResult::Status::Infeasible);

    std::cout << "\nAll checks passed.\n";
}
