#ifndef DRT_HPP
#define DRT_HPP

// Feasibility analysis for the Digraph Real-Time (DRT) task model:
// M. Stigge, P. Ekberg, N. Guan, W. Yi, "The Digraph Real-Time Task Model".
// Header-only: just #include "drt.hpp".

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <cstdio>
#include <optional>
#include <ostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace drt {

using Time = std::int64_t;

struct Edge {
    std::size_t to;  // index of the target vertex
    Time        p;   // minimum inter-release separation p(u, v)
};

struct Vertex {
    std::string       id;
    Time              e;    // worst-case execution time e(v)
    Time              d;    // relative deadline d(v)
    std::vector<Edge> out;  // outgoing edges
};

// One DRT task, i.e. the graph G(T). Edges refer to vertices by index, so
// adding vertices never invalidates existing edges.
class Task {
public:
    explicit Task(std::string name = "") : name_(std::move(name)) {}

    std::size_t add_vertex(std::string id, Time e, Time d);
    void        add_edge(std::size_t from, std::size_t to, Time p);

    const std::string&         name() const { return name_; }
    const std::vector<Vertex>& vertices() const { return vertices_; }

    // Frame Separation property: d(u) <= p(u, v) for every edge (u, v).
    bool frame_separated() const;
    // e_sum(T): sum of e(v) over all vertices (Lemma V.2).
    Time e_sum() const;

private:
    std::string         name_;
    std::vector<Vertex> vertices_;
};

// U(T) = e / p as an exact fraction in lowest terms; 0/1 if G(T) has no cycle.
struct Utilization {
    Time   e = 0;
    Time   p = 1;
    double value() const { return static_cast<double>(e) / static_cast<double>(p); }
};

// U(T) via utilization triples (Section V-A, Definition V.7).
Utilization utilization(const Task& task);

// dbf_T as a step function, exact for all t <= horizon.
// Frame-separated tasks use demand triples (Section IV), all others demand
// quadruples over marked paths (Section VI).
class DemandBoundFunction {
public:
    DemandBoundFunction() = default;
    DemandBoundFunction(Time horizon, std::vector<std::pair<Time, Time>> steps)
        : horizon_(horizon), steps_(std::move(steps)) {}

    Time operator()(Time t) const;
    Time horizon() const { return horizon_; }
    // (d, e): dbf_T steps up to e at t = d. Both coordinates strictly increase.
    const std::vector<std::pair<Time, Time>>& steps() const { return steps_; }

private:
    Time                               horizon_ = 0;
    std::vector<std::pair<Time, Time>> steps_;
};

DemandBoundFunction demand_bound_function(const Task& task, Time horizon);

struct FeasibilityResult {
    enum class Status { Feasible, Infeasible, Inconclusive };
    Status              status      = Status::Inconclusive;
    double              utilization = 0.0;  // U(tau)
    Time                bound       = 0;    // every t <= bound was checked
    std::optional<Time> witness;            // smallest t_f with dbf(t_f) > t_f
    std::string         reason;
};

// Preemptive uniprocessor feasibility, equivalently EDF schedulability
// (Proposition III.4, Lemma V.3, Theorems V.4 and VI.6).
FeasibilityResult check_feasibility(const std::vector<Task>& tasks);

const char* status_name(FeasibilityResult::Status s);

// Writes the task graphs, their utilizations, dbf_T step functions (up to
// max(dbf_horizon, bound)) and the feasibility result as JSON, e.g. for
// draw_drt.py.
void write_json(std::ostream& os, const std::vector<Task>& tasks, Time dbf_horizon);

// ======================================================== implementation


// ---------------------------------------------------------------- Task

inline std::size_t Task::add_vertex(std::string id, Time e, Time d) {
    if (e <= 0 || d <= 0) throw std::invalid_argument("e(v) and d(v) must be positive");
    vertices_.push_back({std::move(id), e, d, {}});
    return vertices_.size() - 1;
}

inline void Task::add_edge(std::size_t from, std::size_t to, Time p) {
    if (from >= vertices_.size() || to >= vertices_.size())
        throw std::out_of_range("add_edge: no such vertex");
    if (p <= 0) throw std::invalid_argument("p(u, v) must be positive");
    vertices_[from].out.push_back({to, p});
}

inline bool Task::frame_separated() const {
    for (const auto& u : vertices_)
        for (const auto& edge : u.out)
            if (u.d > edge.p) return false;
    return true;
}

inline Time Task::e_sum() const {
    Time sum = 0;
    for (const auto& v : vertices_) sum += v.e;
    return sum;
}

// ---------------------------------------------------------------- U(T), Section V-A

namespace detail {

struct EP {
    Time e, p;
};

// Keep only non-dominated (e, p): drop x if some y has y.e >= x.e and y.p <= x.p.
inline void prune(std::vector<EP>& v) {
    std::sort(v.begin(), v.end(), [](const EP& a, const EP& b) {
        return a.p < b.p || (a.p == b.p && a.e > b.e);
    });
    std::vector<EP> kept;
    Time best = -1;
    for (const auto& x : v)
        if (x.e > best) {
            kept.push_back(x);
            best = x.e;
        }
    v.swap(kept);
}

}  // namespace detail

inline Utilization utilization(const Task& task) {
    const auto&       V = task.vertices();
    const std::size_t n = V.size();
    Utilization       best;  // 0/1

    // Lemma V.5: some simple cycle attains U(T). Rotating a simple cycle so it
    // starts at its smallest vertex index s, it only visits vertices >= s and
    // has length <= n - s.
    for (std::size_t s = 0; s < n; ++s) {
        // cur[v]: non-dominated utilization triples <e, p, (s, v)> of the current length
        std::vector<std::vector<detail::EP>> cur(n), next(n);
        cur[s].push_back({0, 0});
        for (std::size_t k = 1; k <= n - s; ++k) {
            for (auto& x : next) x.clear();
            for (std::size_t u = s; u < n; ++u)
                for (const auto& ep : cur[u])
                    for (const auto& edge : V[u].out)
                        if (edge.to >= s)
                            next[edge.to].push_back({ep.e + V[u].e, ep.p + edge.p});
            // Domination among paths of the same length, start and end is safe:
            // any cycle completing the dominated one is no denser.
            for (std::size_t v = s; v < n; ++v) detail::prune(next[v]);
            for (const auto& ep : next[s])  // start == end: a cycle, U = e / p
                if (ep.e * best.p > best.e * ep.p) best = {ep.e, ep.p};
            std::swap(cur, next);
        }
    }
    const Time g = std::gcd(best.e, best.p);
    best.e /= g;
    best.p /= g;
    return best;
}

// ---------------------------------------------------------------- dbf_T

namespace detail {

using Pairs = std::vector<std::pair<Time, Time>>;  // demand pairs as (d, e)

// Frame-separated tasks (Section IV): demand triples <e, d, v>, keeping only
// non-dominated ones (Definition IV.7, justified by Lemma IV.9), explored in
// increasing d. Extending always increases d, so once a triple is popped no
// triple that dominates it can still appear.
inline Pairs demand_pairs_triples(const Task& task, Time H) {
    const auto& V = task.vertices();
    // frontier[v]: d -> e with e strictly increasing in d
    std::vector<std::map<Time, Time>> frontier(V.size());
    using Item = std::tuple<Time, Time, std::size_t>;  // (d, e, v)
    std::priority_queue<Item, std::vector<Item>, std::greater<>> queue;

    auto try_add = [&](Time e, Time d, std::size_t v) {
        if (d > H) return;
        auto& f  = frontier[v];
        auto  it = f.upper_bound(d);
        if (it != f.begin() && std::prev(it)->second >= e) return;   // dominated
        it = f.lower_bound(d);
        while (it != f.end() && it->second <= e) it = f.erase(it);  // now dominated
        f[d] = e;
        queue.emplace(d, e, v);
    };

    for (std::size_t v = 0; v < V.size(); ++v) try_add(V[v].e, V[v].d, v);  // 0-paths

    while (!queue.empty()) {
        const auto [d, e, u] = queue.top();
        queue.pop();
        const auto it = frontier[u].find(d);
        if (it == frontier[u].end() || it->second != e) continue;  // dominated meanwhile
        for (const auto& edge : V[u].out)  // Figure 5, lines 6-7
            try_add(e + V[edge.to].e, d - V[u].d + edge.p + V[edge.to].d, edge.to);
    }

    Pairs pairs;
    for (const auto& f : frontier)
        for (const auto& [d, e] : f) pairs.emplace_back(d, e);
    return pairs;
}

// Arbitrary deadlines (Section VI): demand quadruples <e, d, p, v> of marked
// paths. <e1,d1,p1,v> dominates <e2,d2,p2,v> if e1 >= e2, d1 <= d2, p1 <= p2;
// both extension rules preserve this, so dominated quadruples can be dropped.
inline Pairs demand_pairs_quadruples(const Task& task, Time H) {
    struct Quad {
        Time e, d, p;
        bool alive;
    };
    const auto& V     = task.vertices();
    Time        d_min = std::numeric_limits<Time>::max();
    for (const auto& v : V) d_min = std::min(d_min, v.d);

    std::vector<std::vector<Quad>> store(V.size());
    using Item = std::tuple<Time, Time, std::size_t, std::size_t>;  // (p, d, v, index)
    std::priority_queue<Item, std::vector<Item>, std::greater<>> queue;

    auto try_add = [&](Time e, Time d, Time p, std::size_t v) {
        auto& s = store[v];
        for (const auto& q : s)
            if (q.alive && q.e >= e && q.d <= d && q.p <= p) return;
        for (auto& q : s)
            if (q.alive && e >= q.e && d <= q.d && p <= q.p) q.alive = false;
        s.push_back({e, d, p, true});
        queue.emplace(p, d, v, s.size() - 1);
    };

    for (std::size_t v = 0; v < V.size(); ++v)  // the initial vertex is always marked
        if (V[v].d <= H) try_add(V[v].e, V[v].d, 0, v);

    while (!queue.empty()) {
        const auto [p, d, u, i] = queue.top();
        queue.pop();
        const Quad q = store[u][i];  // copy: try_add may reallocate store[u]
        if (!q.alive) continue;
        for (const auto& edge : V[u].out) {
            const std::size_t w  = edge.to;
            const Time        pw = q.p + edge.p;  // release of w relative to the first job
            // w marked: its job counts towards the demand pair
            const Time dw = std::max(q.d, pw + V[w].d);
            if (dw <= H) try_add(q.e + V[w].e, dw, pw, w);
            // w unmarked: only useful if a later marked job could still fit in H
            if (pw + 1 + d_min <= H) try_add(q.e, q.d, pw, w);
        }
    }

    Pairs pairs;
    for (const auto& s : store)
        for (const auto& q : s)
            if (q.alive) pairs.emplace_back(q.d, q.e);
    return pairs;
}

// dbf_T(t) = max { e | <e, d> demand pair, d <= t }: keep the upper staircase.
inline std::vector<std::pair<Time, Time>> to_steps(Pairs pairs) {
    std::sort(pairs.begin(), pairs.end());
    std::vector<std::pair<Time, Time>> steps;
    Time best = 0;
    for (const auto& [d, e] : pairs) {
        if (e <= best) continue;
        if (!steps.empty() && steps.back().first == d)
            steps.back().second = e;
        else
            steps.emplace_back(d, e);
        best = e;
    }
    return steps;
}

}  // namespace detail

inline DemandBoundFunction demand_bound_function(const Task& task, Time horizon) {
    return {horizon, detail::to_steps(task.frame_separated()
                                          ? detail::demand_pairs_triples(task, horizon)
                                          : detail::demand_pairs_quadruples(task, horizon))};
}

inline Time DemandBoundFunction::operator()(Time t) const {
    if (t > horizon_) throw std::out_of_range("dbf queried beyond its horizon");
    const auto it = std::upper_bound(steps_.begin(), steps_.end(),
                                     std::make_pair(t, std::numeric_limits<Time>::max()));
    return it == steps_.begin() ? 0 : std::prev(it)->second;
}

// ---------------------------------------------------------------- feasibility

namespace detail {
using i128 = __int128;  // GCC/Clang extension; keeps U(tau) exact

inline i128 gcd128(i128 a, i128 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) {
        const i128 r = a % b;
        a = b;
        b = r;
    }
    return a;
}
}  // namespace detail

inline FeasibilityResult check_feasibility(const std::vector<Task>& tasks) {
    using Status = FeasibilityResult::Status;
    FeasibilityResult r;
    if (tasks.empty()) {
        r.status = Status::Feasible;
        r.reason = "empty task set";
        return r;
    }

    // U(tau) = num / den exactly
    using detail::i128;
    i128 num = 0, den = 1;
    Time e_total = 0;
    for (const auto& task : tasks) {
        const Utilization u = utilization(task);
        num                 = num * u.p + i128(u.e) * den;
        den *= u.p;
        const i128 g = detail::gcd128(num, den);
        num /= g;
        den /= g;
        e_total += task.e_sum();
    }
    r.utilization = static_cast<double>(num) / static_cast<double>(den);

    if (num > den) {
        r.status = Status::Infeasible;
        r.reason = "U(tau) > 1: the most dense cycles overload the processor";
        return r;
    }
    if (num == den) {
        r.status = Status::Inconclusive;
        r.reason = "U(tau) = 1: Lemma V.3 gives no finite bound for t_f";
        return r;
    }

    // Lemma V.3: a witness satisfies t_f < e_sum / (1 - U), i.e. t_f <= bound
    const i128 bound = (i128(e_total) * den - 1) / (den - num);
    if (bound > i128(std::numeric_limits<Time>::max() / 4)) {
        r.status = Status::Inconclusive;
        r.reason = "bound for t_f too large (U(tau) too close to 1)";
        return r;
    }
    r.bound = static_cast<Time>(bound);

    // dbf(t) only changes at step points of the individual dbf_T: sweep them in order.
    std::vector<std::pair<Time, Time>> events;  // (t, increase of dbf at t)
    for (const auto& task : tasks) {
        const auto f    = demand_bound_function(task, r.bound);
        Time       prev = 0;
        for (const auto& [d, e] : f.steps()) {
            events.emplace_back(d, e - prev);
            prev = e;
        }
    }
    std::sort(events.begin(), events.end());

    Time demand = 0;
    for (std::size_t i = 0; i < events.size();) {
        const Time t = events[i].first;
        for (; i < events.size() && events[i].first == t; ++i) demand += events[i].second;
        if (demand > t) {
            r.status  = Status::Infeasible;
            r.witness = t;
            r.reason  = "dbf(t) > t (Proposition III.4)";
            return r;
        }
    }
    r.status = Status::Feasible;
    r.reason = "dbf(t) <= t for all t <= bound (Proposition III.4, Lemma V.3)";
    return r;
}

// ---------------------------------------------------------------- JSON export

inline const char* status_name(FeasibilityResult::Status s) {
    switch (s) {
        case FeasibilityResult::Status::Feasible:   return "FEASIBLE";
        case FeasibilityResult::Status::Infeasible: return "INFEASIBLE";
        default:                                    return "INCONCLUSIVE";
    }
}

namespace detail {
inline void json_string(std::ostream& os, const std::string& s) {
    os << '"';
    for (const char c : s) {
        switch (c) {
            case '"':  os << "\\\""; break;
            case '\\': os << "\\\\"; break;
            case '\n': os << "\\n"; break;
            case '\t': os << "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof buf, "\\u%04x", static_cast<unsigned>(c));
                    os << buf;
                } else {
                    os << c;
                }
        }
    }
    os << '"';
}
}  // namespace detail

inline void write_json(std::ostream& os, const std::vector<Task>& tasks, Time dbf_horizon) {
    const FeasibilityResult r       = check_feasibility(tasks);
    const Time              horizon = std::max(dbf_horizon, r.bound);

    os << "{\n  \"tasks\": [";
    for (std::size_t i = 0; i < tasks.size(); ++i) {
        const Task&       t = tasks[i];
        const Utilization u = utilization(t);
        os << (i ? "," : "") << "\n    {\n      \"name\": ";
        detail::json_string(os, t.name());
        os << ",\n      \"frame_separated\": " << (t.frame_separated() ? "true" : "false")
           << ",\n      \"utilization\": {\"e\": " << u.e << ", \"p\": " << u.p << "}"
           << ",\n      \"vertices\": [";
        const auto& V = t.vertices();
        for (std::size_t v = 0; v < V.size(); ++v) {
            os << (v ? "," : "") << "\n        {\"id\": ";
            detail::json_string(os, V[v].id);
            os << ", \"e\": " << V[v].e << ", \"d\": " << V[v].d << "}";
        }
        os << "\n      ],\n      \"edges\": [";
        bool first = true;
        for (std::size_t v = 0; v < V.size(); ++v)
            for (const auto& edge : V[v].out) {
                os << (first ? "" : ",") << "\n        {\"from\": " << v << ", \"to\": " << edge.to
                   << ", \"p\": " << edge.p << "}";
                first = false;
            }
        os << "\n      ],\n      \"dbf\": {\"horizon\": " << horizon << ", \"steps\": [";
        const auto f = demand_bound_function(t, horizon);
        for (std::size_t k = 0; k < f.steps().size(); ++k)
            os << (k ? ", " : "") << "[" << f.steps()[k].first << ", " << f.steps()[k].second << "]";
        os << "]}\n    }";
    }
    os << "\n  ],\n  \"feasibility\": {\n    \"status\": \"" << status_name(r.status) << "\""
       << ",\n    \"utilization\": " << r.utilization << ",\n    \"bound\": " << r.bound
       << ",\n    \"witness\": ";
    if (r.witness) os << *r.witness; else os << "null";
    os << ",\n    \"reason\": ";
    detail::json_string(os, r.reason);
    os << "\n  }\n}\n";
}

}  // namespace drt

#endif
