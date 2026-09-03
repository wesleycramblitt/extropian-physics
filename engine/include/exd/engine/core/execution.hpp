#pragma once

// ─────────────────────────────────────────────────────
// Execution graph (implementation_spec §38, §54, §58).
//
// The execution graph represents computational
// dependencies (ComputeVelocity → ComputeFlux →
// ComputeResidual → SolvePressure → UpdateVelocity …).
// Nodes are steps over state/domains; edges encode
// producer→consumer dependencies; execution follows a
// topological order.  The execution layer owns
// scheduling, synchronization, placement (CPU/GPU),
// memory lifetime and — with a GPU backend — streams and
// events.
//
// Domain steps, coupling exchanges, diagnostics and
// output sinks all register as graph nodes, so a coupled
// run is one DAG (CouplingManager + CoupledSimulation
// remain the reference staggered/implicit drivers and
// surface here as graph composition helpers).
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace exd::engine::core {

enum class NodeKind : uint8_t
{
    Compute,     // domain advance / operator evaluation
    Exchange,    // coupling exchange
    Output,      // persistent output
    Diagnostics, // diagnostic evaluation
    External,    // external consumer/plugin step
};

constexpr const char* to_string(NodeKind k)
{
    switch (k)
    {
    case NodeKind::Compute: return "compute";
    case NodeKind::Exchange: return "exchange";
    case NodeKind::Output: return "output";
    case NodeKind::Diagnostics: return "diagnostics";
    case NodeKind::External: return "external";
    }
    return "?";
}

/// One execution-graph node: a step callback plus its dependency set.
/// The callback returns ok/false through `status` (no exceptions).
struct GraphNode
{
    std::string id;
    NodeKind kind = NodeKind::Compute;
    std::function<bool(double t, double dt, ModelStatus&)> step;  // returns ok
    std::vector<std::string> depends_on;                          // node ids
};

/// Dependency-ordered execution DAG (§38).  Validated (cycle detection,
/// missing dependency rejection) before first execute; serial CPU
/// execution today, async/stream scheduling with a GPU backend.
class ExecutionGraph
{
public:
    bool add_node(GraphNode node, ModelStatus& status);
    const GraphNode* node(std::string_view id) const;

    /// Topologically sort and validate; call once before execute.
    bool finalize(ModelStatus& status);

    /// Execute all nodes in topological order at time t.
    bool execute(double t, double dt, ModelStatus& status);

    size_t node_count() const { return nodes_.size(); }
    const std::vector<std::string>& order() const { return order_; }

    /// Remove the node (by id); invalidates order until finalize again.
    bool remove_node(std::string_view id, ModelStatus& status);

private:
    bool topo_sort(ModelStatus& status);
    std::vector<GraphNode> nodes_;
    std::vector<std::string> order_;
    bool finalized_ = false;
};

inline bool ExecutionGraph::add_node(GraphNode node, ModelStatus& status)
{
    for (auto& n : nodes_)
    {
        if (n.id == node.id)
        {
            status.ok = false;
            status.error = "execution graph: duplicate node '" + node.id + "'";
            return false;
        }
    }
    if (node.id.empty())
    {
        status.ok = false;
        status.error = "execution graph: node id must not be empty";
        return false;
    }
    nodes_.push_back(std::move(node));
    finalized_ = false;
    return true;
}

inline const GraphNode* ExecutionGraph::node(std::string_view id) const
{
    for (auto& n : nodes_)
        if (n.id == id)
            return &n;
    return nullptr;
}

inline bool ExecutionGraph::remove_node(std::string_view id, ModelStatus& status)
{
    for (size_t i = 0; i < nodes_.size(); ++i)
    {
        if (nodes_[i].id == id)
        {
            nodes_.erase(nodes_.begin() + static_cast<std::ptrdiff_t>(i));
            finalized_ = false;
            return true;
        }
    }
    status.ok = false;
    status.error = "execution graph: no node '" + std::string(id) + "'";
    return false;
}

inline bool ExecutionGraph::topo_sort(ModelStatus& status)
{
    // Kahn's algorithm with deterministic tie-breaking by insertion order.
    const size_t n = nodes_.size();
    std::vector<size_t> indeg(n, 0);
    // map id -> index
    std::vector<std::pair<std::string, size_t>> idmap;
    idmap.reserve(n);
    for (size_t i = 0; i < n; ++i) idmap.emplace_back(nodes_[i].id, i);
    auto index_of = [&](std::string_view id) -> int {
        for (auto& [name, idx] : idmap)
            if (name == id)
                return static_cast<int>(idx);
        return -1;
    };
    for (size_t j = 0; j < n; ++j)
    {
        for (auto& dep : nodes_[j].depends_on)
        {
            int di = index_of(dep);
            if (di < 0)
            {
                status.ok = false;
                status.error = "execution graph: node '" + nodes_[j].id +
                               "' depends on unknown node '" + dep + "'";
                return false;
            }
            ++indeg[j];   // incoming edge (dep → j): count for the DEPENDENT node
        }
    }
    std::vector<size_t> ready;
    for (size_t i = 0; i < n; ++i)
        if (indeg[i] == 0)
            ready.push_back(i);
    order_.clear();
    while (!ready.empty())
    {
        const size_t cur = ready.front();
        ready.erase(ready.begin());
        order_.push_back(nodes_[cur].id);
        for (size_t j = 0; j < n; ++j)
        {
            for (auto& dep : nodes_[j].depends_on)
            {
                if (dep == nodes_[cur].id)
                {
                    --indeg[j];
                    if (indeg[j] == 0)
                    {
                        // deterministic insertion: keep ready sorted by index
                        size_t pos = 0;
                        while (pos < ready.size() && ready[pos] < j) ++pos;
                        ready.insert(ready.begin() + static_cast<std::ptrdiff_t>(pos), j);
                    }
                }
            }
        }
    }
    if (order_.size() != n)
    {
        status.ok = false;
        status.error = "execution graph: cycle detected (" +
                       std::to_string(n - order_.size()) + " nodes unreachable)";
        return false;
    }
    return true;
}

inline bool ExecutionGraph::finalize(ModelStatus& status)
{
    if (!topo_sort(status)) return false;
    finalized_ = true;
    return true;
}

inline bool ExecutionGraph::execute(double t, double dt, ModelStatus& status)
{
    if (!finalized_ && !finalize(status)) return false;
    for (auto& id : order_)
    {
        GraphNode* node = nullptr;
        for (auto& n : nodes_)
            if (n.id == id)
            {
                node = &n;
                break;
            }
        if (!node)
        {
            status.ok = false;
            status.error = "execution graph: internal inconsistency at '" + id + "'";
            return false;
        }
        if (!node->step || !node->step(t, dt, status))
        {
            status.ok = false;
            if (status.error.empty())
                status.error = "execution graph: node '" + id + "' failed";
            return false;
        }
    }
    return true;
}

} // namespace exd::engine::core
