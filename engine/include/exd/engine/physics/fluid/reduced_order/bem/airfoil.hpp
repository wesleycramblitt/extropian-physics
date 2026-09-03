#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace exd::engine::physics::fluid::reduced_order::bem
{

/// Single Re-tagged polar table for an airfoil.
struct AirfoilPolar
{
    std::string name;                 // e.g. "naca0012"
    double re = 0.0;                  // design Re; 0 = Re-independent fallback
    std::vector<double> alpha_deg;    // strictly increasing
    std::vector<double> cl;
    std::vector<double> cd;

    [[nodiscard]] bool valid() const
    {
        return !name.empty() && alpha_deg.size() >= 2 &&
               alpha_deg.size() == cl.size() && cl.size() == cd.size();
    }

    /// Linear interpolation in alpha; flat clamp beyond table ends.
    [[nodiscard]] std::pair<double, double> evaluate(double alpha) const
    {
        if (alpha_deg.empty()) return {0.0, 0.0};
        if (alpha <= alpha_deg.front()) return {cl.front(), cd.front()};
        if (alpha >= alpha_deg.back())  return {cl.back(), cd.back()};

        auto it = std::lower_bound(alpha_deg.begin(), alpha_deg.end(), alpha);
        std::size_t i = static_cast<std::size_t>(it - alpha_deg.begin());
        if (i == 0) i = 1;
        const std::size_t im1 = i - 1;
        const double a0 = alpha_deg[im1];
        const double a1 = alpha_deg[i];
        const double t = (alpha - a0) / (a1 - a0 + std::numeric_limits<double>::epsilon());
        const double cl_val = cl[im1] + t * (cl[i] - cl[im1]);
        const double cd_val = cd[im1] + t * (cd[i] - cd[im1]);
        return {cl_val, cd_val};
    }
};

/// Database of airfoil polars. Supports nearest-Re selection and CSV loading.
class PolarDatabase
{
public:
    void add(AirfoilPolar p)
    {
        auto key = std::make_pair(p.name, p.re);
        storage_[key] = std::move(p);
    }

    void add_builtin_polars();

    bool load_csv(const std::string& path);
    bool load_directory(const std::string& path);

    /// Nearest-Re selection. A polar with re == 0 has distance +inf and is only
    /// used when no Re-tagged polar exists for that name.
    const AirfoilPolar* find(const std::string& name, double re) const
    {
        const AirfoilPolar* best = nullptr;
        const AirfoilPolar* zero = nullptr;
        double best_dist = std::numeric_limits<double>::infinity();

        for (const auto& [key, polar] : storage_)
        {
            if (key.first != name) continue;
            if (polar.re == 0.0)
            {
                zero = &polar;
                continue;
            }
            const double dist = std::fabs(polar.re - re);
            if (dist < best_dist)
            {
                best_dist = dist;
                best = &polar;
            }
        }
        return best ? best : zero;
    }

    bool has(const std::string& name) const
    {
        for (const auto& [key, polar] : storage_)
            if (key.first == name) return true;
        return false;
    }

    std::vector<std::string> airfoil_names() const
    {
        std::vector<std::string> names;
        for (const auto& [key, polar] : storage_)
        {
            if (std::find(names.begin(), names.end(), key.first) == names.end())
                names.push_back(key.first);
        }
        std::sort(names.begin(), names.end());
        return names;
    }

    bool empty() const noexcept { return storage_.empty(); }
    explicit operator bool() const { return !empty(); }

private:
    using Key = std::pair<std::string, double>;
    struct KeyHash
    {
        std::size_t operator()(const Key& k) const noexcept
        {
            std::size_t h1 = std::hash<std::string>{}(k.first);
            std::size_t h2 = std::hash<double>{}(k.second);
            return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
        }
    };
    std::unordered_map<Key, AirfoilPolar, KeyHash> storage_;

    void add_builtin(std::string name, double re,
                     std::vector<double> alpha_deg,
                     std::vector<double> cl,
                     std::vector<double> cd)
    {
        AirfoilPolar p;
        p.name = std::move(name);
        p.re = re;
        p.alpha_deg = std::move(alpha_deg);
        p.cl = std::move(cl);
        p.cd = std::move(cd);
        add(std::move(p));
    }
};

} // namespace exd::engine::physics::fluid::reduced_order::bem
