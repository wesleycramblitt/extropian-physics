#include "bem_internal.hpp"

#include <exd/geometry/turbine.hpp>

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>
#include <vector>

namespace exd::physics::fluid::reduced_order::bem
{

namespace
{

using namespace exd::geometry;

struct Vec2 { double x = 0.0, y = 0.0; };

Vec2 operator+(const Vec2& a, const Vec2& b) { return {a.x + b.x, a.y + b.y}; }
Vec2 operator-(const Vec2& a, const Vec2& b) { return {a.x - b.x, a.y - b.y}; }
Vec2 operator*(const Vec2& a, double s) { return {a.x * s, a.y * s}; }

double length(const Vec2& v) { return std::sqrt(v.x * v.x + v.y * v.y); }

Vec2 to_vec(math::Vec2f v) { return {static_cast<double>(v.x), static_cast<double>(v.y)}; }

Vec2 lerp(const Vec2& a, const Vec2& b, double f) { return a + (b - a) * f; }

MonotoneCubicSpline make_spline(const std::vector<math::Vec2f>& points)
{
    std::vector<float> xs, ys;
    xs.reserve(points.size());
    ys.reserve(points.size());
    for (const auto& p : points)
    {
        xs.push_back(p.x);
        ys.push_back(p.y);
    }
    return MonotoneCubicSpline(std::move(xs), std::move(ys));
}

struct ResolvedRow
{
    const BladeRow* row = nullptr;
    uint32_t index = 0;
    double rpm = 0.0;
};

std::optional<ResolvedRow> resolve_rotor_row(const TurbineDefinition& turbine,
                                             const BEMSolverConfig& config,
                                             std::vector<std::string>& warnings)
{
    ResolvedRow result;
    int rotor_count = 0;
    for (std::size_t i = 0; i < turbine.blade_rows.size(); ++i)
    {
        if (turbine.blade_rows[i].type == BladeRowType::Rotor)
        {
            ++rotor_count;
            result.row = &turbine.blade_rows[i];
            result.index = static_cast<uint32_t>(i);
        }
        else
        {
            warnings.push_back("non-rotor blade row ignored");
        }
    }

    if (rotor_count == 0)
        return std::nullopt;

    if (rotor_count > 1)
        return std::nullopt;

    if (config.row_index != result.index)
    {
        warnings.push_back("config.row_index does not match the single Rotor row; row_index ignored");
    }

    result.rpm = static_cast<double>(result.row->rotational_speed.value);
    return result;
}

double stagger_at_span(const BladeRow& row, double span, std::vector<std::string>& warnings)
{
    const auto& sections = row.sections;
    if (sections.empty())
    {
        warnings.push_back("BladeRow.sections empty; zero-stagger defaults");
        return 0.0;
    }

    std::vector<std::pair<double, double>> pts;
    pts.reserve(sections.size());
    for (const auto& s : sections)
        pts.emplace_back(static_cast<double>(s.span), static_cast<double>(s.stagger.value));
    std::sort(pts.begin(), pts.end(), [](const auto& a, const auto& b){ return a.first < b.first; });

    if (span <= pts.front().first)
    {
        if (span < pts.front().first)
            warnings.push_back("stagger extrapolated at low span");
        return pts.front().second;
    }
    if (span >= pts.back().first)
    {
        if (span > pts.back().first)
            warnings.push_back("stagger extrapolated at high span");
        return pts.back().second;
    }

    auto it = std::lower_bound(pts.begin(), pts.end(), span,
        [](const auto& a, double v){ return a.first < v; });
    const std::size_t i = static_cast<std::size_t>(it - pts.begin());
    const std::size_t im1 = i - 1;
    const double f = (span - pts[im1].first) / (pts[i].first - pts[im1].first);
    return pts[im1].second + f * (pts[i].second - pts[im1].second);
}

std::string airfoil_at_span(const std::vector<AirfoilAssignment>& assignments, double span)
{
    if (assignments.empty()) return "naca0012";

    const AirfoilAssignment* best = nullptr;
    for (const auto& a : assignments)
    {
        if (a.span <= span && (!best || a.span > best->span))
            best = &a;
    }
    if (!best) best = &assignments.front();
    return best->airfoil;
}

} // namespace

GeometryResult build_blade_geometry(const TurbineDefinition& turbine,
                                    const BEMSolverConfig& config)
{
    GeometryResult out;
    auto& warnings = out.warnings;

    if (turbine.flow_path.shroud_points.empty() || turbine.flow_path.hub_points.empty())
    {
        out.error = "empty flow path";
        return out;
    }

    auto maybe_row = resolve_rotor_row(turbine, config, warnings);
    if (!maybe_row)
    {
        int rotor_count = 0;
        for (const auto& row : turbine.blade_rows)
            if (row.type == BladeRowType::Rotor) ++rotor_count;
        if (rotor_count == 0)
            out.error = "no Rotor row in TurbineDefinition";
        else
            out.error = std::to_string(rotor_count) + " Rotor rows; single-rotor solver";
        return out;
    }

    const auto& row = *maybe_row;
    if (row.rpm <= 0.0)
    {
        out.error = "rpm <= 0";
        return out;
    }

    const Vec2 le_hub = to_vec(row.row->leading_edge_hub);
    const Vec2 le_shroud = to_vec(row.row->leading_edge_shroud);
    const Vec2 te_hub = to_vec(row.row->trailing_edge_hub);
    const Vec2 te_shroud = to_vec(row.row->trailing_edge_shroud);

    const double z_r = 0.5 * (0.5 * (le_hub.x + le_shroud.x) + 0.5 * (te_hub.x + te_shroud.x));

    const double r_hub = 0.5 * (le_hub.y + te_hub.y);
    double r_tip = 0.5 * (le_shroud.y + te_shroud.y);

    auto shroud_spline = make_spline(turbine.flow_path.shroud_points);
    auto hub_spline = make_spline(turbine.flow_path.hub_points);

    if (row.row->tip_feature == TipFeature::Clearance)
    {
        const double r_shroud_at_zr = static_cast<double>(shroud_spline.evaluate(static_cast<float>(z_r)));
        const double clearance = static_cast<double>(turbine.flow_path.tip_clearance.value);
        r_tip = std::min(r_tip, r_shroud_at_zr - clearance);
    }

    if (r_tip <= r_hub)
    {
        out.error = "R_tip <= R_hub";
        return out;
    }

    const Vec2 le_mid = lerp(le_hub, le_shroud, 0.5);
    const Vec2 te_mid = lerp(te_hub, te_shroud, 0.5);
    if (length(te_mid - le_mid) <= 0.0)
    {
        out.error = "chord <= 0";
        return out;
    }

    const uint32_t n = config.element_count;
    const double dr = (r_tip - r_hub) / static_cast<double>(n);

    BladeGeometry geo;
    geo.z_r = z_r;
    geo.r_hub = r_hub;
    geo.r_tip = r_tip;
    geo.blade_count = static_cast<double>(row.row->blade_count.value);
    geo.rpm = row.rpm;
    geo.omega = row.rpm * 2.0 * M_PI / 60.0;
    geo.shroud_spline = std::move(shroud_spline);
    geo.hub_spline = std::move(hub_spline);
    geo.elements.reserve(n);

    for (uint32_t i = 0; i < n; ++i)
    {
        const double r = r_hub + (static_cast<double>(i) + 0.5) * dr;
        const double span = (r - r_hub) / (r_tip - r_hub);

        const Vec2 le = lerp(le_hub, le_shroud, span);
        const Vec2 te = lerp(te_hub, te_shroud, span);
        const double chord = length(te - le);
        if (chord <= 0.0)
        {
            out.error = "chord <= 0 at span " + std::to_string(span);
            return out;
        }

        const double beta = stagger_at_span(*row.row, span, warnings);
        BladeElementInput e;
        e.r = r;
        e.dr = dr;
        e.chord = chord;
        e.beta_deg = beta;
        e.span = span;
        e.blade_count = geo.blade_count;
        e.r_hub = r_hub;
        e.r_tip = r_tip;
        e.airfoil = airfoil_at_span(config.airfoils, span);
        geo.elements.push_back(e);
    }

    out.ok = true;
    out.geometry = std::move(geo);
    return out;
}

} // namespace exd::physics::fluid::reduced_order::bem
