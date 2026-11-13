#!/usr/bin/env python3
"""
Build the standalone, zero-server `modelling_dashboard/` showcase site.

Reads real CSVs, Markdown evaluation reports, and rendered PNG/SVG figures from
`output/`, copies curated high-DPI PNG maps (`--dpi=300`) and Python-minified
vector SVGs into `modelling_dashboard/assets/figures/`, and writes a synchronous
`modelling_dashboard/data.js` bundle safe for direct `file://` browser viewing.
"""

import csv
import json
import os
import re
import shutil

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(BASE_DIR, "output")
DASH_DIR = os.path.join(BASE_DIR, "modelling_dashboard")
FIG_DEST_DIR = os.path.join(DASH_DIR, "assets", "figures")

SVG_MAX_BYTES = 1_200_000  # Keep SVG for pure statistical charts; use 300 DPI PNG for dense 90k maps/scatters


def coerce_value(val):
    if val is None:
        return None
    s = val.strip()
    if s == "" or s.upper() == "NA":
        return None
    try:
        if re.fullmatch(r"-?\d+", s):
            return int(s)
        return round(float(s), 6)
    except ValueError:
        return s


def read_csv_rows(rel_path):
    full_path = os.path.join(BASE_DIR, rel_path)
    if not os.path.exists(full_path):
        print(f"  [WARN] Missing CSV: {rel_path}")
        return []
    with open(full_path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        return [{k: coerce_value(v) for k, v in row.items()} for row in reader]


def read_text_excerpt(rel_path, max_lines=160):
    full_path = os.path.join(BASE_DIR, rel_path)
    if not os.path.exists(full_path):
        return ""
    with open(full_path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()[:max_lines]
    return "\n".join(lines)


def _round_float_match(m):
    num_str = m.group(1)
    rounded = f"{float(num_str):.2f}".rstrip("0").rstrip(".")
    return rounded if rounded != "-0" else "0"


def optimize_svg(svg_text):
    """Minify SVG safely: strip comments, collapse tag gaps, round attribute floats only."""
    cleaned = re.sub(r"<!--.*?-->", "", svg_text, flags=re.DOTALL)
    # Round long floats strictly inside tag markup (<...>), preserving <text> content verbatim
    parts = re.split(r"(<[^>]+>)", cleaned)
    out_parts = []
    for part in parts:
        if part.startswith("<") and part.endswith(">"):
            part = re.sub(r"(?<![A-Za-z0-9_#.-])(-?\d+\.\d{3,})", _round_float_match, part)
            part = re.sub(r"\s{2,}", " ", part)
        out_parts.append(part)
    minified = "".join(out_parts)
    minified = re.sub(r">\s+<", "><", minified).strip()
    return minified


def stage_figure(rel_base_no_ext, alias_name, allow_svg=True):
    """Copy either minified SVG (for statistical charts) or 300-DPI PNG (for dense maps)."""
    os.makedirs(FIG_DEST_DIR, exist_ok=True)
    svg_src = os.path.join(BASE_DIR, rel_base_no_ext + ".svg")
    png_src = os.path.join(BASE_DIR, rel_base_no_ext + ".png")

    use_svg = False
    if allow_svg and os.path.exists(svg_src) and os.path.getsize(svg_src) <= SVG_MAX_BYTES:
        use_svg = True
    elif not os.path.exists(png_src) and os.path.exists(svg_src):
        use_svg = True

    if use_svg:
        dest_file = alias_name + ".svg"
        dest_path = os.path.join(FIG_DEST_DIR, dest_file)
        with open(svg_src, "r", encoding="utf-8") as f:
            raw_svg = f.read()
        opt_svg = optimize_svg(raw_svg)
        with open(dest_path, "w", encoding="utf-8") as f:
            f.write(opt_svg)
        raw_kb = os.path.getsize(svg_src) / 1024.0
        opt_kb = os.path.getsize(dest_path) / 1024.0
        print(f"  [SVG] {dest_file}: {raw_kb:.0f} KB -> {opt_kb:.0f} KB (minified)")
        return {
            "path": f"assets/figures/{dest_file}",
            "format": "SVG (Vector Minified)",
            "size_kb": round(opt_kb, 1),
        }

    if os.path.exists(png_src):
        dest_file = alias_name + ".png"
        dest_path = os.path.join(FIG_DEST_DIR, dest_file)
        shutil.copy2(png_src, dest_path)
        size_kb = os.path.getsize(dest_path) / 1024.0
        print(f"  [PNG] {dest_file}: {size_kb:.0f} KB (300 DPI)")
        return {
            "path": f"assets/figures/{dest_file}",
            "format": "PNG (300 DPI)",
            "size_kb": round(size_kb, 1),
        }

    print(f"  [WARN] Missing figure asset: {rel_base_no_ext}.(png|svg)")
    return None


def main():
    print("Building modelling_dashboard/ showcase bundle...")
    os.makedirs(FIG_DEST_DIR, exist_ok=True)

    figures = {
        # Case Study 1: Nitrate (0117) + Non-Linear GAM
        "no3_england": stage_figure("output/results/det_0117_orn_rsn/fig1_england_stream_0117_map", "cs1_no3_england", allow_svg=False),
        "no3_humber": stage_figure("output/results/det_0117_orn_rsn/fig2_humber_basin_0117_map", "cs1_no3_humber", allow_svg=False),
        "no3_thames_severn": stage_figure("output/results/det_0117_orn_rsn/fig3_thames_severn_0117_maps", "cs1_no3_thames_severn", allow_svg=False),
        "no3_diag": stage_figure("output/results/det_0117_orn_rsn/fig4_cv_and_diagnostics_0117", "cs1_no3_diagnostics", allow_svg=True),
        "no3_gam_england": stage_figure("output/results/det_0117_orn_rsn_trend_gam/fig1_england_stream_0117_map", "cs1_no3_gam_england", allow_svg=False),
        "no3_gam_diag": stage_figure("output/results/det_0117_orn_rsn_trend_gam/fig4_cv_and_diagnostics_0117", "cs1_no3_gam_diagnostics", allow_svg=True),

        # Case Study 2: Ecological Quality (WHPT ASPT), Risk & Robust Z-Score Normalization
        "whpt_england": stage_figure("output/results/inverts_whpt_aspt_orn_rsn/fig1_england_stream_whpt_aspt_map", "cs2_whpt_england", allow_svg=False),
        "whpt_thames_severn": stage_figure("output/results/inverts_whpt_aspt_orn_rsn/fig3_thames_severn_whpt_aspt_maps", "cs2_whpt_thames_severn", allow_svg=False),
        "whpt_diag": stage_figure("output/results/inverts_whpt_aspt_orn_rsn/fig4_cv_and_diagnostics_whpt_aspt", "cs2_whpt_diagnostics", allow_svg=True),
        "no3_z_england": stage_figure("output/results/det_0117_orn_rsn/fig1_england_stream_0117_map_robust_zscore", "cs2_no3_z_england", allow_svg=False),
        "no3_z_humber": stage_figure("output/results/det_0117_orn_rsn/fig2_humber_basin_0117_map_robust_zscore", "cs2_no3_z_humber", allow_svg=False),
        "no3_z_thames_severn": stage_figure("output/results/det_0117_orn_rsn/fig3_thames_severn_0117_maps_robust_zscore", "cs1_no3_z_thames_severn", allow_svg=False),
        "no3_exceed_england": stage_figure("output/results/det_0117_orn_rsn/fig_exceedance_england_0117", "cs2_no3_exceed_england", allow_svg=False),
        "no3_exceed_basins": stage_figure("output/results/det_0117_orn_rsn/fig_exceedance_basins_0117", "cs2_no3_exceed_basins", allow_svg=False),
        "po4_exceed_england": stage_figure("output/results/det_0180_orn_rsn/fig_exceedance_england_0180", "cs2_po4_exceed_england", allow_svg=False),
        "opt_design_england": stage_figure("output/design_0117/results/det_0117_orn_rsn/fig_optimal_design_england_0117", "cs2_opt_design_england", allow_svg=False),
        "opt_design_humber": stage_figure("output/design_0117/results/det_0117_orn_rsn/fig_optimal_design_humber_0117", "cs2_opt_design_humber", allow_svg=False),
        "opt_design_thames": stage_figure("output/design_0117/results/det_0117_orn_rsn/fig_optimal_design_thames_0117", "cs2_opt_design_thames", allow_svg=False),
        "opt_design_severn": stage_figure("output/design_0117/results/det_0117_orn_rsn/fig_optimal_design_severn_0117", "cs2_opt_design_severn", allow_svg=False),

        # Case Study 3: Spatio-Temporal Trends & Borrowing
        "st_year_grid": stage_figure("output/results/det_0117_orn_rsn_st/fig5_england_by_year_0117_map", "cs3_st_year_grid", allow_svg=False),
        "st_season_grid": stage_figure("output/results/det_0117_orn_rsn_st/fig5_england_by_season_0117_map", "cs3_st_season_grid", allow_svg=False),
        "st_year_trends": stage_figure("output/results/det_0117_orn_rsn_st/fig6_spatiotemporal_trends_by_year_0117", "cs3_st_year_trends", allow_svg=False),
        "st_season_trends": stage_figure("output/results/det_0117_orn_rsn_st/fig6_spatiotemporal_trends_by_season_0117", "cs3_st_season_trends", allow_svg=False),
        "st_season_trans": stage_figure("output/results/det_0117_orn_rsn_st/fig7_spatiotemporal_transitions_by_season_0117", "cs3_st_season_transitions", allow_svg=False),

        # Case Study 4: Multivariate Co-Kriging & Multi-Source RSN+SSN Fusion
        "cokrig_xcorr": stage_figure("output/results/cokrig_0117_0180_orn_rsn/fig_cokrig_cross_correlation", "cs4_cokrig_cross_corr", allow_svg=True),
        "cokrig_0117": stage_figure("output/results/cokrig_0117_0180_orn_rsn/fig_cokrig_england_0117", "cs4_cokrig_0117", allow_svg=False),
        "cokrig_0180": stage_figure("output/results/cokrig_0117_0180_orn_rsn/fig_cokrig_england_0180", "cs4_cokrig_0180", allow_svg=False),
        "rss_england": stage_figure("output/results/det_0117_orn_rss/fig1_england_stream_0117_map", "cs4_rss_england", allow_svg=False),
        "rss_diag": stage_figure("output/results/det_0117_orn_rss/fig4_cv_and_diagnostics_0117", "cs4_rss_diagnostics", allow_svg=True),

        # Case Study 5: Basin Audits, Descriptive Stats, N:P Comparison & Station Snapping
        "desc_dist": stage_figure("output/results/det_0117_orn_rsn/descriptives/fig_desc1_distributions_and_representativeness", "cs5_desc_dist", allow_svg=False),
        "desc_grad": stage_figure("output/results/det_0117_orn_rsn/descriptives/fig_desc2_stratified_gradients", "cs5_desc_grad", allow_svg=True),
        "comp_cv": stage_figure("output/comparisons/det_0117_orn_rsn_vs_det_0180_orn_rsn/fig_comp1_cv_performance", "cs5_comp_cv", allow_svg=True),
        "comp_forest": stage_figure("output/comparisons/det_0117_orn_rsn_vs_det_0180_orn_rsn/fig_comp3_parameter_forest", "cs5_comp_forest", allow_svg=True),
        "comp_maps": stage_figure("output/comparisons/det_0117_orn_rsn_vs_det_0180_orn_rsn/fig_comp4_side_by_side_maps", "cs5_comp_maps", allow_svg=False),
        "comp_rel": stage_figure("output/comparisons/det_0117_orn_rsn_vs_det_0180_orn_rsn/fig_comp5_reach_relationships", "cs5_comp_rel", allow_svg=False),
        "snap_overview": stage_figure("output/station_snapping/orn_rsn/fig_snap_overview_orn", "cs5_snap_overview", allow_svg=False),
        "snap_thames": stage_figure("output/station_snapping/orn_rsn/fig_snap_catchment_3055_river_thames", "cs5_snap_thames", allow_svg=False),
        "snap_humber": stage_figure("output/station_snapping/orn_rsn/fig_snap_catchment_80_river_humber", "cs5_snap_humber", allow_svg=False),

        # Case Study 6: RSN + SSN + Citizen Science WaterBlitz Co-Kriging (Equal vs Favoured RSN+SSN)
        "wb_diff_0117": stage_figure("output/comparisons/cs6_equal_vs_favoured/fig_cs6_diff_maps_0117", "cs6_diff_maps_0117", allow_svg=False),
        "wb_diff_0180": stage_figure("output/comparisons/cs6_equal_vs_favoured/fig_cs6_diff_maps_0180", "cs6_diff_maps_0180", allow_svg=False),
        "wb_comp_svg": stage_figure("output/comparisons/cs6_equal_vs_favoured/fig_cs6_equal_vs_favoured_comparison", "cs6_equal_vs_favoured_comp", allow_svg=True),
        "wb_eq_0117": stage_figure("output/cs6_equal/results/cokrig_0117_0180_orn_rsn+ssn+waterblitz_tobit/fig_cokrig_england_0117", "cs6_eq_england_0117", allow_svg=False),
        "wb_eq_0180": stage_figure("output/cs6_equal/results/cokrig_0117_0180_orn_rsn+ssn+waterblitz_tobit/fig_cokrig_england_0180", "cs6_eq_england_0180", allow_svg=False),
        "wb_eq_xcorr": stage_figure("output/cs6_equal/results/cokrig_0117_0180_orn_rsn+ssn+waterblitz_tobit/fig_cokrig_cross_correlation", "cs6_eq_xcorr", allow_svg=True),
        "wb_wt_0117": stage_figure("output/cs6_weighted/results/cokrig_0117_0180_orn_rsn+ssn+waterblitz_tobit/fig_cokrig_england_0117", "cs6_wt_england_0117", allow_svg=False),
        "wb_wt_0180": stage_figure("output/cs6_weighted/results/cokrig_0117_0180_orn_rsn+ssn+waterblitz_tobit/fig_cokrig_england_0180", "cs6_wt_england_0180", allow_svg=False),
        "wb_wt_xcorr": stage_figure("output/cs6_weighted/results/cokrig_0117_0180_orn_rsn+ssn+waterblitz_tobit/fig_cokrig_cross_correlation", "cs6_wt_xcorr", allow_svg=True),
        "wb_snap_overview": stage_figure("output/station_snapping/orn_rsn+ssn+waterblitz_tobit/fig_snap_overview_orn", "cs6_wb_snap_overview", allow_svg=False),

        # Case Study 7: Nitrate (0117) Data Transformations (Raw vs Log vs Sqrt vs Box-Cox)
        "tr_diff_map": stage_figure("output/comparisons/cs7_transforms_0117/fig_cs7_transform_diff_0117", "cs7_tr_diff_map", allow_svg=False),
        "tr_raw_map": stage_figure("output/results/det_0117_orn_rsn_trans_raw/fig1_england_stream_0117_map", "cs7_tr_raw_map", allow_svg=False),
        "tr_sqrt_map": stage_figure("output/results/det_0117_orn_rsn_trans_sqrt/fig1_england_stream_0117_map", "cs7_tr_sqrt_map", allow_svg=False),
        "tr_boxcox_map": stage_figure("output/results/det_0117_orn_rsn_trans_boxcox_lam0.50/fig1_england_stream_0117_map", "cs7_tr_boxcox_map", allow_svg=False),
        "tr_side_by_side": stage_figure("output/comparisons/cs7_transforms_0117/fig_comp4_side_by_side_maps", "cs7_tr_side_by_side", allow_svg=False),
        "tr_comp_svg": stage_figure("output/comparisons/cs7_transforms_0117/fig_cs7_transform_comparison", "cs7_tr_comp_svg", allow_svg=True),
        "tr_pred_dens": stage_figure("output/comparisons/cs7_transforms_0117/fig_comp2_predictions_and_uncertainty", "cs7_tr_pred_dens", allow_svg=False),
        "tr_param_forest": stage_figure("output/comparisons/cs7_transforms_0117/fig_comp3_parameter_forest", "cs7_tr_param_forest", allow_svg=False),
    }

    payload = {
        "meta": {
            "title": "Spatial Stream Network Water Quality & Ecological Modelling Dashboard",
            "network": "OS Open Rivers Network v2 (ORN v2)",
            "reaches_england_rsn": 90876,
            "reaches_england_rss": 91118,
            "reaches_england_rsn_ssn_wb": 95773,
            "obs_rsn_wq": 1210,
            "obs_rss_wq": 1384,
            "obs_rsn_whpt": 1160,
            "obs_rsn_ssn_wb": 8346,
        },
        "figures": figures,
        "tables": {
            "cv_0117": read_csv_rows("output/results/det_0117_orn_rsn/cv_comparison_0117.csv"),
            "cv_0117_gam": read_csv_rows("output/results/det_0117_orn_rsn_trend_gam/cv_comparison_0117.csv"),
            "cv_0180": read_csv_rows("output/results/det_0180_orn_rsn/cv_comparison_0180.csv"),
            "cv_whpt": read_csv_rows("output/results/inverts_whpt_aspt_orn_rsn/cv_comparison_whpt_aspt.csv"),
            "cv_rss": read_csv_rows("output/results/det_0117_orn_rss/cv_comparison_0117.csv"),
            "params_0117": [
                r for r in read_csv_rows("output/results/det_0117_orn_rsn/parameters_0117.csv")
                if r.get("Domain") == "All England RSN Network" and ("Joint" in str(r.get("Model", "")) or "Combined" in str(r.get("Model", "")))
            ],
            "exceed_0117": read_csv_rows("output/results/det_0117_orn_rsn/exceedance_summary_0117.csv"),
            "exceed_0180": read_csv_rows("output/results/det_0180_orn_rsn/exceedance_summary_0180.csv"),
            "cs2_design_eng": [
                r for r in read_csv_rows("output/design_0117/results/det_0117_orn_rsn/optimal_monitoring_design_0117.csv")
                if r.get("Domain") == "All England RSN Network" and int(r.get("Rank", 999)) <= 20
            ],
            "cs2_design_summary": [
                r for r in read_csv_rows("output/design_0117/results/det_0117_orn_rsn/optimal_monitoring_design_0117.csv")
                if int(r.get("Rank", 0)) == 50
            ],
            "st_year_cv": [
                r for r in read_csv_rows("output/results/det_0117_orn_rsn_st/cv_summary_by_year_0117.csv")
                if r.get("Domain") == "All England RSN Network" and "Combined" in str(r.get("Model", ""))
            ],
            "st_season_cv": [
                r for r in read_csv_rows("output/results/det_0117_orn_rsn_st/cv_summary_by_season_0117.csv")
                if r.get("Domain") == "All England RSN Network" and "Combined" in str(r.get("Model", ""))
            ],
            "cokrig_cv": read_csv_rows("output/results/cokrig_0117_0180_orn_rsn/cv_univariate_vs_cokriging.csv"),
            "cokrig_cov": read_csv_rows("output/results/cokrig_0117_0180_orn_rsn/cokriging_cross_covariance.csv"),
            "cokrig_exceed": read_csv_rows("output/results/cokrig_0117_0180_orn_rsn/exceedance_summary_cokriging.csv"),
            "desc_bias": read_csv_rows("output/results/det_0117_orn_rsn/descriptives/monitoring_representativeness.csv"),
            "desc_strat": read_csv_rows("output/results/det_0117_orn_rsn/descriptives/stratified_results_summary.csv")[:16],
            "comp_tests": read_csv_rows("output/comparisons/det_0117_orn_rsn_vs_det_0180_orn_rsn/statistical_tests_summary.csv")[:18],
            "snapping_top": read_csv_rows("output/station_snapping/orn_rsn/catchment_snapping_summary.csv")[:12],
            "cs6_reach_comp": read_csv_rows("output/comparisons/cs6_equal_vs_favoured/cs6_reach_comparison_summary.csv"),
            "cs6_cv_comp": read_csv_rows("output/comparisons/cs6_equal_vs_favoured/cs6_cv_comparison.csv"),
            "cs6_cov_comp": read_csv_rows("output/comparisons/cs6_equal_vs_favoured/cs6_cross_covariance_comparison.csv"),
            "cs7_transform_summary": read_csv_rows("output/comparisons/cs7_transforms_0117/cs7_transform_summary.csv"),
            "cs7_cv_comp": read_csv_rows("output/comparisons/cs7_transforms_0117/cv_metrics_comparison.csv"),
            "cs7_pred_dist": read_csv_rows("output/comparisons/cs7_transforms_0117/prediction_distribution_summary.csv"),
        },
        "reports": {
            "cokrig_eval": read_text_excerpt("output/results/cokrig_0117_0180_orn_rsn/cokriging_evaluation_report.md", max_lines=80),
            "cs6_equal_eval": read_text_excerpt("output/cs6_equal/results/cokrig_0117_0180_orn_rsn+ssn+waterblitz_tobit/cokriging_evaluation_report.md", max_lines=80),
            "cs6_weighted_eval": read_text_excerpt("output/cs6_weighted/results/cokrig_0117_0180_orn_rsn+ssn+waterblitz_tobit/cokriging_evaluation_report.md", max_lines=80),
            "cs7_transforms_eval": read_text_excerpt("output/comparisons/cs7_transforms_0117/comparison_report.md", max_lines=80),
        },
    }

    data_js_path = os.path.join(DASH_DIR, "data.js")
    with open(data_js_path, "w", encoding="utf-8") as f:
        f.write("window.MODELLING_DASHBOARD_DATA = ")
        json.dump(payload, f, indent=2)
        f.write(";\n")
    print(f"Wrote {data_js_path} ({os.path.getsize(data_js_path) / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
