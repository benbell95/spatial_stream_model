#!/usr/bin/env python3
"""Unit tests for shared Python utilities in scripts/s3n_py_utils.py and preprocessing logic."""

import math
import os
import struct
import sys
import tempfile
import unittest

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SCRIPTS_DIR = os.path.join(_REPO_ROOT, "scripts")
if _SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIR)

import s3n_py_utils as u  # noqa: E402
import format_custom_data as fcd  # noqa: E402
import importlib  # noqa: E402
prep_orn = importlib.import_module("01_preprocess_orn")  # noqa: E402


def _make_wkb_linestring(coords, has_z=False):
    bo = "<"
    wkb_type = 1002 if has_z else 2
    parts = [struct.pack(bo + "BII", 1, wkb_type, len(coords))]
    for pt in coords:
        if has_z:
            parts.append(struct.pack(bo + "3d", pt[0], pt[1], pt[2]))
        else:
            parts.append(struct.pack(bo + "2d", pt[0], pt[1]))
    return b"".join(parts)


def _make_gpkg_blob(wkb_bytes, bbox=None, srs_id=27700):
    if bbox is None:
        header = struct.pack("<2sBBi", b"GP", 0, 0x01, srs_id)
    else:
        minx, maxx, miny, maxy = bbox
        flags = 0x01 | (1 << 1)  # little-endian, 2D envelope (4 doubles)
        header = struct.pack("<2sBBi4d", b"GP", 0, flags, srs_id, minx, maxx, miny, maxy)
    return header + wkb_bytes


class TestCoordsAndDates(unittest.TestCase):
    def test_wgs84_to_bng_greenwich(self):
        # Royal Observatory Greenwich ~ (51.4769, -0.0005) -> E~538850, N~177400
        e, n = u.wgs84_to_bng(51.4769, -0.0005)
        self.assertTrue(537000 < e < 541000, f"Unexpected Easting: {e}")
        self.assertTrue(175000 < n < 180000, f"Unexpected Northing: {n}")

    def test_transform_coords_batch_identity_and_roundtrip(self):
        pts = [(400000.0, 300000.0), (530000.0, 180000.0)]
        self.assertEqual(u.transform_coords_batch(pts, 27700, 27700), pts)
        wgs = u.transform_coords_batch(pts, 27700, 4326)
        back = u.transform_coords_batch(wgs, 4326, 27700)
        for (x0, y0), (x1, y1) in zip(pts, back):
            self.assertAlmostEqual(x0, x1, delta=2.0)
            self.assertAlmostEqual(y0, y1, delta=2.0)

    def test_parse_datetime_iso(self):
        self.assertEqual(u.parse_datetime_iso("2024-03-15"), "2024-03-15T12:00:00")
        self.assertEqual(u.parse_datetime_iso("15/01/2023", "09:30"), "2023-01-15T09:30:00")
        self.assertEqual(u.parse_datetime_iso("2022-08-05T14:22:00Z"), "2022-08-05T14:22:00")
        self.assertEqual(u.parse_datetime_iso(""), "2023-06-01T12:00:00")

    def test_parse_month_arg(self):
        self.assertEqual(u.parse_month_arg("jan"), 1)
        self.assertEqual(u.parse_month_arg("December"), 12)
        self.assertEqual(u.parse_month_arg("07"), 7)
        with self.assertRaises(ValueError):
            u.parse_month_arg("13")


class TestMeasurementTokensAndTobit(unittest.TestCase):
    def test_parse_raw_measurement_token(self):
        p = u.parse_raw_measurement_token("12.5")
        self.assertEqual(p["type"], "point")
        self.assertAlmostEqual(p["val"], 12.5)

        left = u.parse_raw_measurement_token("<0.4 mg/l")
        self.assertEqual(left["type"], "left")
        self.assertEqual(left["b"], 0.4)

        right = u.parse_raw_measurement_token(">50")
        self.assertEqual(right["type"], "right")
        self.assertEqual(right["a"], 50.0)

        iv = u.parse_raw_measurement_token("10 - 25")
        self.assertEqual(iv["type"], "interval")
        self.assertEqual((iv["a"], iv["b"]), (10.0, 25.0))

        iv_rev = u.parse_raw_measurement_token("20 to 5")
        self.assertEqual((iv_rev["a"], iv_rev["b"]), (5.0, 20.0))

        self.assertIsNone(u.parse_raw_measurement_token("NO_WATER_AVAILABLE"))
        self.assertIsNone(u.parse_raw_measurement_token(""))

        sci_left = u.parse_raw_measurement_token("<1e-3")
        self.assertEqual(sci_left["type"], "left")
        self.assertAlmostEqual(sci_left["b"], 1e-3, places=8)

        sci_iv = u.parse_raw_measurement_token("1.5e-2 - 3.5e-2")
        self.assertEqual(sci_iv["type"], "interval")
        self.assertAlmostEqual(sci_iv["a"], 0.015, places=8)
        self.assertAlmostEqual(sci_iv["b"], 0.035, places=8)

        sci_right = u.parse_raw_measurement_token(">2.5e2")
        self.assertEqual(sci_right["type"], "right")
        self.assertAlmostEqual(sci_right["a"], 250.0, places=6)

    def test_interval_tobit_model(self):
        counts = {
            ("point", 5.0, 5.0): 20,
            ("point", 10.0, 10.0): 20,
            ("point", 15.0, 15.0): 15,
            ("interval", 5.0, 15.0): 25,
            ("left", 0.0, 2.0): 5,
            ("right", 30.0, None): 5,
        }
        mod = u.IntervalTobitModel("Nitrate")
        mod.fit_from_token_counts(counts)
        self.assertTrue(mod.fitted)

        pt = mod.predict("point", 8.0, 8.0, 8.0)
        self.assertEqual((pt[0], pt[5]), (8.0, "exact_point"))

        iv = mod.predict("interval", 5.0, 15.0)
        self.assertTrue(5.0 <= iv[0] <= 15.0)
        self.assertEqual(iv[5], "tobit_interval")
        self.assertGreater(iv[4], 0.0)

        lf = mod.predict("left", 0.0, 2.0)
        self.assertTrue(0.0 < lf[0] <= 2.0)

        rt = mod.predict("right", 30.0, None)
        self.assertGreaterEqual(rt[0], 30.0)
        self.assertIsNone(rt[3])

        lf0 = mod.predict("left", 0.0, 0.0)
        self.assertTrue(math.isfinite(lf0[0]))

        # Zero-anchored [0, b] interval must behave as left-censored (< b), not fail log(0)
        counts_zero_iv = dict(counts)
        counts_zero_iv[("interval", 0.0, 3.0)] = 10
        mod2 = u.IntervalTobitModel("ZeroIntervalTest")
        mod2.fit_from_token_counts(counts_zero_iv)
        self.assertTrue(mod2.fitted)
        iv_zero = mod2.predict("interval", 0.0, 3.0)
        self.assertEqual(iv_zero[5], "tobit_interval")
        self.assertEqual(iv_zero[2], 0.0)
        self.assertTrue(0.0 < iv_zero[0] <= 3.0)

    def test_match_known_determinand(self):
        d_n = fcd.match_known_determinand("Nitrate (mg/L)")
        self.assertEqual(d_n["code"], "0117")
        self.assertEqual(d_n["unit"], "mg/L")
        d_p = fcd.match_known_determinand("Phosphate (mg/L)")
        self.assertEqual(d_p["code"], "0180")
        d_ph = fcd.match_known_determinand("pH")
        self.assertEqual(d_ph["code"], "0061")
        self.assertEqual(fcd.match_known_determinand("Field pH")["code"], "0061")
        self.assertEqual(fcd.match_known_determinand("Water Temp (C)")["code"], "0076")
        self.assertEqual(fcd.match_known_determinand("Organic Carbon")["code"], "0301")
        # Avoid false substring collision on "ph" in unrelated names like "phosphorus"/"graphite"
        self.assertNotEqual(fcd.match_known_determinand("graphite_index").get("code"), "0061")

    def test_compute_alternative_range_result(self):
        tok = u.parse_raw_measurement_token("4 - 16")
        self.assertAlmostEqual(u.compute_alternative_range_result(tok, "geometric")[0], 8.0, places=4)
        self.assertAlmostEqual(u.compute_alternative_range_result(tok, "midpoint")[0], 10.0, places=4)
        self.assertAlmostEqual(u.compute_alternative_range_result(tok, "lower")[0], 4.0, places=4)
        self.assertAlmostEqual(u.compute_alternative_range_result(tok, "upper")[0], 16.0, places=4)

    def test_tobit_no_double_imputation_guard(self):
        # Verify 01_preprocess_orn preserves Tobit-imputed left-censored records via meas_method flag
        with open(os.path.join(_SCRIPTS_DIR, "01_preprocess_orn.py"), "r", encoding="utf-8") as f:
            src = f.read()
        self.assertIn('"tobit" not in meth_val', src)


class TestSourcesAndRangeDetection(unittest.TestCase):
    def test_resolve_source_tag(self):
        reg = u.load_source_registry(_REPO_ROOT)
        self.assertEqual(u.resolve_source_tag("rsn", reg, base_dir=_REPO_ROOT)[0], "rsn")
        self.assertEqual(u.resolve_source_tag("rss", reg, base_dir=_REPO_ROOT)[0], "rss")
        key, tokens, abbr, _ = u.resolve_source_tag("waterblitz_tobit+rsn", reg, base_dir=_REPO_ROOT)
        self.assertEqual(key, "rsn+waterblitz_tobit")
        self.assertEqual(tokens, ["rsn", "waterblitz_tobit"])
        self.assertIn("RSN", abbr)

        key_hyph, tokens_hyph, _, _ = u.resolve_source_tag("rsn,waterblitz-tobit", reg, base_dir=_REPO_ROOT)
        self.assertEqual(key_hyph, "rsn+waterblitz_tobit")
        self.assertEqual(tokens_hyph, ["rsn", "waterblitz_tobit"])

        # Range-mode suffix on a registered source without dedicated subfolder should resolve to single base constituent
        key_rsn_tobit, toks_rsn_tobit, _, _ = u.resolve_source_tag("rsn_tobit", reg, base_dir=_REPO_ROOT)
        self.assertEqual(key_rsn_tobit, "rsn_tobit")
        self.assertEqual(toks_rsn_tobit, ["rsn"])

    def test_detect_range_mode_avoids_false_substring_match(self):
        sources_reg, _ = u.load_source_registry(_REPO_ROOT)
        self.assertEqual(u.detect_range_mode("waterblitz_tobit", sources_reg, base_dir=_REPO_ROOT), "tobit")
        self.assertEqual(u.detect_range_mode("custom-midpoint", sources_reg, base_dir=_REPO_ROOT), "midpoint")
        self.assertIsNone(u.detect_range_mode("custom_midpoint_survey", sources_reg, base_dir=_REPO_ROOT))
        self.assertIsNone(u.detect_range_mode("rsn", sources_reg, cli_mode="none", base_dir=_REPO_ROOT))


class TestGeometryCensoringAndTransforms(unittest.TestCase):
    def test_gpkg_envelope_and_geometry(self):
        pts_in = [(400000.0, 300000.0, 52.5), (400500.0, 300400.0, 48.0)]
        wkb3d = _make_wkb_linestring(pts_in, has_z=True)
        gpb = _make_gpkg_blob(wkb3d, bbox=(400000.0, 400500.0, 300000.0, 300400.0))

        self.assertEqual(u.parse_gpkg_envelope(gpb), (400000.0, 400500.0, 300000.0, 300400.0))
        self.assertIsNone(u.parse_gpkg_envelope(wkb3d))

        pts_out, zs = u.parse_gpkg_geometry(gpb)
        self.assertEqual(pts_out, [(400000.0, 300000.0), (400500.0, 300400.0)])
        self.assertEqual(zs, [52.5, 48.0])

        pts_raw, zs_raw = u.parse_gpkg_geometry(wkb3d)
        self.assertEqual(pts_raw, pts_out)
        self.assertEqual(zs_raw, zs)

    def test_point_to_polyline_and_metrics(self):
        poly = [(0.0, 0.0), (100.0, 0.0), (100.0, 100.0)]
        d, along = u.point_to_polyline(50.0, 10.0, poly)
        self.assertAlmostEqual(d, 10.0, places=5)
        self.assertAlmostEqual(along, 50.0, places=5)

        tot, sinu, coords = u.compute_polyline_metrics(poly)
        self.assertAlmostEqual(tot, 200.0, places=5)
        self.assertGreater(sinu, 1.3)
        self.assertAlmostEqual(coords[2], 100.0, places=5)
        self.assertAlmostEqual(coords[3], 0.0, places=5)

    def test_ros_censor_imputation(self):
        raw_by_site = {
            "S1": [(0.2, True), (0.5, True), (1.0, False), (2.0, False), (4.0, False), (8.0, False)],
            "S2": [(0.3, True), (1.5, False)],
        }
        hl = u.impute_censored_by_site(raw_by_site, method="half_lod")
        self.assertAlmostEqual(hl["S1"][0], 0.1)

        ros = u.impute_censored_by_site(raw_by_site, method="ros")
        self.assertTrue(0.0 < ros["S1"][0] <= 0.2)
        self.assertAlmostEqual(ros["S1"][3], 2.0)

        ros_site = u.impute_censored_by_site(raw_by_site, method="ros_site")
        self.assertTrue(0.0 < ros_site["S2"][0] <= 0.3)

    def test_transforms_and_runtime_config(self):
        vals = [1.0, 2.0, 4.0, 8.0, 16.0]
        lam = u.estimate_boxcox_lambda(vals)
        self.assertTrue(-0.5 <= lam <= 0.5)
        self.assertAlmostEqual(u.apply_forward_transform(math.e, "log"), 1.0, places=5)
        self.assertAlmostEqual(u.apply_forward_transform(100.0, "log_10"), 2.0, places=5)
        self.assertAlmostEqual(u.apply_forward_transform(0.0, "log_1p"), 0.0, places=5)
        self.assertAlmostEqual(u.apply_forward_transform(9.0, "sqrt"), 3.0, places=5)
        self.assertAlmostEqual(u.apply_forward_transform(2.0, "inverse"), -0.5, places=5)
        self.assertAlmostEqual(u.apply_forward_transform(1.0, "box_cox", lam=0.5), 0.0, places=5)
        self.assertAlmostEqual(u.apply_forward_transform(0.0, "yeo_johnson", lam=0.5), 0.0, places=5)

        with tempfile.TemporaryDirectory() as tmpd:
            cfg = u.build_runtime_config(
                base_dir=_REPO_ROOT,
                gis_dir=os.path.join(tmpd, "gis"),
                data_dir=os.path.join(tmpd, "data"),
                out_dir=os.path.join(tmpd, "out"),
                create_out_dir=True,
            )
            self.assertEqual(cfg.gis_dir, os.path.join(tmpd, "gis"))
            self.assertTrue(os.path.isdir(cfg.processed_dir))

    def test_compute_afv_reach_weight(self):
        w_len = u.compute_afv_reach_weight("length", us_accum_km=12.5, shreve=4, elevation=200.0, is_aquifer=1)
        w_shr = u.compute_afv_reach_weight("shreve", us_accum_km=12.5, shreve=4, elevation=200.0, is_aquifer=1)
        w_hyd = u.compute_afv_reach_weight("hydro_yield", us_accum_km=12.5, shreve=4, elevation=200.0, is_aquifer=1)
        self.assertAlmostEqual(w_len, 12.5, places=5)
        self.assertAlmostEqual(w_shr, 4.0, places=5)
        self.assertAlmostEqual(w_hyd, 12.5 * (1.0 + 0.002 * 200.0) * 1.5, places=5)

    def test_yeojohnson_negative_jacobian_sign(self):
        # Negative-only sample skewed towards zero should behave sensibly and remain finite
        # Previously wrong Jacobian sign (-log(-y+1)) pushed lambda to extremes
        neg_vals = [-1.0, -2.0, -4.0, -8.0, -16.0]
        lam = u.estimate_yeojohnson_lambda(neg_vals)
        self.assertTrue(-2.0 <= lam <= 3.0)
        # Positive skewed sample should prefer Box-Cox-like lambda near 0 (not ~2.0)
        pos_skew = [1.0, 2.0, 4.0, 8.0, 16.0]
        lam_pos = u.estimate_yeojohnson_lambda(pos_skew)
        self.assertTrue(-1.0 <= lam_pos <= 1.0)

    def test_summarize_reach_site_obs_canonical_mean_and_se(self):
        # Site A has 3 samples at along_m=100m, Site B has 1 sample at along_m=500m -> weighted along_m = 200m -> obs_up_dist_km = 10 - 0.2 = 9.8
        s_list = [
            ("SA", 10.0, 100.0, [10.0, 30.0, 50.0], [0.0, 0.3, 0.4], ["exact_point", "tobit_interval"]),
            ("SB", 20.0, 500.0, [90.0], [0.0], ["exact_point"]),
        ]
        row = prep_orn.summarize_reach_site_obs(s_list, d2m_km=10.0, trans_type="logit", eff_scale=100.0)
        self.assertEqual(row[0], 1)
        self.assertAlmostEqual(row[2], 10.0, places=2)
        self.assertAlmostEqual(row[3], 9.8, places=4)
        # Canonical trans_obs_mean = mean of logit(y / 100), not logit(mean(y) / 100)
        expected_trans = sum(u.apply_forward_transform(v, "logit", scale=100.0) for v in [10.0, 30.0, 50.0, 90.0]) / 4.0
        self.assertAlmostEqual(row[6], expected_trans, places=4)
        # Pooled SE of sample mean = sqrt(0.3^2 + 0.4^2) / 4 = 0.5 / 4 = 0.125
        self.assertAlmostEqual(row[7], 0.125, places=5)
        self.assertEqual(row[8], "tobit_interval")
        self.assertEqual(row[9], 4)

    def test_format_custom_dataset_bbox_in_bng(self):
        import contextlib
        import io
        self.assertEqual(u.parse_datetime_iso("2023-06-01 02:30:00 PM"), "2023-06-01T14:30:00")
        self.assertEqual(u.parse_datetime_iso("2023-06-01 12:30:00 AM"), "2023-06-01T00:30:00")
        with tempfile.TemporaryDirectory() as tmpd:
            in_csv = os.path.join(tmpd, "raw_bng.csv")
            out_dir = os.path.join(tmpd, "formatted")
            with open(in_csv, "w", encoding="utf-8") as f:
                f.write("site_id,easting,northing,date,nitrate\n")
                f.write("S_IN,450000,350000,2024-06-01,12.5\n")
                f.write("S_OUT,200000,100000,2024-06-01,9.0\n")
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                res = fcd.format_custom_dataset(
                    input_path=in_csv,
                    output_dir=out_dir,
                    source_name="bng_bbox_test",
                    bbox="400000,300000,500000,400000",
                    verbose=True,
                )
            self.assertIsNotNone(res)
            self.assertEqual(res["n_samples"], 1)
            self.assertIn("Valid Samples Written", buf.getvalue())
            with open(os.path.join(out_dir, "sites.csv"), "r", encoding="utf-8") as f:
                sites_txt = f.read()
            self.assertIn("S_IN", sites_txt)
            self.assertNotIn("S_OUT", sites_txt)

    def test_gdal_library_candidates(self):
        mac = prep_orn._gdal_library_candidates("darwin")
        self.assertIn("libgdal.dylib", mac)
        self.assertIn("/opt/homebrew/lib/libgdal.dylib", mac)
        lin = prep_orn._gdal_library_candidates("linux")
        self.assertIn("libgdal.so", lin)
        win = prep_orn._gdal_library_candidates("win32")
        self.assertIn("gdal.dll", win)


if __name__ == "__main__":
    unittest.main()
