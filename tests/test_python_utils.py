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
        self.assertAlmostEqual(u.apply_forward_transform(9.0, "sqrt"), 3.0, places=5)
        self.assertAlmostEqual(u.apply_forward_transform(0.0, "yeojohnson", lam=0.5), 0.0, places=5)

        with tempfile.TemporaryDirectory() as tmpd:
            cfg = u.build_runtime_config(
                base_dir=_REPO_ROOT,
                gis_dir=os.path.join(tmpd, "gis"),
                data_dir=os.path.join(tmpd, "data"),
                out_dir=os.path.join(tmpd, "out"),
                create_out_dir=True,
            )
            self.assertEqual(cfg.gis_dir, os.path.join(tmpd, "gis"))
            self.assertTrue(os.path.isdir(cfg.output_dir))


if __name__ == "__main__":
    unittest.main()
