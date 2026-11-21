#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

inline uint16_t bswap16(uint16_t v) {
    return static_cast<uint16_t>((v >> 8) | (v << 8));
}

inline uint32_t bswap32(uint32_t v) {
    return ((v & 0x000000FFU) << 24) |
           ((v & 0x0000FF00U) << 8)  |
           ((v & 0x00FF0000U) >> 8)  |
           ((v & 0xFF000000U) >> 24);
}

inline uint64_t bswap64(uint64_t v) {
    return ((v & 0x00000000000000FFULL) << 56) |
           ((v & 0x000000000000FF00ULL) << 40) |
           ((v & 0x0000000000FF0000ULL) << 24) |
           ((v & 0x00000000FF000000ULL) << 8)  |
           ((v & 0x000000FF00000000ULL) >> 8)  |
           ((v & 0x0000FF0000000000ULL) >> 24) |
           ((v & 0x00FF000000000000ULL) >> 40) |
           ((v & 0xFF00000000000000ULL) >> 56);
}

inline bool read_u32(const uint8_t* p, size_t len, size_t pos, bool little_endian, uint32_t& out) {
    if (pos + 4 > len) return false;
    uint32_t raw = 0;
    std::memcpy(&raw, p + pos, 4);
#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
    out = little_endian ? bswap32(raw) : raw;
#else
    out = little_endian ? raw : bswap32(raw);
#endif
    return true;
}

inline bool read_f64(const uint8_t* p, size_t len, size_t pos, bool little_endian, double& out) {
    if (pos + 8 > len) return false;
    uint64_t raw = 0;
    std::memcpy(&raw, p + pos, 8);
#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
    uint64_t u = little_endian ? bswap64(raw) : raw;
#else
    uint64_t u = little_endian ? raw : bswap64(raw);
#endif
    std::memcpy(&out, &u, 8);
    return true;
}

inline bool parse_gpb_envelope_internal(
    const uint8_t* blob,
    size_t len,
    double& minx,
    double& maxx,
    double& miny,
    double& maxy
) {
    if (len < 40) return false;
    if (blob[0] != 0x47 || blob[1] != 0x50) return false;
    uint8_t flags = blob[3];
    if (flags & 0x10U) return false;
    uint8_t env_ind = static_cast<uint8_t>((flags >> 1) & 0x07U);
    if (env_ind < 1 || env_ind > 4) return false;
    bool le = (flags & 0x01U) != 0;
    if (!read_f64(blob, len, 8, le, minx) ||
        !read_f64(blob, len, 16, le, maxx) ||
        !read_f64(blob, len, 24, le, miny) ||
        !read_f64(blob, len, 32, le, maxy)) {
        return false;
    }
    if (std::isnan(minx) || std::isnan(maxx) || std::isnan(miny) || std::isnan(maxy)) {
        return false;
    }
    return true;
}

inline bool decode_wkb_type_info(
    uint32_t wkb_type,
    uint32_t& base_type,
    int& ndim,
    bool& has_z,
    size_t& hdr_len
) {
    bool has_srid = (wkb_type & 0x20000000U) != 0U;
    uint32_t low16 = wkb_type & 0xFFFFU;
    uint32_t iso_dim = low16 / 1000U;
    has_z = ((wkb_type & 0x80000000U) != 0U) || (low16 >= 1000U && low16 < 4000U && (iso_dim == 1U || iso_dim == 3U));
    bool has_m = ((wkb_type & 0x40000000U) != 0U) || (low16 >= 1000U && low16 < 4000U && (iso_dim == 2U || iso_dim == 3U));
    ndim = 2 + (has_z ? 1 : 0) + (has_m ? 1 : 0);
    base_type = low16 % 1000U;
    hdr_len = has_srid ? 9 : 5;
    return true;
}

bool decode_wkb_linestring_coords(
    const uint8_t* blob,
    size_t len,
    std::vector<double>& xs,
    std::vector<double>& ys,
    int32_t& z_count,
    double& z_first,
    double& z_last,
    double& z_sum
) {
    xs.clear();
    ys.clear();
    z_count = 0;
    z_first = std::numeric_limits<double>::quiet_NaN();
    z_last = std::numeric_limits<double>::quiet_NaN();
    z_sum = 0.0;

    if (!blob || len < 9) return false;
    size_t offset = 0;
    if (blob[0] == 0x47 && blob[1] == 0x50) {
        static const size_t env_sizes[8] = {0, 32, 48, 48, 64, 0, 0, 0};
        uint8_t env_ind = static_cast<uint8_t>((blob[3] >> 1) & 0x07U);
        offset = 8 + env_sizes[env_ind];
    }
    if (offset + 9 > len) return false;

    bool le = (blob[offset] == 1);
    uint32_t wkb_type = 0;
    if (!read_u32(blob, len, offset + 1, le, wkb_type)) return false;

    uint32_t base_type = 0;
    int ndim = 2;
    bool has_z = false;
    size_t hdr_len = 5;
    decode_wkb_type_info(wkb_type, base_type, ndim, has_z, hdr_len);

    auto append_coords = [&](size_t coord_pos, uint32_t npts, int g_dim, bool g_z, bool g_le) -> bool {
        size_t needed = static_cast<size_t>(npts) * static_cast<size_t>(g_dim) * 8ULL;
        if (coord_pos + needed > len) return false;
        size_t old_sz = xs.size();
        xs.resize(old_sz + npts);
        ys.resize(old_sz + npts);
        size_t stride = static_cast<size_t>(g_dim) * 8ULL;
        for (uint32_t k = 0; k < npts; ++k) {
            size_t base = coord_pos + static_cast<size_t>(k) * stride;
            double x = 0.0, y = 0.0;
            read_f64(blob, len, base, g_le, x);
            read_f64(blob, len, base + 8, g_le, y);
            xs[old_sz + k] = x;
            ys[old_sz + k] = y;
            if (g_z) {
                double zv = 0.0;
                read_f64(blob, len, base + 16, g_le, zv);
                if (!std::isnan(zv) && zv > -500.0 && zv < 9000.0) {
                    if (z_count == 0) z_first = zv;
                    z_last = zv;
                    z_sum += zv;
                    ++z_count;
                }
            }
        }
        return true;
    };

    if (base_type == 5) {
        if (offset + hdr_len + 4 > len) return false;
        uint32_t num_geoms = 0;
        if (!read_u32(blob, len, offset + hdr_len, le, num_geoms)) return false;
        size_t pos = offset + hdr_len + 4;
        for (uint32_t g = 0; g < num_geoms; ++g) {
            if (pos + 9 > len) return false;
            bool g_le = (blob[pos] == 1);
            uint32_t g_type = 0;
            if (!read_u32(blob, len, pos + 1, g_le, g_type)) return false;
            uint32_t g_base = 0;
            int g_dim = 2;
            bool g_z = false;
            size_t g_hdr = 5;
            decode_wkb_type_info(g_type, g_base, g_dim, g_z, g_hdr);
            if (pos + g_hdr + 4 > len) return false;
            uint32_t npts = 0;
            if (!read_u32(blob, len, pos + g_hdr, g_le, npts)) return false;
            size_t coord_pos = pos + g_hdr + 4;
            if (!append_coords(coord_pos, npts, g_dim, g_z, g_le)) return false;
            pos = coord_pos + static_cast<size_t>(npts) * static_cast<size_t>(g_dim) * 8ULL;
        }
        return true;
    }

    if (base_type == 2) {
        if (offset + hdr_len + 4 > len) return false;
        uint32_t npts = 0;
        if (!read_u32(blob, len, offset + hdr_len, le, npts)) return false;
        size_t coord_pos = offset + hdr_len + 4;
        return append_coords(coord_pos, npts, ndim, has_z, le);
    }

    return false;
}

inline int64_t cell_coord(double v, double cell_size) {
    return static_cast<int64_t>(std::floor(v / cell_size));
}

inline uint64_t pack_2d_key(int64_t gx, int64_t gy) {
    return (static_cast<uint64_t>(static_cast<uint32_t>(gx)) << 32) |
           static_cast<uint64_t>(static_cast<uint32_t>(gy));
}

struct FileGridKey {
    int32_t f;
    int32_t gx;
    int32_t gy;
    bool operator==(const FileGridKey& o) const {
        return f == o.f && gx == o.gx && gy == o.gy;
    }
    bool operator<(const FileGridKey& o) const {
        if (f != o.f) return f < o.f;
        if (gx != o.gx) return gx < o.gx;
        return gy < o.gy;
    }
};

struct FileGridHash {
    size_t operator()(const FileGridKey& k) const {
        uint64_t h = (static_cast<uint64_t>(static_cast<uint32_t>(k.f)) * 0x9e3779b97f4a7c15ULL) ^
                     (static_cast<uint64_t>(static_cast<uint32_t>(k.gx)) << 32) ^
                     static_cast<uint64_t>(static_cast<uint32_t>(k.gy));
        h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ULL;
        h = (h ^ (h >> 27)) * 0x94d049bb133111ebULL;
        return static_cast<size_t>(h ^ (h >> 31));
    }
};

struct GridSpan {
    int32_t start;
    int32_t end;
};

struct FlatGrid2D {
    std::vector<int32_t> items;
    std::unordered_map<uint64_t, GridSpan> spans;

    void build(int32_t n, const double* xs, const double* ys, double cell_size) {
        items.clear();
        spans.clear();
        if (n <= 0) return;

        struct Entry {
            uint64_t key;
            int32_t idx;
        };
        std::vector<Entry> entries(n);
        for (int32_t i = 0; i < n; ++i) {
            entries[i] = {pack_2d_key(cell_coord(xs[i], cell_size), cell_coord(ys[i], cell_size)), i};
        }
        std::sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) {
            if (a.key != b.key) return a.key < b.key;
            return a.idx < b.idx;
        });

        items.resize(n);
        spans.reserve(static_cast<size_t>(n) / 2 + 16);
        int32_t pos = 0;
        while (pos < n) {
            uint64_t k = entries[pos].key;
            int32_t start = pos;
            while (pos < n && entries[pos].key == k) {
                items[pos] = entries[pos].idx;
                ++pos;
            }
            spans.emplace(k, GridSpan{start, pos});
        }
    }
};

inline void point_to_polyline_d2(
    double px,
    double py,
    const std::vector<double>& xs,
    const std::vector<double>& ys,
    double& best_d2,
    double& best_along
) {
    best_d2 = 1e30;
    best_along = 0.0;
    double cum_len = 0.0;
    size_t n_pts = xs.size();
    for (size_t i = 0; i + 1 < n_pts; ++i) {
        double x1 = xs[i], y1 = ys[i];
        double x2 = xs[i + 1], y2 = ys[i + 1];
        double dx = x2 - x1;
        double dy = y2 - y1;
        double seg_len2 = dx * dx + dy * dy;
        double seg_len = std::sqrt(seg_len2);
        double t = 0.0;
        if (seg_len2 > 0.0) {
            t = ((px - x1) * dx + (py - y1) * dy) / seg_len2;
            if (t < 0.0) t = 0.0;
            else if (t > 1.0) t = 1.0;
        }
        double projx = x1 + t * dx;
        double projy = y1 + t * dy;
        double rx = px - projx;
        double ry = py - projy;
        double d2 = rx * rx + ry * ry;
        if (d2 < best_d2) {
            best_d2 = d2;
            best_along = cum_len + t * seg_len;
        }
        cum_len += seg_len;
    }
}

void build_children_csr(
    int32_t n,
    const int32_t* parent_dn,
    std::vector<int32_t>& offsets,
    std::vector<int32_t>& children
) {
    offsets.assign(static_cast<size_t>(n) + 1, 0);
    for (int32_t i = 0; i < n; ++i) {
        int32_t dn = parent_dn[i];
        if (dn >= 0 && dn < n && dn != i) {
            offsets[static_cast<size_t>(dn) + 1]++;
        }
    }
    for (int32_t i = 0; i < n; ++i) {
        offsets[static_cast<size_t>(i) + 1] += offsets[i];
    }
    children.assign(offsets[n], 0);
    std::vector<int32_t> cursor(offsets.begin(), offsets.begin() + n);
    for (int32_t i = 0; i < n; ++i) {
        int32_t dn = parent_dn[i];
        if (dn >= 0 && dn < n && dn != i) {
            children[cursor[dn]++] = i;
        }
    }
}

} // namespace

extern "C" {

int32_t s3n_geom_get_omp_max_threads() {
#ifdef _OPENMP
    return static_cast<int32_t>(omp_get_max_threads());
#else
    return 1;
#endif
}

void s3n_geom_set_omp_threads(int32_t n_threads) {
#ifdef _OPENMP
    if (n_threads > 0) {
        omp_set_num_threads(static_cast<int>(n_threads));
    }
#else
    (void)n_threads;
#endif
}

void s3n_parse_wkb_metrics_batch(
    int32_t n_blobs,
    const uint8_t* blob_buf,
    const int64_t* blob_offsets,
    double* tot_len_m,
    double* sinuosity,
    double* start_x,
    double* start_y,
    double* mid_x,
    double* mid_y,
    double* end_x,
    double* end_y,
    double* minx,
    double* maxx,
    double* miny,
    double* maxy,
    double* z_first,
    double* z_last,
    double* z_mean,
    int32_t* n_pts,
    int32_t* z_count
) {
    if (n_blobs <= 0) return;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<double> xs;
        std::vector<double> ys;
        std::vector<double> seg_lens;
        xs.reserve(128);
        ys.reserve(128);
        seg_lens.reserve(128);

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 256)
#endif
        for (int32_t i = 0; i < n_blobs; ++i) {
            int64_t b_start = blob_offsets[i];
            int64_t b_end = blob_offsets[i + 1];
            size_t blen = (b_end > b_start) ? static_cast<size_t>(b_end - b_start) : 0;

            int32_t zc = 0;
            double zf = std::numeric_limits<double>::quiet_NaN();
            double zl = std::numeric_limits<double>::quiet_NaN();
            double zs = 0.0;
            bool ok = (blen > 0) && decode_wkb_linestring_coords(
                blob_buf + b_start, blen, xs, ys, zc, zf, zl, zs
            );

            z_count[i] = zc;
            z_first[i] = zf;
            z_last[i] = zl;
            z_mean[i] = (zc > 0) ? (zs / static_cast<double>(zc)) : std::numeric_limits<double>::quiet_NaN();

            size_t np = ok ? xs.size() : 0;
            n_pts[i] = static_cast<int32_t>(np);

            if (np == 0) {
                tot_len_m[i] = 0.0;
                sinuosity[i] = 1.0;
                start_x[i] = start_y[i] = 0.0;
                mid_x[i] = mid_y[i] = 0.0;
                end_x[i] = end_y[i] = 0.0;
                minx[i] = maxx[i] = miny[i] = maxy[i] = 0.0;
                continue;
            }

            if (np == 1) {
                double x0 = xs[0], y0 = ys[0];
                tot_len_m[i] = 0.0;
                sinuosity[i] = 1.0;
                start_x[i] = mid_x[i] = end_x[i] = x0;
                start_y[i] = mid_y[i] = end_y[i] = y0;
                minx[i] = maxx[i] = x0;
                miny[i] = maxy[i] = y0;
                continue;
            }

            size_t n_seg = np - 1;
            seg_lens.resize(n_seg);
            double tot = 0.0;
            double x_lo = xs[0], x_hi = xs[0];
            double y_lo = ys[0], y_hi = ys[0];
            for (size_t s = 0; s < n_seg; ++s) {
                double x2 = xs[s + 1], y2 = ys[s + 1];
                if (x2 < x_lo) x_lo = x2;
                if (x2 > x_hi) x_hi = x2;
                if (y2 < y_lo) y_lo = y2;
                if (y2 > y_hi) y_hi = y2;
                double d = std::hypot(x2 - xs[s], y2 - ys[s]);
                seg_lens[s] = d;
                tot += d;
            }

            double half = tot * 0.5;
            double cur = 0.0;
            double mx = xs[0], my = ys[0];
            for (size_t s = 0; s < n_seg; ++s) {
                double d = seg_lens[s];
                if (cur + d >= half && d > 0.0) {
                    double t = (half - cur) / d;
                    mx = xs[s] + t * (xs[s + 1] - xs[s]);
                    my = ys[s] + t * (ys[s + 1] - ys[s]);
                    break;
                }
                cur += d;
            }

            double straight = std::hypot(xs[np - 1] - xs[0], ys[np - 1] - ys[0]);
            double sinu = (straight > 1e-3) ? std::max(1.0, tot / straight) : 1.0;

            tot_len_m[i] = tot;
            sinuosity[i] = sinu;
            start_x[i] = xs[0];
            start_y[i] = ys[0];
            mid_x[i] = mx;
            mid_y[i] = my;
            end_x[i] = xs[np - 1];
            end_y[i] = ys[np - 1];
            minx[i] = x_lo;
            maxx[i] = x_hi;
            miny[i] = y_lo;
            maxy[i] = y_hi;
        }
    }
}

void s3n_snap_stations_wkb_batch(
    int32_t n_stations,
    const double* st_x,
    const double* st_y,
    double cell_size,
    double snap_pad,
    int32_t n_blobs,
    const uint8_t* blob_buf,
    const int64_t* blob_offsets,
    int64_t chunk_base_idx,
    double* best_dist_m,
    int64_t* best_reach_global_idx,
    double* best_along_m
) {
    if (n_stations <= 0 || n_blobs <= 0) return;

    FlatGrid2D st_grid;
    st_grid.build(n_stations, st_x, st_y, cell_size);

    std::vector<double> init_d2(n_stations);
    for (int32_t s = 0; s < n_stations; ++s) {
        double d = best_dist_m[s];
        init_d2[s] = (d >= 1e15) ? 1e30 : (d * d);
    }

    int n_threads = 1;
#ifdef _OPENMP
    n_threads = std::max(1, omp_get_max_threads());
#endif
    std::vector<double> all_d2(static_cast<size_t>(n_threads) * n_stations, 1e30);
    std::vector<int64_t> all_idx(static_cast<size_t>(n_threads) * n_stations, -1);
    std::vector<double> all_along(static_cast<size_t>(n_threads) * n_stations, 0.0);

    for (int t = 0; t < n_threads; ++t) {
        size_t base = static_cast<size_t>(t) * n_stations;
        std::copy(init_d2.begin(), init_d2.end(), all_d2.begin() + base);
        std::copy(best_reach_global_idx, best_reach_global_idx + n_stations, all_idx.begin() + base);
        std::copy(best_along_m, best_along_m + n_stations, all_along.begin() + base);
    }

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        size_t tbase = static_cast<size_t>(tid) * n_stations;
        double* loc_d2 = all_d2.data() + tbase;
        int64_t* loc_idx = all_idx.data() + tbase;
        double* loc_along = all_along.data() + tbase;

        std::vector<double> xs;
        std::vector<double> ys;
        xs.reserve(128);
        ys.reserve(128);

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 128)
#endif
        for (int32_t i = 0; i < n_blobs; ++i) {
            int64_t b_start = blob_offsets[i];
            int64_t b_end = blob_offsets[i + 1];
            if (b_end <= b_start) continue;
            const uint8_t* blob = blob_buf + b_start;
            size_t blen = static_cast<size_t>(b_end - b_start);

            double minx = 0.0, maxx = 0.0, miny = 0.0, maxy = 0.0;
            bool have_env = parse_gpb_envelope_internal(blob, blen, minx, maxx, miny, maxy);
            bool decoded = false;
            int32_t zc = 0;
            double zf = 0.0, zl = 0.0, zs = 0.0;

            if (!have_env) {
                if (!decode_wkb_linestring_coords(blob, blen, xs, ys, zc, zf, zl, zs) || xs.size() < 2) {
                    continue;
                }
                decoded = true;
                minx = maxx = xs[0];
                miny = maxy = ys[0];
                for (size_t k = 1; k < xs.size(); ++k) {
                    if (xs[k] < minx) minx = xs[k];
                    if (xs[k] > maxx) maxx = xs[k];
                    if (ys[k] < miny) miny = ys[k];
                    if (ys[k] > maxy) maxy = ys[k];
                }
            }

            int64_t gx0 = cell_coord(minx - snap_pad, cell_size);
            int64_t gx1 = cell_coord(maxx + snap_pad, cell_size);
            int64_t gy0 = cell_coord(miny - snap_pad, cell_size);
            int64_t gy1 = cell_coord(maxy + snap_pad, cell_size);

            bool any_cand = false;
            for (int64_t gx = gx0; gx <= gx1 && !any_cand; ++gx) {
                for (int64_t gy = gy0; gy <= gy1; ++gy) {
                    if (st_grid.spans.find(pack_2d_key(gx, gy)) != st_grid.spans.end()) {
                        any_cand = true;
                        break;
                    }
                }
            }
            if (!any_cand) continue;

            if (!decoded) {
                if (!decode_wkb_linestring_coords(blob, blen, xs, ys, zc, zf, zl, zs) || xs.size() < 2) {
                    continue;
                }
                minx = maxx = xs[0];
                miny = maxy = ys[0];
                for (size_t k = 1; k < xs.size(); ++k) {
                    if (xs[k] < minx) minx = xs[k];
                    if (xs[k] > maxx) maxx = xs[k];
                    if (ys[k] < miny) miny = ys[k];
                    if (ys[k] > maxy) maxy = ys[k];
                }
            }

            int64_t global_i = chunk_base_idx + static_cast<int64_t>(i);
            for (int64_t gx = gx0; gx <= gx1; ++gx) {
                for (int64_t gy = gy0; gy <= gy1; ++gy) {
                    auto it = st_grid.spans.find(pack_2d_key(gx, gy));
                    if (it == st_grid.spans.end()) continue;
                    for (int32_t pos = it->second.start; pos < it->second.end; ++pos) {
                        int32_t s = st_grid.items[pos];
                        double px = st_x[s];
                        double py = st_y[s];
                        double bx = (px < minx) ? (minx - px) : ((px > maxx) ? (px - maxx) : 0.0);
                        double by = (py < miny) ? (miny - py) : ((py > maxy) ? (py - maxy) : 0.0);
                        double bbox_d2 = bx * bx + by * by;
                        if (bbox_d2 > loc_d2[s] * (1.0 + 1e-12) + 1e-12) {
                            continue;
                        }
                        double d2 = 1e30, along = 0.0;
                        point_to_polyline_d2(px, py, xs, ys, d2, along);
                        if (d2 < loc_d2[s] || (d2 == loc_d2[s] && (loc_idx[s] < 0 || global_i < loc_idx[s]))) {
                            loc_d2[s] = d2;
                            loc_idx[s] = global_i;
                            loc_along[s] = along;
                        }
                    }
                }
            }
        }
    }

    for (int32_t s = 0; s < n_stations; ++s) {
        double bd2 = init_d2[s];
        int64_t bidx = best_reach_global_idx[s];
        double balong = best_along_m[s];
        for (int t = 0; t < n_threads; ++t) {
            size_t off = static_cast<size_t>(t) * n_stations + s;
            double td2 = all_d2[off];
            int64_t tidx = all_idx[off];
            if (td2 < bd2 || (td2 == bd2 && tidx >= 0 && (bidx < 0 || tidx < bidx))) {
                bd2 = td2;
                bidx = tidx;
                balong = all_along[off];
            }
        }
        best_dist_m[s] = (bd2 >= 1e29) ? 1e18 : std::sqrt(bd2);
        best_reach_global_idx[s] = bidx;
        best_along_m[s] = balong;
    }
}

void s3n_snap_endpoints(
    int32_t n,
    const int32_t* file_idx,
    const double* start_x,
    const double* start_y,
    const double* end_x,
    const double* end_y,
    double cell_size,
    double max_dist,
    int32_t* dn_idx
) {
    if (n <= 0) return;

    struct Entry {
        FileGridKey key;
        int32_t idx;
    };
    std::vector<Entry> entries(n);
    for (int32_t i = 0; i < n; ++i) {
        entries[i] = {
            FileGridKey{
                file_idx[i],
                static_cast<int32_t>(cell_coord(start_x[i], cell_size)),
                static_cast<int32_t>(cell_coord(start_y[i], cell_size))
            },
            i
        };
    }
    std::sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) {
        if (a.key < b.key) return true;
        if (b.key < a.key) return false;
        return a.idx < b.idx;
    });

    std::vector<int32_t> items(n);
    std::unordered_map<FileGridKey, GridSpan, FileGridHash> spans;
    spans.reserve(static_cast<size_t>(n) + 16);
    int32_t pos = 0;
    while (pos < n) {
        FileGridKey k = entries[pos].key;
        int32_t start = pos;
        while (pos < n && entries[pos].key == k) {
            items[pos] = entries[pos].idx;
            ++pos;
        }
        spans.emplace(k, GridSpan{start, pos});
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256)
#endif
    for (int32_t i = 0; i < n; ++i) {
        if (dn_idx[i] >= 0) continue;
        int32_t f = file_idx[i];
        double ex = end_x[i], ey = end_y[i];
        int32_t gx0 = static_cast<int32_t>(cell_coord(ex, cell_size));
        int32_t gy0 = static_cast<int32_t>(cell_coord(ey, cell_size));
        double best_d = max_dist;
        int32_t best_cand = -1;
        for (int32_t gx = gx0 - 1; gx <= gx0 + 1; ++gx) {
            for (int32_t gy = gy0 - 1; gy <= gy0 + 1; ++gy) {
                auto it = spans.find(FileGridKey{f, gx, gy});
                if (it == spans.end()) continue;
                for (int32_t p = it->second.start; p < it->second.end; ++p) {
                    int32_t cand = items[p];
                    if (cand == i) continue;
                    double d = std::hypot(ex - start_x[cand], ey - start_y[cand]);
                    if (d < best_d || (d == best_d && (best_cand < 0 || cand < best_cand))) {
                        best_d = d;
                        best_cand = cand;
                    }
                }
            }
        }
        dn_idx[i] = best_cand;
    }
}

int32_t s3n_break_cycles(
    int32_t n,
    const int32_t* raw_parent,
    const int32_t* catch_id,
    int32_t mode,
    int32_t* out_parent_dn,
    int32_t* out_topo_order
) {
    if (n <= 0) return 0;
    int32_t cycles_broken = 0;

    if (mode == 0) {
        for (int32_t i = 0; i < n; ++i) {
            int32_t p = raw_parent[i];
            out_parent_dn[i] = (p >= 0 && p < n && p != i) ? p : -1;
        }
        std::vector<int8_t> visited_state(n, -1); // -1 = unvisited, 0 = visiting, 1 = done
        std::vector<int32_t> chain;
        chain.reserve(256);

        for (int32_t i = 0; i < n; ++i) {
            int32_t curr = i;
            chain.clear();
            while (curr >= 0 && visited_state[curr] == -1) {
                visited_state[curr] = 0;
                chain.push_back(curr);
                int32_t nxt = out_parent_dn[curr];
                if (nxt >= 0 && visited_state[nxt] == 0) {
                    out_parent_dn[curr] = -1;
                    ++cycles_broken;
                    break;
                }
                curr = nxt;
            }
            for (int32_t u : chain) {
                visited_state[u] = 1;
            }
        }

        std::vector<int32_t> ch_off, ch_idx;
        build_children_csr(n, out_parent_dn, ch_off, ch_idx);

        std::vector<int32_t> mouths;
        mouths.reserve(n / 10 + 16);
        for (int32_t i = 0; i < n; ++i) {
            if (out_parent_dn[i] < 0) mouths.push_back(i);
        }
        std::sort(mouths.begin(), mouths.end(), [&](int32_t a, int32_t b) {
            if (catch_id[a] != catch_id[b]) return catch_id[a] < catch_id[b];
            return a < b;
        });

        std::vector<int32_t> stack;
        stack.reserve(n);
        int32_t tpos = 0;
        for (int32_t m_id : mouths) {
            stack.push_back(m_id);
            while (!stack.empty()) {
                int32_t u = stack.back();
                stack.pop_back();
                out_topo_order[tpos++] = u;
                for (int32_t p = ch_off[u]; p < ch_off[u + 1]; ++p) {
                    stack.push_back(ch_idx[p]);
                }
            }
        }
        return cycles_broken;
    }

    // mode == 1 (ORN mouth-first traversal + disconnected cycle breaking)
    std::vector<int32_t> temp_parent(n, -1);
    std::vector<int32_t> mouths;
    mouths.reserve(n / 10 + 16);
    for (int32_t i = 0; i < n; ++i) {
        int32_t p = raw_parent[i];
        if (p < 0 || p >= n || p == i) {
            mouths.push_back(i);
        } else {
            temp_parent[i] = p;
        }
    }

    std::vector<int32_t> ch_off, ch_idx;
    build_children_csr(n, temp_parent.data(), ch_off, ch_idx);

    std::sort(mouths.begin(), mouths.end(), [&](int32_t a, int32_t b) {
        if (catch_id[a] != catch_id[b]) return catch_id[a] < catch_id[b];
        return a < b;
    });

    std::vector<uint8_t> visited(n, 0);
    std::fill(out_parent_dn, out_parent_dn + n, -1);
    std::vector<int32_t> stack;
    stack.reserve(n);
    int32_t tpos = 0;

    auto traverse_from_mouth = [&](int32_t m_id) {
        if (visited[m_id]) return;
        visited[m_id] = 1;
        out_parent_dn[m_id] = -1;
        stack.push_back(m_id);
        while (!stack.empty()) {
            int32_t u = stack.back();
            stack.pop_back();
            out_topo_order[tpos++] = u;
            for (int32_t ci = ch_off[u]; ci < ch_off[u + 1]; ++ci) {
                int32_t p = ch_idx[ci];
                if (!visited[p]) {
                    visited[p] = 1;
                    out_parent_dn[p] = u;
                    stack.push_back(p);
                }
            }
        }
    };

    for (int32_t m_id : mouths) {
        traverse_from_mouth(m_id);
    }

    std::vector<int32_t> rem(n);
    for (int32_t i = 0; i < n; ++i) rem[i] = i;
    std::sort(rem.begin(), rem.end(), [&](int32_t a, int32_t b) {
        if (catch_id[a] != catch_id[b]) return catch_id[a] < catch_id[b];
        return a < b;
    });
    for (int32_t r_id : rem) {
        if (!visited[r_id]) {
            ++cycles_broken;
            traverse_from_mouth(r_id);
        }
    }

    return cycles_broken;
}

int32_t s3n_check_reach_orientation(
    int32_t n,
    const int32_t* parent_dn,
    double* start_x,
    double* start_y,
    double* end_x,
    double* end_y,
    double* z_first,
    double* z_last,
    int32_t* out_flipped
) {
    if (n <= 0) return 0;
    std::vector<double> sx0(start_x, start_x + n);
    std::vector<double> sy0(start_y, start_y + n);
    std::vector<double> ex0(end_x, end_x + n);
    std::vector<double> ey0(end_y, end_y + n);

    int32_t flipped_cnt = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:flipped_cnt) schedule(static)
#endif
    for (int32_t i = 0; i < n; ++i) {
        out_flipped[i] = 0;
        int32_t dn = parent_dn[i];
        if (dn < 0 || dn >= n || dn == i) continue;

        double dn_sx = sx0[dn], dn_sy = sy0[dn];
        double dn_ex = ex0[dn], dn_ey = ey0[dn];

        double d_end_to_dn = std::min(
            std::hypot(ex0[i] - dn_sx, ey0[i] - dn_sy),
            std::hypot(ex0[i] - dn_ex, ey0[i] - dn_ey)
        );
        double d_start_to_dn = std::min(
            std::hypot(sx0[i] - dn_sx, sy0[i] - dn_sy),
            std::hypot(sx0[i] - dn_ex, sy0[i] - dn_ey)
        );

        if (d_start_to_dn + 10.0 < d_end_to_dn) {
            out_flipped[i] = 1;
            start_x[i] = ex0[i];
            start_y[i] = ey0[i];
            end_x[i] = sx0[i];
            end_y[i] = sy0[i];
            if (z_first != nullptr && z_last != nullptr) {
                std::swap(z_first[i], z_last[i]);
            }
            flipped_cnt += 1;
        }
    }
    return flipped_cnt;
}

void s3n_accumulate_tree_hydro(
    int32_t n,
    const int32_t* parent_dn,
    const int32_t* topo_order,
    const double* length_km,
    const double* raw_d2m_km,
    const int32_t* raw_shreve,
    const int32_t* raw_strahler,
    const double* raw_us_accum_km,
    double* out_d2m_km,
    int32_t* out_shreve,
    int32_t* out_strahler,
    double* out_us_accum_km
) {
    if (n <= 0) return;

    // Downstream-to-upstream pass: Dist2Mth_km
    for (int32_t idx = 0; idx < n; ++idx) {
        int32_t u = topo_order[idx];
        double r_d2m = raw_d2m_km[u];
        if (!std::isnan(r_d2m) && r_d2m > 0.0) {
            out_d2m_km[u] = r_d2m;
        } else {
            int32_t dn = parent_dn[u];
            double base_d2m = (dn >= 0 && dn < n) ? out_d2m_km[dn] : 0.0;
            out_d2m_km[u] = base_d2m + length_km[u];
        }
    }

    std::vector<int32_t> ch_off, ch_idx;
    build_children_csr(n, parent_dn, ch_off, ch_idx);

    // Upstream-to-downstream pass: Shreve, Strahler, US_Accum_km
    for (int32_t idx = n - 1; idx >= 0; --idx) {
        int32_t u = topo_order[idx];
        int32_t cs = ch_off[u];
        int32_t ce = ch_off[u + 1];

        // Shreve
        if (raw_shreve[u] > 0) {
            out_shreve[u] = raw_shreve[u];
        } else if (cs == ce) {
            out_shreve[u] = 1;
        } else {
            int32_t s_sum = 0;
            for (int32_t ci = cs; ci < ce; ++ci) s_sum += out_shreve[ch_idx[ci]];
            out_shreve[u] = s_sum;
        }

        // Strahler
        if (raw_strahler[u] > 0) {
            out_strahler[u] = raw_strahler[u];
        } else if (cs == ce) {
            out_strahler[u] = 1;
        } else {
            int32_t mx_s = 0;
            int32_t cnt_mx = 0;
            for (int32_t ci = cs; ci < ce; ++ci) {
                int32_t st = out_strahler[ch_idx[ci]];
                if (st > mx_s) {
                    mx_s = st;
                    cnt_mx = 1;
                } else if (st == mx_s) {
                    cnt_mx += 1;
                }
            }
            out_strahler[u] = (cnt_mx >= 2) ? (mx_s + 1) : mx_s;
        }

        // US_Accum_km
        double r_afv = raw_us_accum_km[u];
        if (!std::isnan(r_afv) && r_afv > 0.0) {
            out_us_accum_km[u] = std::max(1e-4, r_afv);
        } else {
            double acc = length_km[u];
            for (int32_t ci = cs; ci < ce; ++ci) acc += out_us_accum_km[ch_idx[ci]];
            out_us_accum_km[u] = acc;
        }
    }
}

void s3n_compute_log_afv_tree(
    int32_t n,
    const int32_t* parent_dn,
    const int32_t* topo_order,
    const double* us_accum_km,
    const int32_t* shreve,
    const double* elev,
    const int32_t* is_aq,
    int32_t afv_scheme,
    double* out_log_afv
) {
    if (n <= 0) return;

    std::vector<int32_t> ch_off, ch_idx;
    build_children_csr(n, parent_dn, ch_off, ch_idx);

    std::fill(out_log_afv, out_log_afv + n, 0.0);

    auto reach_weight = [&](int32_t p) -> double {
        double us_km = std::max(1e-4, us_accum_km[p]);
        if (afv_scheme == 1) {
            int32_t sh = (shreve != nullptr) ? shreve[p] : 1;
            return std::max(1.0, static_cast<double>(sh > 0 ? sh : 1));
        }
        if (afv_scheme == 2) {
            double el = (elev != nullptr && !std::isnan(elev[p])) ? std::max(0.0, elev[p]) : 0.0;
            double aq = (is_aq != nullptr && is_aq[p] > 0) ? static_cast<double>(is_aq[p]) : 0.0;
            return std::max(1e-4, us_km * (1.0 + 0.002 * el) * (1.0 + 0.50 * aq));
        }
        return us_km;
    };

    for (int32_t idx = 0; idx < n; ++idx) {
        int32_t u = topo_order[idx];
        if (parent_dn[u] < 0) out_log_afv[u] = 0.0;
        int32_t cs = ch_off[u];
        int32_t ce = ch_off[u + 1];
        if (cs >= ce) continue;

        double sum_w = 0.0;
        for (int32_t ci = cs; ci < ce; ++ci) {
            sum_w += reach_weight(ch_idx[ci]);
        }
        if (sum_w <= 0.0) continue;
        for (int32_t ci = cs; ci < ce; ++ci) {
            int32_t p = ch_idx[ci];
            double w = reach_weight(p);
            out_log_afv[p] = out_log_afv[u] + std::log(w / sum_w);
        }
    }
}

void s3n_impute_tree_index(
    int32_t n,
    const int32_t* parent_dn,
    const int32_t* topo_order,
    int32_t n_passes,
    int32_t* cov_src
) {
    if (n <= 0 || n_passes <= 0) return;

    std::vector<int32_t> ch_off, ch_idx;
    build_children_csr(n, parent_dn, ch_off, ch_idx);

    for (int32_t pass = 0; pass < n_passes; ++pass) {
        for (int32_t idx = 0; idx < n; ++idx) {
            int32_t u = topo_order[idx];
            if (cov_src[u] < 0) {
                int32_t dn = parent_dn[u];
                if (dn >= 0 && dn < n && cov_src[dn] >= 0) {
                    cov_src[u] = cov_src[dn];
                }
            }
        }
        for (int32_t idx = n - 1; idx >= 0; --idx) {
            int32_t u = topo_order[idx];
            if (cov_src[u] < 0) {
                for (int32_t ci = ch_off[u]; ci < ch_off[u + 1]; ++ci) {
                    int32_t p = ch_idx[ci];
                    if (cov_src[p] >= 0) {
                        cov_src[u] = cov_src[p];
                        break;
                    }
                }
            }
        }
    }
}

void s3n_spatial_hash_query_batch(
    int32_t n_ref,
    const double* ref_x,
    const double* ref_y,
    double cell_size,
    int32_t n_q,
    const double* q_x,
    const double* q_y,
    int32_t* out_ref_idx
) {
    if (n_q <= 0) return;
    if (n_ref <= 0) {
        std::fill(out_ref_idx, out_ref_idx + n_q, -1);
        return;
    }

    FlatGrid2D ref_grid;
    ref_grid.build(n_ref, ref_x, ref_y, cell_size);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 512)
#endif
    for (int32_t i = 0; i < n_q; ++i) {
        double mx = q_x[i];
        double my = q_y[i];
        int64_t gx0 = cell_coord(mx, cell_size);
        int64_t gy0 = cell_coord(my, cell_size);

        double best_d2 = 1e18;
        int32_t best_idx = -1;
        const int64_t radii[2] = {1, 3};
        for (int r_i = 0; r_i < 2; ++r_i) {
            int64_t rad = radii[r_i];
            for (int64_t gx = gx0 - rad; gx <= gx0 + rad; ++gx) {
                for (int64_t gy = gy0 - rad; gy <= gy0 + rad; ++gy) {
                    auto it = ref_grid.spans.find(pack_2d_key(gx, gy));
                    if (it == ref_grid.spans.end()) continue;
                    for (int32_t pos = it->second.start; pos < it->second.end; ++pos) {
                        int32_t ref_i = ref_grid.items[pos];
                        double dx = mx - ref_x[ref_i];
                        double dy = my - ref_y[ref_i];
                        double d2 = dx * dx + dy * dy;
                        if (d2 < best_d2 || (d2 == best_d2 && (best_idx < 0 || ref_i < best_idx))) {
                            best_d2 = d2;
                            best_idx = ref_i;
                        }
                    }
                }
            }
            if (best_idx >= 0) break;
        }
        out_ref_idx[i] = best_idx;
    }
}

} // extern "C"
