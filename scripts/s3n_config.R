# ==============================================================================
# s3n_config.R
#
# Shared path resolution, directory configs, JSON parsing, data source
# registry, water network metadata resolution, and determinand catalog helpers.
# Automatically sourced by scripts/s3n_plot_functions.R.
# ==============================================================================

# ==============================================================================
# RELATIVE PATH & CUSTOM DIRECTORY RESOLUTION HELPERS
# ==============================================================================

resolve_base_dir <- function() {
  if (file.exists("scripts/s3n_plot_functions.R") || dir.exists("configs/networks")) {
    return(".")
  }
  if (file.exists("s3n_plot_functions.R") && dir.exists("../configs/networks")) {
    return("..")
  }
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg  <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_p  <- sub("^--file=", "", file_arg[1])
    cand_root <- dirname(dirname(script_p))
    if (nzchar(cand_root) && dir.exists(cand_root)) {
      return(cand_root)
    }
  }
  "."
}

rel_path <- function(base_dir, ...) {
  parts <- c(list(...))
  parts <- as.character(unlist(parts))
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return(base_dir)
  p <- if (is.null(base_dir) || !nzchar(base_dir) || base_dir == ".") {
    do.call(file.path, as.list(parts))
  } else {
    do.call(file.path, as.list(c(base_dir, parts)))
  }
  sub("^\\./", "", p)
}

is_abs_path <- function(p) {
  grepl("^(/|\\\\\\\\|[A-Za-z]:[/\\\\])", as.character(p)[1])
}

format_rel_display_path <- function(p, base_dir = ".") {
  p_clean <- sub("^\\./", "", p)
  if (!is.null(base_dir) && nzchar(base_dir) && base_dir != ".") {
    b_clean <- gsub("[/\\\\]+$", "", base_dir)
    for (sep in c("/", "\\")) {
      prefix <- paste0(b_clean, sep)
      if (startsWith(p_clean, prefix)) {
        p_clean <- substr(p_clean, nchar(prefix) + 1L, nchar(p_clean))
        break
      }
    }
  }
  p_clean
}

resolve_python_bin <- function() {
  for (cand in c("python3", "python")) {
    p <- unname(Sys.which(cand))
    if (nzchar(p) && !grepl("WindowsApps", p, ignore.case = TRUE)) {
      return(list(bin = p, prefix = character(0)))
    }
  }
  if (.Platform$OS.type == "windows") {
    py_launcher <- unname(Sys.which("py"))
    if (nzchar(py_launcher)) {
      return(list(bin = py_launcher, prefix = "-3"))
    }
  }
  list(bin = "python3", prefix = character(0))
}

run_python_script <- function(script_path, script_args = character(0), stdout = "", stderr = "") {
  py <- resolve_python_bin()
  system2(py$bin, args = c(py$prefix, shQuote(script_path), script_args), stdout = stdout, stderr = stderr)
}

# ==============================================================================
# TITLE & TEXT WRAPPING ENGINE
# ==============================================================================

wrap_title <- function(text, max_len = 65) {
  if (is.null(text) || length(text) == 0 || all(is.na(text))) return("")
  if (length(text) > 1) {
    return(vapply(text, function(t) wrap_title(t, max_len = max_len), character(1), USE.NAMES = FALSE))
  }
  str <- as.character(text)[1]
  if (is.na(str) || !nzchar(str)) return("")
  if (is.null(max_len) || is.na(max_len) || max_len <= 0) return(str)

  lines <- strsplit(str, "\n", fixed = TRUE)[[1]]
  wrapped_lines <- unlist(lapply(lines, function(l) {
    if (nchar(l) <= max_len) return(l)
    w <- strwrap(l, width = max_len)
    if (length(w) > 0) paste(w, collapse = "\n") else l
  }))
  paste(wrapped_lines, collapse = "\n")
}

# ==============================================================================
# DIRECTORY CONFIGURATION HELPERS
# ==============================================================================

read_json_safe <- function(filepath) {
  if (!file.exists(filepath)) return(NULL)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    return(tryCatch(jsonlite::fromJSON(filepath, simplifyVector = FALSE), error = function(e) NULL))
  }
  NULL
}

parse_json_string_field <- function(lines, key, default = "") {
  pat <- paste0('"', key, '"\\s*:\\s*"([^"]*)"')
  hit <- grep(pat, lines, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0('^.*', pat, '.*$'), "\\1", hit[1])
}

parse_json_int_field <- function(lines, key, default = NA_integer_) {
  pat <- paste0('"', key, '"\\s*:\\s*([0-9]+)')
  hit <- grep(pat, lines, value = TRUE)
  if (length(hit) == 0) return(default)
  as.integer(sub(paste0('^.*', pat, '.*$'), "\\1", hit[1]))
}

find_directory_config_file <- function(base_dir = ".") {
  rel_path(base_dir, "configs", "directories.json")
}

load_directory_config <- function(base_dir = ".", auto_create = FALSE) {
  cfg_file <- find_directory_config_file(base_dir)
  if (!file.exists(cfg_file) && isTRUE(auto_create)) {
    dir.create(dirname(cfg_file), recursive = TRUE, showWarnings = FALSE)
    writeLines(c('{',
                 '  "data_dir": "data",',
                 '  "gis_dir": "gis",',
                 '  "output_dir": "output"',
                 '}'), cfg_file)
  }

  raw_cfg <- list(data_dir = "data", gis_dir = "gis", output_dir = "output", processed_dir = "")
  j_obj <- read_json_safe(cfg_file)
  if (is.list(j_obj)) {
    pick_j <- function(k, default = "") {
      if (!is.null(j_obj[[k]]) && nzchar(as.character(j_obj[[k]]))) as.character(j_obj[[k]]) else default
    }
    raw_cfg$data_dir      <- pick_j("data_dir", "data")
    raw_cfg$gis_dir       <- pick_j("gis_dir", "gis")
    raw_cfg$output_dir    <- pick_j("output_dir", "output")
    raw_cfg$processed_dir <- pick_j("processed_dir", "")
  } else if (file.exists(cfg_file)) {
    txt <- tryCatch(paste(readLines(cfg_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n"), error = function(e) "")
    get_field <- function(k, default = "") {
      pat <- paste0('"', k, '"\\s*:\\s*"([^"]*)"')
      if (grepl(pat, txt)) sub(paste0('^.*?', pat, '.*$'), "\\1", txt) else default
    }
    raw_cfg$data_dir      <- get_field("data_dir", "data")
    raw_cfg$gis_dir       <- get_field("gis_dir", "gis")
    raw_cfg$output_dir    <- get_field("output_dir", "output")
    raw_cfg$processed_dir <- get_field("processed_dir", "")
  }

  resolve_p <- function(p, default) {
    p <- if (is.null(p) || is.na(p) || !nzchar(trimws(p))) default else trimws(p)
    p <- path.expand(p)
    if (!grepl("^(/|[A-Za-z]:)", p)) {
      p <- rel_path(base_dir, p)
    }
    sub("^\\./", "", p)
  }

  data_d <- resolve_p(raw_cfg$data_dir, "data")
  gis_d  <- resolve_p(raw_cfg$gis_dir, "gis")
  out_d  <- resolve_p(raw_cfg$output_dir, "output")
  proc_d <- if (nzchar(raw_cfg$processed_dir)) resolve_p(raw_cfg$processed_dir, rel_path(data_d, "processed")) else rel_path(data_d, "processed")

  if (isTRUE(auto_create) && nzchar(out_d)) {
    dir.create(out_d, recursive = TRUE, showWarnings = FALSE)
  }

  list(
    config_file   = if (file.exists(cfg_file)) cfg_file else NULL,
    data_dir      = data_d,
    gis_dir       = gis_d,
    output_dir    = out_d,
    processed_dir = proc_d
  )
}

parse_custom_dir_args <- function(args = character(0), base_dir = ".") {
  cfg_dirs <- load_directory_config(base_dir = base_dir, auto_create = TRUE)
  gis_dir  <- cfg_dirs$gis_dir
  data_dir <- cfg_dirs$data_dir
  out_dir  <- cfg_dirs$output_dir
  has_custom_out  <- FALSE
  has_custom_data <- FALSE
  proc_cli        <- ""

  env_gis  <- Sys.getenv("S3N_GIS_DIR", unset = "")
  env_data <- Sys.getenv("S3N_DATA_DIR", unset = "")
  env_out  <- Sys.getenv("S3N_OUTPUT_DIR", unset = "")
  env_proc <- Sys.getenv("S3N_PROCESSED_DIR", unset = "")
  if (nzchar(env_gis))  gis_dir  <- env_gis
  if (nzchar(env_data)) { data_dir <- env_data; has_custom_data <- TRUE }
  if (nzchar(env_out))  { out_dir <- env_out; has_custom_out <- TRUE }

  for (a in args) {
    if (grepl("^--gis-dir=", a)) {
      gis_dir <- sub("^--gis-dir=", "", a)
    } else if (grepl("^--data-dir=", a)) {
      data_dir <- sub("^--data-dir=", "", a)
      has_custom_data <- TRUE
    } else if (grepl("^--processed-dir=", a)) {
      proc_cli <- sub("^--processed-dir=", "", a)
    } else if (grepl("^--(out-dir|output-dir)=", a)) {
      out_dir <- sub("^--(out-dir|output-dir)=", "", a)
      has_custom_out <- TRUE
    }
  }

  resolve_arg <- function(p, default) {
    p <- if (is.null(p) || is.na(p) || !nzchar(trimws(p))) default else trimws(p)
    p <- path.expand(p)
    if (!is_abs_path(p)) {
      p <- rel_path(base_dir, p)
    }
    sub("^\\./", "", p)
  }

  gis_dir  <- resolve_arg(gis_dir, "gis")
  data_dir <- resolve_arg(data_dir, "data")
  out_dir  <- resolve_arg(out_dir, "output")
  if (nzchar(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  wq_dir <- rel_path(data_dir, "NCEA_River_Surveillance_Network_ARD_Open")
  proc_root <- if (nzchar(proc_cli)) {
    resolve_arg(proc_cli, rel_path(data_dir, "processed"))
  } else if (nzchar(env_proc)) {
    resolve_arg(env_proc, rel_path(data_dir, "processed"))
  } else if (has_custom_data && data_dir != cfg_dirs$data_dir) {
    rel_path(data_dir, "processed")
  } else {
    cfg_dirs$processed_dir
  }

  list(
    gis_dir        = gis_dir,
    data_dir       = data_dir,
    wq_dir         = wq_dir,
    proc_root      = proc_root,
    processed_dir  = proc_root,
    out_dir        = out_dir,
    has_custom_out = has_custom_out,
    config_file    = cfg_dirs$config_file
  )
}

# ==============================================================================
# CENTRALIZED DETERMINAND METADATA CATALOG (configs/determinands.json)
# ==============================================================================

.s3n_det_catalog_cache <- new.env(parent = emptyenv())

load_determinand_catalog <- function(base_dir = NULL, force_reload = FALSE) {
  if (is.null(base_dir) || !nzchar(base_dir)) base_dir <- resolve_base_dir()
  cache_key <- normalizePath(base_dir, mustWork = FALSE)
  if (!isTRUE(force_reload) && exists(cache_key, envir = .s3n_det_catalog_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .s3n_det_catalog_cache, inherits = FALSE))
  }

  cat_path <- rel_path(base_dir, "configs", "determinands.json")
  catalog <- list()

  j_obj <- read_json_safe(cat_path)
  if (is.list(j_obj) && length(j_obj) > 0) {
    for (k in names(j_obj)) {
      item <- j_obj[[k]]
      if (is.list(item)) {
        catalog[[k]] <- c(
          short = if (!is.null(item$short)) as.character(item$short) else k,
          desc  = if (!is.null(item$name))  as.character(item$name)  else k,
          name  = if (!is.null(item$name))  as.character(item$name)  else k,
          unit  = if (!is.null(item$unit))  as.character(item$unit)  else ""
        )
      }
    }
  } else if (file.exists(cat_path)) {
    txt <- paste(readLines(cat_path, warn = FALSE, encoding = "UTF-8"), collapse = " ")
    m <- gregexpr('"([^"]+)"\\s*:\\s*\\{([^{}]+)\\}', txt, perl = TRUE)
    blocks <- regmatches(txt, m)[[1]]
    for (blk in blocks) {
      k <- sub('^"([^"]+)".*$', "\\1", blk)
      get_f <- function(fk, def = "") {
        pat <- paste0('"', fk, '"\\s*:\\s*"([^"]*)"')
        if (grepl(pat, blk)) sub(paste0('^.*', pat, '.*$'), "\\1", blk) else def
      }
      nm <- get_f("name", k)
      sh <- get_f("short", nm)
      un <- get_f("unit", "")
      catalog[[k]] <- c(short = sh, desc = nm, name = nm, unit = un)
    }
  }

  assign(cache_key, catalog, envir = .s3n_det_catalog_cache)
  catalog
}

get_determinand_meta <- function(det_code, base_dir = NULL, default_unit = "") {
  catalog <- load_determinand_catalog(base_dir = base_dir)
  key <- as.character(det_code)[1]
  if (key %in% names(catalog)) return(catalog[[key]])
  lk <- tolower(key)
  hit <- which(tolower(names(catalog)) == lk)
  if (length(hit) > 0) return(catalog[[hit[1]]])
  c(short = key, desc = key, name = key, unit = default_unit)
}

# ==============================================================================
# DATA SOURCE MANAGEMENT (REGISTRY, MULTI-SOURCE & CUSTOM DATASETS)
# ==============================================================================

parse_source_config_file <- function(filepath) {
  if (!file.exists(filepath)) return(NULL)
  j_obj <- read_json_safe(filepath)
  if (is.list(j_obj)) {
    def_nm <- sub("\\.json$", "", basename(filepath))
    nm <- tolower(as.character(if (!is.null(j_obj$name)) j_obj$name else def_nm))
    return(list(
      name         = nm,
      abbr         = as.character(if (!is.null(j_obj$abbr)) j_obj$abbr else toupper(def_nm)),
      label        = as.character(if (!is.null(j_obj$label)) j_obj$label else nm),
      description  = as.character(if (!is.null(j_obj$description)) j_obj$description else ""),
      data_dir     = as.character(if (!is.null(j_obj$data_dir)) j_obj$data_dir else ""),
      sites_file   = as.character(if (!is.null(j_obj$sites_file)) j_obj$sites_file else ""),
      alias_target = tolower(as.character(if (!is.null(j_obj$alias_target)) j_obj$alias_target else "")),
      range_mode   = tolower(as.character(if (!is.null(j_obj$range_mode)) j_obj$range_mode else "")),
      sources      = if (!is.null(j_obj$sources)) tolower(as.character(unlist(j_obj$sources))) else character(0)
    ))
  }

  lines <- readLines(filepath, warn = FALSE)
  txt <- paste(lines, collapse = " ")
  get_field <- function(key, default = "") {
    pat <- paste0('"', key, '"\\s*:\\s*"([^"]*)"')
    if (grepl(pat, txt)) sub(paste0('^.*?', pat, '.*$'), "\\1", txt) else default
  }
  get_array <- function(key, default = character(0)) {
    pat <- paste0('"', key, '"\\s*:\\s*\\[([^]]*)\\]')
    if (grepl(pat, txt)) {
      m <- regmatches(txt, regexec(pat, txt))[[1]]
      if (length(m) > 1) {
        items <- unlist(regmatches(m[2], gregexpr('"([^"]*)"', m[2])))
        gsub('"', '', items)
      } else default
    } else default
  }
  list(
    name         = tolower(get_field("name", sub("\\.json$", "", basename(filepath)))),
    abbr         = get_field("abbr", toupper(sub("\\.json$", "", basename(filepath)))),
    label        = get_field("label", get_field("name")),
    description  = get_field("description", ""),
    data_dir     = get_field("data_dir", ""),
    sites_file   = get_field("sites_file", ""),
    alias_target = tolower(get_field("alias_target", "")),
    range_mode   = tolower(get_field("range_mode", "")),
    sources      = get_array("sources")
  )
}

.s3n_registry_cache <- new.env(parent = emptyenv())

load_source_registry <- function(base_dir = NULL, force_reload = FALSE) {
  if (is.null(base_dir) || !nzchar(base_dir)) base_dir <- resolve_base_dir()
  cache_key <- normalizePath(base_dir, mustWork = FALSE)
  if (!isTRUE(force_reload) && exists(cache_key, envir = .s3n_registry_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .s3n_registry_cache, inherits = FALSE))
  }
  cfg_dir <- rel_path(base_dir, "configs", "sources")
  sources <- list()
  composites <- list()

  sources$rsn <- list(
    name        = "rsn",
    abbr        = "RSN",
    label       = "River Surveillance Network (RSN)",
    description = "Environment Agency River Surveillance Network (1,382 stations)",
    data_dir    = "data/NCEA_River_Surveillance_Network_ARD_Open",
    sites_file  = "RSN Water Quality Sample Points.csv"
  )
  sources$ssn <- list(
    name        = "ssn",
    abbr        = "SSN",
    label       = "Small Streams Network (SSN)",
    description = "Environment Agency Small Streams Network (313 headwater stations)",
    data_dir    = "data/NCEA_Small_Streams_Network_ARD_Open",
    sites_file  = "SSN Water Quality Sample Points.csv"
  )

  composites$rss <- list(
    name        = "rss",
    abbr        = "RSS",
    label       = "Combined River Surveillance & Small Streams (RSS)",
    description = "Integrated River Surveillance Network and Small Streams Network",
    sources     = c("rsn", "ssn")
  )
  composites$all <- list(
    name        = "all",
    abbr        = "ALL",
    label       = "All Registered Monitoring Networks",
    description = "Union of all active monitoring networks",
    sources     = c("rsn", "ssn")
  )

  if (!dir.exists(cfg_dir) || length(list.files(cfg_dir, pattern = "\\.json$")) == 0) {
    py_script <- rel_path(base_dir, "scripts", "01_preprocess_orn.py")
    if (file.exists(py_script)) {
      suppressWarnings(tryCatch({
        run_python_script(py_script, "--list-sources", stdout = FALSE, stderr = FALSE)
      }, error = function(e) NULL))
    }
  }

  if (dir.exists(cfg_dir)) {
    j_files <- list.files(cfg_dir, pattern = "\\.json$", full.names = TRUE)
    for (jf in j_files) {
      fname <- basename(jf)
      if (fname == "composites.json") {
        j_comp <- read_json_safe(jf)
        if (is.list(j_comp) && is.list(j_comp$composites)) {
          for (cn in names(j_comp$composites)) {
            ci <- j_comp$composites[[cn]]
            if (is.list(ci) && length(ci$sources) > 0) {
              composites[[tolower(cn)]] <- list(
                name        = tolower(cn),
                abbr        = if (!is.null(ci$abbr)) as.character(ci$abbr) else toupper(cn),
                label       = if (!is.null(ci$label)) as.character(ci$label) else cn,
                description = if (!is.null(ci$description)) as.character(ci$description) else "",
                sources     = tolower(as.character(unlist(ci$sources)))
              )
            }
          }
          next
        }
        txt_c <- paste(readLines(jf, warn = FALSE), collapse = " ")
        m_blocks <- gregexpr('"([a-zA-Z0-9_]+)"\\s*:\\s*\\{[^{}]*"sources"\\s*:\\s*\\[[^\\]]*\\][^{}]*\\}', txt_c, perl = TRUE)
        blocks <- regmatches(txt_c, m_blocks)[[1]]
        for (blk in blocks) {
          c_name <- tolower(sub('^"([a-zA-Z0-9_]+)".*$', "\\1", blk))
          c_abbr <- if (grepl('"abbr"\\s*:\\s*"([^"]*)"', blk)) sub('^.*"abbr"\\s*:\\s*"([^"]*)".*$', "\\1", blk) else toupper(c_name)
          c_lbl  <- if (grepl('"label"\\s*:\\s*"([^"]*)"', blk)) sub('^.*"label"\\s*:\\s*"([^"]*)".*$', "\\1", blk) else c_name
          c_desc <- if (grepl('"description"\\s*:\\s*"([^"]*)"', blk)) sub('^.*"description"\\s*:\\s*"([^"]*)".*$', "\\1", blk) else ""
          s_inner <- sub('^.*"sources"\\s*:\\s*\\[([^\\]]*)\\].*$', "\\1", blk)
          s_parts <- trimws(gsub('"', "", unlist(strsplit(s_inner, ","))))
          s_parts <- tolower(s_parts[nzchar(s_parts)])
          if (nzchar(c_name) && length(s_parts) > 0) {
            composites[[c_name]] <- list(name = c_name, abbr = c_abbr, label = c_lbl, description = c_desc, sources = s_parts)
          }
        }
        next
      }
      parsed <- parse_source_config_file(jf)
      if (!is.null(parsed) && nzchar(parsed$name)) {
        sources[[parsed$name]] <- parsed
      }
    }
  }

  cfg_dirs  <- load_directory_config(base_dir = base_dir)
  cust_root <- rel_path(cfg_dirs$data_dir, "custom")
  if (dir.exists(cust_root)) {
    sub_dirs <- list.dirs(cust_root, recursive = FALSE, full.names = FALSE)
    for (sd in sub_dirs) {
      s_key <- tolower(sd)
      if (!s_key %in% names(sources)) {
        sources[[s_key]] <- list(
          name        = s_key,
          abbr        = toupper(s_key),
          label       = paste("Custom Dataset:", sd),
          description = paste("User-supplied custom dataset in", file.path(format_rel_display_path(cfg_dirs$data_dir, base_dir), "custom", sd)),
          data_dir    = file.path(format_rel_display_path(cfg_dirs$data_dir, base_dir), "custom", sd),
          sites_file  = "sites.csv"
        )
      }
    }
  }

  reg_out <- list(sources = sources, composites = composites)
  assign(cache_key, reg_out, envir = .s3n_registry_cache)
  reg_out
}

resolve_source_tag <- function(source = "rsn", base_dir = NULL) {
  if (is.null(source) || is.na(source) || !nzchar(trimws(as.character(source)))) return("rsn")
  if (is.null(base_dir) || !nzchar(base_dir)) base_dir <- resolve_base_dir()
  raw_orig <- trimws(as.character(source))
  if (!grepl("[,+]", raw_orig) && grepl("^[^()]*\\(([A-Za-z0-9_]+)\\)[^()]*$", raw_orig)) {
    inner <- tolower(sub("^[^()]*\\(([A-Za-z0-9_]+)\\)[^()]*$", "\\1", raw_orig))
    if (inner %in% c("rsn", "ssn", "rss")) return(inner)
    reg <- load_source_registry(base_dir)
    for (k in names(reg$sources)) {
      if (tolower(reg$sources[[k]]$abbr) == inner || k == inner) return(k)
    }
    for (k in names(reg$composites)) {
      if (tolower(reg$composites[[k]]$abbr) == inner || k == inner) return(k)
    }
  }
  raw_s <- tolower(raw_orig)
  tokens <- trimws(unlist(strsplit(raw_s, "[,+]+")))
  tokens <- tokens[nzchar(tokens)]
  reg <- load_source_registry(base_dir)

  norm_one <- function(tk) {
    tk_under <- gsub("-", "_", tk)
    cand <- if (tk %in% names(reg$sources) || tk %in% names(reg$composites) || dir.exists(file.path(base_dir, "data", "custom", tk))) {
      tk
    } else if (tk_under %in% names(reg$sources) || tk_under %in% names(reg$composites) || dir.exists(file.path(base_dir, "data", "custom", tk_under))) {
      tk_under
    } else {
      gsub("[^a-z0-9_]+", "_", tk_under)
    }
    if (cand %in% names(reg$sources) && !is.null(reg$sources[[cand]]$alias_target) && nzchar(reg$sources[[cand]]$alias_target)) {
      reg$sources[[cand]]$alias_target
    } else if (paste0(cand, "_tobit") %in% names(reg$sources)) {
      paste0(cand, "_tobit")
    } else {
      cand
    }
  }

  if (length(tokens) > 1) {
    norm_tokens <- unique(unlist(lapply(tokens, function(tk) {
      cnd <- norm_one(tk)
      if (cnd == "rss") c("rsn", "ssn")
      else if (cnd %in% names(reg$composites)) reg$composites[[cnd]]$sources
      else cnd
    })))
    norm_tokens <- sort(unique(norm_tokens))
    if (setequal(norm_tokens, c("rsn", "ssn"))) return("rss")
    return(paste(norm_tokens, collapse = "+"))
  }

  s <- norm_one(tokens[1])
  if (s %in% names(reg$composites)) s else s
}

format_source_abbr <- function(source = "rsn", base_dir = NULL) {
  s <- resolve_source_tag(source, base_dir = base_dir)
  if (s == "rsn") return("RSN")
  if (s == "ssn") return("SSN")
  if (s == "rss") return("RSS")

  reg <- load_source_registry(base_dir)
  if (s %in% names(reg$composites)) return(reg$composites[[s]]$abbr)
  if (s %in% names(reg$sources)) return(reg$sources[[s]]$abbr)

  s_tokens <- trimws(unlist(strsplit(s, "[,+]+")))
  s_tokens <- s_tokens[nzchar(s_tokens)]
  if (length(s_tokens) > 1) {
    tok_abbrs <- vapply(s_tokens, function(tk) {
      format_source_abbr(tk, base_dir = base_dir)
    }, character(1))
    return(paste(tok_abbrs, collapse = "+"))
  }

  if (grepl("_(tobit|geometric|midpoint|lower|upper)$", s)) {
    base_s <- sub("_(tobit|geometric|midpoint|lower|upper)$", "", s)
    rm     <- sub(".*_(tobit|geometric|midpoint|lower|upper)$", "\\1", s)
    rm_str <- tools::toTitleCase(rm)
    base_a <- if (base_s %in% names(reg$sources)) reg$sources[[base_s]]$abbr else toupper(base_s)
    if (grepl("\\[", base_a)) return(base_a)
    return(paste0(base_a, " [", rm_str, "]"))
  }

  toupper(s)
}

format_source_label <- function(source = "rsn", base_dir = NULL) {
  s <- resolve_source_tag(source, base_dir = base_dir)
  switch(s,
    rsn = "River Surveillance Network (RSN)",
    ssn = "Small Streams Network (SSN)",
    rss = "Combined River Surveillance & Small Streams (RSS)",
    {
      reg <- load_source_registry(base_dir)
      if (s %in% names(reg$composites)) return(reg$composites[[s]]$label)
      if (s %in% names(reg$sources)) return(reg$sources[[s]]$label)

      s_tokens <- trimws(unlist(strsplit(s, "[,+]+")))
      s_tokens <- s_tokens[nzchar(s_tokens)]
      if (length(s_tokens) > 1) {
        tok_labels <- vapply(s_tokens, function(tk) {
          format_source_label(tk, base_dir = base_dir)
        }, character(1))
        return(paste("Combined Monitoring Sources:", paste(tok_labels, collapse = " + ")))
      }

      if (grepl("_(tobit|geometric|midpoint|lower|upper)$", s)) {
        base_s <- sub("_(tobit|geometric|midpoint|lower|upper)$", "", s)
        rm     <- sub(".*_(tobit|geometric|midpoint|lower|upper)$", "\\1", s)
        rm_str <- tools::toTitleCase(rm)
        base_l <- if (base_s %in% names(reg$sources)) reg$sources[[base_s]]$label else tools::toTitleCase(gsub("_", " ", base_s))
        if (grepl("\\[", base_l)) return(base_l)
        return(paste0(base_l, " [", rm_str, "]"))
      }

      paste0("Custom Monitoring Source (", toupper(s), ")")
    }
  )
}

format_source_folder_tag <- function(source = "rsn") {
  resolve_source_tag(source)
}

parse_source_arg <- function(args) {
  src_arg <- grep("^--source=", args, value = TRUE)
  if (length(src_arg) > 0) {
    return(resolve_source_tag(sub("^--source=", "", src_arg[1])))
  }
  src_idx <- which(args == "--source")
  if (length(src_idx) > 0 && src_idx[1] < length(args)) {
    return(resolve_source_tag(args[src_idx[1] + 1]))
  }
  if (any(grepl("^--(custom-samples|custom-sites|custom-dir)=", args))) {
    cust_name_arg <- grep("^--custom-name=", args, value = TRUE)
    if (length(cust_name_arg) > 0) {
      return(resolve_source_tag(sub("^--custom-name=", "", cust_name_arg[1])))
    }
    return("custom")
  }
  "rsn"
}

resolve_data_source_dirs <- function(data_dir, source = "rsn", base_dir = NULL, custom_dir = NULL, custom_samples = NULL, custom_sites = NULL) {
  if (is.null(base_dir) || !nzchar(base_dir)) base_dir <- resolve_base_dir()
  src_tag  <- resolve_source_tag(source, base_dir = base_dir)
  reg      <- load_source_registry(base_dir = base_dir)

  rsn_dir  <- rel_path(data_dir, "NCEA_River_Surveillance_Network_ARD_Open")
  ssn_dir  <- rel_path(data_dir, "NCEA_Small_Streams_Network_ARD_Open")

  source_tokens <- if (src_tag == "rss") {
    c("rsn", "ssn")
  } else if (src_tag %in% names(reg$composites)) {
    reg$composites[[src_tag]]$sources
  } else if (grepl("[+]", src_tag)) {
    unlist(strsplit(src_tag, "[+]"))
  } else {
    src_tag
  }

  active_dirs <- list()
  for (tk in source_tokens) {
    if (tk == "rsn") {
      active_dirs[["rsn"]] <- rsn_dir
    } else if (tk == "ssn") {
      active_dirs[["ssn"]] <- ssn_dir
    } else if (!is.null(custom_dir) && nzchar(custom_dir) && dir.exists(custom_dir)) {
      active_dirs[[tk]] <- custom_dir
    } else if (tk %in% names(reg$sources)) {
      cfg_d <- reg$sources[[tk]]$data_dir
      cand_d <- if (is.null(cfg_d) || !nzchar(cfg_d)) {
        rel_path(data_dir, "custom", tk)
      } else if (dir.exists(cfg_d)) {
        cfg_d
      } else {
        rel_path(base_dir, cfg_d)
      }
      active_dirs[[tk]] <- if (dir.exists(cand_d)) cand_d else rel_path(data_dir, "custom", tk)
    } else {
      active_dirs[[tk]] <- rel_path(data_dir, "custom", tk)
    }
  }

  list(
    source         = src_tag,
    source_abbr    = format_source_abbr(src_tag, base_dir = base_dir),
    source_label   = format_source_label(src_tag, base_dir = base_dir),
    folder_tag     = format_source_folder_tag(src_tag),
    rsn_dir        = rsn_dir,
    ssn_dir        = ssn_dir,
    active_dirs    = active_dirs,
    custom_samples = custom_samples,
    custom_sites   = custom_sites
  )
}

apply_range_mode_to_source <- function(source_tag, range_mode = "auto", base_dir = NULL) {
  rm <- tolower(trimws(as.character(range_mode)))
  if (!nzchar(rm) || rm %in% c("auto", "none")) {
    return(list(source_tag = source_tag, range_mode = rm))
  }
  valid_modes <- c("tobit", "geometric", "midpoint", "lower", "upper")
  if (!(rm %in% valid_modes)) return(list(source_tag = source_tag, range_mode = rm))
  if (grepl(paste0("_", rm, "$"), source_tag)) return(list(source_tag = source_tag, range_mode = rm))
  if (grepl("[+]", source_tag) && any(endsWith(strsplit(source_tag, "[+]")[[1]], paste0("_", rm)))) {
    return(list(source_tag = source_tag, range_mode = rm))
  }
  if (source_tag == "rsn" || source_tag == "ssn" || source_tag == "rss") {
    return(list(source_tag = source_tag, range_mode = rm))
  }
  base_cand <- sub("_(tobit|geometric|midpoint|lower|upper)$", "", source_tag)
  eff_tag <- resolve_source_tag(paste0(base_cand, "_", rm), base_dir = base_dir)
  list(source_tag = eff_tag, range_mode = rm)
}

print_available_sources <- function(base_dir = NULL) {
  if (is.null(base_dir) || !nzchar(base_dir)) base_dir <- resolve_base_dir()
  reg <- load_source_registry(base_dir = base_dir)

  cat("Available Monitoring Data Sources:\n\n")
  cat(sprintf("  %-18s | %-6s | %-44s | %-12s\n", "Source", "Abbr", "Description", "Type"))
  cat(strrep("-", 88), "\n")
  cat(sprintf("  %-18s | %-6s | %-44s | %-12s\n", "rsn", "RSN", "River Surveillance Network (1,382 stations)", "YES (Default)"))
  cat(sprintf("  %-18s | %-6s | %-44s | %-12s\n", "ssn", "SSN", "Small Streams Network (313 headwaters)", "Registered"))
  cat(sprintf("  %-18s | %-6s | %-44s | %-12s\n", "rss", "RSS", "Combined River Surveillance & Small Streams", "Composite"))

  for (s_name in names(reg$sources)) {
    if (s_name %in% c("rsn", "ssn")) next
    s_info <- reg$sources[[s_name]]
    desc <- if (nzchar(s_info$description)) s_info$description else s_info$label
    cat(sprintf("  %-18s | %-6s | %-44s | %-12s\n",
                substr(s_name, 1, 18), substr(s_info$abbr, 1, 6), substr(desc, 1, 44), "Configured"))
  }

  cat("\nUsage & Multi-Source Syntax:\n")
  cat("  Rscript scripts/run_s3n_determinand.R 0117                  # Default RSN\n")
  cat("  Rscript scripts/run_s3n_determinand.R 0117 --source=ssn      # Small Streams Network only\n")
  cat("  Rscript scripts/run_s3n_determinand.R 0117 --source=rss      # Combined RSN + SSN\n")
  cat("  Rscript scripts/run_s3n_determinand.R 0117 --source=rsn+ssn  # On-the-fly combination syntax (+ or ,)\n")
  cat("  Rscript scripts/run_s3n_determinand.R 0117 --source=custom_example # Registered custom dataset\n")
  cat("  Rscript scripts/run_s3n_determinand.R 0117 --custom-samples=data/custom/example_study/samples.csv\n\n")
}

# ==============================================================================
# WATER NETWORK METADATA & DIRECTORY RESOLUTION HELPERS
# ==============================================================================

resolve_network_processed_dir <- function(base_dir, network_arg = NULL, data_dir = NULL, gis_dir = NULL, source = "rsn", processed_dir = NULL) {
  src_tag   <- resolve_source_tag(source, base_dir = base_dir)
  src_suff  <- paste0("_", src_tag)
  cfg_dirs  <- load_directory_config(base_dir = base_dir)
  d_dir     <- if (!is.null(data_dir) && nzchar(data_dir)) data_dir else cfg_dirs$data_dir
  g_dir     <- if (!is.null(gis_dir) && nzchar(gis_dir)) gis_dir else cfg_dirs$gis_dir
  proc_root <- if (!is.null(processed_dir) && nzchar(processed_dir)) {
    processed_dir
  } else if (nzchar(cfg_dirs$processed_dir) && (is.null(data_dir) || data_dir == cfg_dirs$data_dir)) {
    cfg_dirs$processed_dir
  } else {
    rel_path(d_dir, "processed")
  }

  if (is.null(network_arg) || !nzchar(trimws(network_arg))) {
    abbr <- "orn"
    net_cli_arg <- "orn"
  } else {
    raw_net_orig <- trimws(network_arg)
    raw_net <- gsub("-", "_", tolower(raw_net_orig))
    if (grepl("\\.(gpkg|parquet|geoparquet|pq|json)$", raw_net)) {
      stem <- sub("\\.(gpkg|parquet|geoparquet|pq|json)$", "", basename(raw_net))
      stem <- sub("\\.geo$", "", stem)
      abbr <- tolower(stem)
      net_cli_arg <- raw_net_orig
    } else {
      abbr <- raw_net
      net_cli_arg <- raw_net_orig
    }
  }

  sub_d <- rel_path(proc_root, paste0(abbr, src_suff))

  if (!file.exists(rel_path(sub_d, "england_rsn_reaches.csv"))) {
    dirs_info <- resolve_data_source_dirs(d_dir, source = source, base_dir = base_dir)
    py_src <- if (length(dirs_info$active_dirs) > 1) {
      paste(names(dirs_info$active_dirs), collapse = "+")
    } else {
      src_tag
    }
    cat(sprintf("Preprocessing water network '%s' (source='%s') via scripts/01_preprocess_orn.py...\n", net_cli_arg, py_src))
    py_args <- c(paste0("--network=", shQuote(net_cli_arg)),
                 paste0("--source=", shQuote(py_src)),
                 paste0("--out-dir=", shQuote(sub_d)))
    if (!is.null(g_dir) && nzchar(g_dir)) {
      py_args <- c(py_args, paste0("--gis-dir=", shQuote(g_dir)))
    }
    if (!is.null(d_dir) && nzchar(d_dir)) {
      py_args <- c(py_args, paste0("--data-dir=", shQuote(d_dir)))
    }
    status <- run_python_script(rel_path(base_dir, "scripts", "01_preprocess_orn.py"), py_args)
    if (status != 0L) {
      stop(sprintf("Failed to preprocess water network '%s' with source '%s'. Run `python3 scripts/01_preprocess_orn.py --list-networks`.", net_cli_arg, src_tag))
    }
  }
  sub_d
}

load_network_metadata <- function(proc_dir) {
  meta_p <- file.path(proc_dir, "network_metadata.json")
  if (!file.exists(meta_p)) {
    return(list(
      abbr         = "orn",
      label        = "OS Open Rivers Network (ORN v2)",
      short_label  = "ORN v2",
      source       = "rsn",
      source_abbr  = "RSN",
      source_label = "River Surveillance Network (RSN)",
      range_mode   = "",
      meas_method  = "",
      humber_cid   = 80L,
      thames_cid   = 3055L,
      severn_cid   = 2786L
    ))
  }
  j_obj <- read_json_safe(meta_p)
  if (is.list(j_obj)) {
    src_key <- as.character(if (!is.null(j_obj$source)) j_obj$source else "rsn")
    cb <- if (is.list(j_obj$core_basins)) j_obj$core_basins else list()
    return(list(
      abbr         = as.character(if (!is.null(j_obj$abbr)) j_obj$abbr else "orn"),
      label        = as.character(if (!is.null(j_obj$label)) j_obj$label else "OS Open Rivers Network (ORN v2)"),
      short_label  = as.character(if (!is.null(j_obj$short_label)) j_obj$short_label else "ORN v2"),
      source       = src_key,
      source_abbr  = as.character(if (!is.null(j_obj$source_abbr)) j_obj$source_abbr else format_source_abbr(src_key)),
      source_label = as.character(if (!is.null(j_obj$source_label)) j_obj$source_label else format_source_label(src_key)),
      range_mode   = as.character(if (!is.null(j_obj$range_mode)) j_obj$range_mode else ""),
      meas_method  = as.character(if (!is.null(j_obj$meas_method)) j_obj$meas_method else ""),
      humber_cid   = if (!is.null(cb$humber)) as.integer(cb$humber) else NA_integer_,
      thames_cid   = if (!is.null(cb$thames)) as.integer(cb$thames) else NA_integer_,
      severn_cid   = if (!is.null(cb$severn)) as.integer(cb$severn) else NA_integer_
    ))
  }
  lines <- readLines(meta_p, warn = FALSE)
  src_key <- parse_json_string_field(lines, "source", "rsn")
  list(
    abbr         = parse_json_string_field(lines, "abbr", "orn"),
    label        = parse_json_string_field(lines, "label", "OS Open Rivers Network (ORN v2)"),
    short_label  = parse_json_string_field(lines, "short_label", "ORN v2"),
    source       = src_key,
    source_abbr  = parse_json_string_field(lines, "source_abbr", format_source_abbr(src_key)),
    source_label = parse_json_string_field(lines, "source_label", format_source_label(src_key)),
    range_mode   = parse_json_string_field(lines, "range_mode", ""),
    meas_method  = parse_json_string_field(lines, "meas_method", ""),
    humber_cid   = parse_json_int_field(lines, "humber", NA_integer_),
    thames_cid   = parse_json_int_field(lines, "thames", NA_integer_),
    severn_cid   = parse_json_int_field(lines, "severn", NA_integer_)
  )
}

extract_network_label_from_tables <- function(cv_tab = NULL, param_tab = NULL, dir_path = NULL) {
  if (!is.null(cv_tab) && "Network_Label" %in% colnames(cv_tab) && nrow(cv_tab) > 0) {
    val <- as.character(cv_tab$Network_Label[1])
    if (!is.na(val) && nzchar(val)) return(val)
  }
  if (!is.null(param_tab) && "Network_Label" %in% colnames(param_tab) && nrow(param_tab) > 0) {
    val <- as.character(param_tab$Network_Label[1])
    if (!is.na(val) && nzchar(val)) return(val)
  }
  if (!is.null(dir_path)) {
    bname <- basename(normalizePath(dir_path, mustWork = FALSE))
    if (grepl("euhydro_thames", bname)) return("EU-Hydro (Thames)")
    if (grepl("euhydro_tweed", bname))  return("EU-Hydro (Tweed)")
    if (grepl("euhydro", bname))        return("EU-Hydro v1.3")
    if (grepl("_orn", bname))           return("ORN v2")
  }
  "ORN v2"
}

# ==============================================================================
# MULTI-SOURCE PREFERENCE, RELIABILITY WEIGHTING & UNCERTAINTY HELPERS
# ==============================================================================

#' Parse --prefer-source=<s1[,s2,...]> CLI flag into ordered lowercase source tokens.
parse_prefer_source_arg <- function(args) {
  hit <- grep("^--prefer-source=", args, value = TRUE)
  if (length(hit) == 0L) return(character(0))
  raw <- sub("^--prefer-source=", "", hit[length(hit)])
  toks <- tolower(trimws(unlist(strsplit(raw, "[,+]+"))))
  unique(toks[nzchar(toks)])
}

#' Parse --source-weight=rsn:1,wb:0.25 or --source-se=rsn:0.01,wb:0.5 CLI flags.
parse_source_kv_arg <- function(args, flag_prefix = "--source-weight=",
                                min_val = 1e-4, normalize_gt1 = FALSE) {
  hits <- grep(paste0("^", flag_prefix), args, value = TRUE)
  if (length(hits) == 0L) return(numeric(0))
  out <- numeric(0)
  for (h in hits) {
    raw <- sub(paste0("^", flag_prefix), "", h)
    pairs <- trimws(unlist(strsplit(raw, "[,;]+")))
    pairs <- pairs[nzchar(pairs)]
    for (pr in pairs) {
      parts <- trimws(unlist(strsplit(pr, "[:=]")))
      if (length(parts) != 2L || !nzchar(parts[1]) || !nzchar(parts[2])) {
        stop(sprintf("Invalid entry '%s' in %s (expected <source>:<numeric_value>).", pr, flag_prefix))
      }
      key <- tolower(parts[1])
      val <- suppressWarnings(as.numeric(parts[2]))
      if (!is.finite(val) || val < min_val) {
        stop(sprintf("Invalid numeric value '%s' for source '%s' in %s (must be >= %g).",
                     parts[2], key, flag_prefix, min_val))
      }
      out[key] <- val
    }
  }
  if (normalize_gt1 && length(out) > 0L && max(out) > 1.0) {
    out <- out / max(out)
  }
  out
}

#' Vectorized source lookup matching exact lowercase name first, then base name without range suffix.
lookup_source_numeric <- function(src_vec, kv_map, default_val = NA_real_) {
  n <- length(src_vec)
  if (n == 0L) return(numeric(0))
  if (length(kv_map) == 0L) return(rep(as.numeric(default_val), n))
  s_low  <- tolower(trimws(as.character(src_vec)))
  s_base <- sub("_(tobit|geometric|midpoint|lower|upper)$", "", s_low)
  res    <- rep(as.numeric(default_val), n)
  nms    <- names(kv_map)
  for (i in seq_len(n)) {
    if (s_low[i] %in% nms) {
      res[i] <- unname(kv_map[s_low[i]])
    } else if (s_base[i] %in% nms) {
      res[i] <- unname(kv_map[s_base[i]])
    }
  }
  res
}

#' Filter co-located samples on identical RivIDs according to --prefer-source priority order.
filter_samples_by_source_priority <- function(df, prefer_sources, reach_col = "RivID", source_col = "data_source") {
  if (length(prefer_sources) == 0L || nrow(df) == 0L || !(reach_col %in% names(df))) {
    return(list(df = df, n_dropped = 0L, n_reaches = 0L))
  }
  src_vec <- if (source_col %in% names(df)) as.character(df[[source_col]]) else rep("unknown", nrow(df))
  rank_map <- stats::setNames(seq_along(prefer_sources), tolower(prefer_sources))
  ranks    <- lookup_source_numeric(src_vec, rank_map, default_val = length(prefer_sources) + 1.0)

  reach_idx <- split(seq_len(nrow(df)), df[[reach_col]])
  keep <- rep(TRUE, nrow(df))
  n_affected <- 0L
  for (ix in reach_idx) {
    r_ix  <- ranks[ix]
    r_min <- min(r_ix)
    drop_local <- r_ix > r_min
    if (any(drop_local)) {
      keep[ix[drop_local]] <- FALSE
      n_affected <- n_affected + 1L
    }
  }
  list(
    df        = df[keep, , drop = FALSE],
    n_dropped = as.integer(sum(!keep)),
    n_reaches = as.integer(n_affected)
  )
}

#' Aggregate sample rows to reach level, supporting --source-weight and --source-se.
aggregate_samples_by_reach_weighted <- function(df, trans_obj, src_weights = numeric(0), src_se = numeric(0)) {
  if (!"data_source" %in% names(df)) df$data_source <- "RSN"
  if (!"meas_uncertainty" %in% names(df)) df$meas_uncertainty <- NA_real_
  if (!"meas_method" %in% names(df)) df$meas_method <- ""
  if (!"is_cens" %in% names(df)) df$is_cens <- FALSE
  if (!"raw_lod" %in% names(df)) df$raw_lod <- df$val

  if (length(src_se) > 0L) {
    se_ov <- lookup_source_numeric(df$data_source, src_se, default_val = NA_real_)
    has_ov <- !is.na(se_ov) & is.finite(se_ov) & se_ov >= 0
    df$meas_uncertainty[has_ov] <- se_ov[has_ov]
  }

  has_sw <- (length(src_weights) > 0L)
  riv_groups <- split(seq_len(nrow(df)), df$RivID)
  n_r <- length(riv_groups)

  riv_ids <- as.integer(names(riv_groups))
  val_out <- numeric(n_r)
  y_out   <- numeric(n_r)
  up_out  <- numeric(n_r)
  cens_o  <- integer(n_r)
  lod_out <- rep(NA_real_, n_r)
  unc_out <- rep(NA_real_, n_r)
  mth_out <- character(n_r)
  wt_out  <- rep(1.0, n_r)

  for (k in seq_len(n_r)) {
    ix <- riv_groups[[k]]
    if (!has_sw) {
      val_out[k] <- mean(df$val[ix])
      y_out[k]   <- mean(df$y_sample[ix])
      up_out[k]  <- mean(df$obs_up_dist_km[ix])
      all_c      <- all(df$is_cens[ix] == 1L)
      cens_o[k]  <- as.integer(all_c)
      if (all_c) lod_out[k] <- mean(trans_obj$fwd(df$raw_lod[ix]), na.rm = TRUE)
      u <- df$meas_uncertainty[ix]
      u_clean <- u[!is.na(u) & is.finite(u) & u >= 0]
      if (length(u_clean) > 0L) unc_out[k] <- sqrt(sum(u_clean^2)) / length(u_clean)
      m <- df$meas_method[ix]
      m_clean <- m[!is.na(m) & nzchar(m)]
      if (length(m_clean) > 0L) mth_out[k] <- m_clean[1]
      wt_out[k] <- 1.0
    } else {
      # Group by data_source first so high sampling frequency within one campaign never outvotes another network
      src_splits <- split(ix, tolower(trimws(as.character(df$data_source[ix]))))
      ns <- length(src_splits)
      s_val  <- numeric(ns); s_y <- numeric(ns); s_up <- numeric(ns)
      s_w    <- numeric(ns); s_u <- rep(NA_real_, ns); s_c <- logical(ns)
      s_lod  <- rep(NA_real_, ns)
      for (j in seq_len(ns)) {
        sj <- src_splits[[j]]
        s_val[j] <- mean(df$val[sj])
        s_y[j]   <- mean(df$y_sample[sj])
        s_up[j]  <- mean(df$obs_up_dist_km[sj])
        s_w[j]   <- lookup_source_numeric(df$data_source[sj[1]], src_weights, default_val = 1.0)
        s_c[j]   <- all(df$is_cens[sj] == 1L)
        if (s_c[j]) s_lod[j] <- mean(trans_obj$fwd(df$raw_lod[sj]), na.rm = TRUE)
        u_j <- df$meas_uncertainty[sj]
        u_clean <- u_j[!is.na(u_j) & is.finite(u_j) & u_j >= 0]
        if (length(u_clean) > 0L) s_u[j] <- sqrt(sum(u_clean^2)) / length(u_clean)
      }
      w_norm <- s_w / sum(s_w)
      val_out[k] <- sum(w_norm * s_val)
      y_out[k]   <- sum(w_norm * s_y)
      up_out[k]  <- sum(w_norm * s_up)
      all_c      <- all(s_c)
      cens_o[k]  <- as.integer(all_c)
      if (all_c) lod_out[k] <- sum(w_norm * s_lod, na.rm = TRUE)
      ok_u <- !is.na(s_u) & is.finite(s_u)
      if (any(ok_u)) {
        w_u <- s_w[ok_u] / sum(s_w[ok_u])
        unc_out[k] <- sqrt(sum((w_u * s_u[ok_u])^2))
      }
      m <- df$meas_method[ix]
      m_clean <- m[!is.na(m) & nzchar(m)]
      if (length(m_clean) > 0L) mth_out[k] <- m_clean[1]
      wt_out[k] <- max(1e-4, min(1.0, max(s_w)))
    }
  }

  data.frame(
    RivID            = riv_ids,
    val              = val_out,
    y_resp           = y_out,
    obs_up_dist_km   = up_out,
    is_censored      = cens_o,
    lod_sc           = lod_out,
    meas_uncertainty = unc_out,
    meas_method      = mth_out,
    source_weight    = wt_out,
    stringsAsFactors = FALSE
  )
}
