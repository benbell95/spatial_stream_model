(function () {
  "use strict";

  function initThemeToggle() {
    const btn = document.getElementById("themeToggleBtn");
    const root = document.documentElement;
    if (!btn) return;
    function apply(theme) {
      root.setAttribute("data-theme", theme);
      btn.textContent = theme === "dark" ? "☀️ Light Mode" : "🌙 Dark Mode";
      try { localStorage.setItem("modelling_dashboard_theme", theme); } catch (e) {}
    }
    let saved = "dark";
    try { saved = localStorage.getItem("modelling_dashboard_theme") || "dark"; } catch (e) {}
    apply(saved);
    btn.addEventListener("click", () => apply(root.getAttribute("data-theme") === "dark" ? "light" : "dark"));
  }

  function fallbackCopy(txt, cb) {
    const ta = document.createElement("textarea");
    ta.value = txt;
    ta.style.position = "fixed";
    ta.style.opacity = "0";
    document.body.appendChild(ta);
    ta.select();
    try { document.execCommand("copy"); cb(); } catch (e) {}
    document.body.removeChild(ta);
  }

  function initCopyButtons() {
    document.querySelectorAll("[data-copy-target]").forEach((btn) => {
      btn.addEventListener("click", () => {
        const target = document.getElementById(btn.dataset.copyTarget);
        if (!target) return;
        const txt = target.textContent;
        const done = () => {
          const orig = btn.textContent;
          btn.textContent = "Copied!";
          setTimeout(() => { btn.textContent = orig; }, 1400);
          if (btn.dataset.copyTarget === "wbOutputCmd" && typeof window.__wbRecordHistory === "function") {
            window.__wbRecordHistory(txt);
          }
        };
        if (navigator.clipboard && navigator.clipboard.writeText) {
          navigator.clipboard.writeText(txt).then(done).catch(() => fallbackCopy(txt, done));
        } else {
          fallbackCopy(txt, done);
        }
      });
    });
  }

  function initInfoTooltips() {
    function closeAll(exceptWrap) {
      document.querySelectorAll(".info-wrap.open").forEach((wrap) => {
        if (wrap !== exceptWrap) wrap.classList.remove("open");
      });
    }

    document.querySelectorAll(".info-tip-btn").forEach((btn) => {
      btn.setAttribute("aria-expanded", "false");
      btn.addEventListener("click", (e) => {
        e.preventDefault();
        e.stopPropagation();
        const wrap = btn.closest(".info-wrap");
        if (!wrap) return;
        const willOpen = !wrap.classList.contains("open");
        closeAll(wrap);
        wrap.classList.toggle("open", willOpen);
        btn.setAttribute("aria-expanded", willOpen ? "true" : "false");
      });
    });

    document.addEventListener("click", (e) => {
      if (!e.target.closest(".info-wrap")) {
        closeAll(null);
        document.querySelectorAll(".info-tip-btn").forEach((b) => b.setAttribute("aria-expanded", "false"));
      }
    });

    document.addEventListener("keydown", (e) => {
      if (e.key === "Escape") {
        closeAll(null);
        document.querySelectorAll(".info-tip-btn").forEach((b) => b.setAttribute("aria-expanded", "false"));
      }
    });
  }

  function initWorkbench() {
    const ids = [
      "wbMode", "wbTargetSel", "wbTargetCustom", "wbDataset",
      "wbSourceSel", "wbSourceCustom", "wbNetworkSel", "wbNetworkCustom",
      "wbTransform", "wbNorm", "wbCensor", "wbTrend", "wbKernel", "wbAfv",
      "wbCombine", "wbVarsSel", "wbVarsCustom", "wbBoot",
      "wbCvMode", "wbThresholdSel", "wbThresholdVal", "wbThresholdDir",
      "wbTemporal", "wbTemporalSeason", "wbTemporalYear", "wbTemporalMonth",
      "wbDesignSites", "wbBasinSel", "wbBasinCustom",
      "wbPreferSel", "wbPreferCustom", "wbWeightSel", "wbWeightCustom",
      "wbSeSel", "wbSeCustom", "wbFormat", "wbPalette", "wbOmp",
      "wbRVersionSel", "wbRVersionCustom",
      "wbPostDir1", "wbCompareMode", "wbPostDir2", "wbPostLabels", "wbPostOut",
      "wbSnapCatchmentSel", "wbSnapCatchmentCustom", "wbSnapOut",
      "wbPlotDir", "wbPlotBasin", "wbPlotVmin", "wbPlotVmax", "wbPlotFigSize", "wbPlotNorm", "wbPlotPalette", "wbPlotFormat"
    ];
    const els = {};
    ids.forEach((id) => { els[id] = document.getElementById(id); });

    const outEl = document.getElementById("wbOutputCmd");
    const targetScriptEl = document.getElementById("wbCmdTargetScript");
    const notesEl = document.getElementById("wbNotesList");
    const wrapBtn = document.getElementById("wbWrapToggleBtn");
    const resetBtn = document.getElementById("wbResetBtn");
    const thresholdWrap = document.getElementById("wbThresholdCustomWrap");
    const winRBar = document.getElementById("wbWinRBar");
    const winFullPathEl = document.getElementById("wbWinFullPath");
    const winFullPathLabel = document.getElementById("wbWinFullPathLabel");
    const winRVersionWrap = document.getElementById("wbWinRVersionWrap");
    const postCardEl = document.getElementById("wbPostCard");
    const snapCardEl = document.getElementById("wbSnapCard");
    const plotCardEl = document.getElementById("wbPlotCard");
    const postDir2Wrap = document.getElementById("wbPostDir2Wrap");
    const activateSnapBtn = document.getElementById("wbActivateSnapBtn");
    const activatePlotBtn = document.getElementById("wbActivatePlotBtn");
    const discoveryCardEl = document.getElementById("wbDiscoveryCard");
    const configureSectionEl = document.getElementById("wbConfigureSection");
    if (!outEl) return;

    const platformHint = (typeof navigator !== "undefined" && ((navigator.userAgentData && navigator.userAgentData.platform) || navigator.platform)) || "";
    const uaHint = (typeof navigator !== "undefined" && navigator.userAgent) || "";
    const isMacOS = /mac/i.test(platformHint) || (/Mac OS X/i.test(uaHint) && !/iPhone|iPad|iPod/i.test(uaHint));
    const isWindows = /win/i.test(platformHint) || /Windows/i.test(uaHint);
    const unixShellLabel = isMacOS ? "Zsh" : "Bash";
    const unixShellBtn = document.querySelector('#wbShellBtns [data-shell="bash"]');
    if (unixShellBtn) unixShellBtn.textContent = unixShellLabel;

    let wrapFlags = false;
    let shellMode = isWindows ? "powershell" : "bash";
    let activeListCmd = "";
    let activeSectionTab = "model";

    const DEFAULTS = {
      wbMode: "det",
      wbTargetSel: "0117",
      wbTargetCustom: "0117",
      wbDataset: "auto",
      wbSourceSel: "rsn",
      wbSourceCustom: "rsn",
      wbNetworkSel: "orn",
      wbNetworkCustom: "",
      wbTransform: "default",
      wbNorm: "none",
      wbCensor: "default",
      wbTrend: "linear",
      wbKernel: "exponential",
      wbAfv: "length",
      wbCombine: "joint",
      wbVarsSel: "",
      wbVarsCustom: "",
      wbBoot: "50",
      wbCvMode: "random",
      wbThresholdSel: "none",
      wbThresholdVal: "10",
      wbThresholdDir: "above",
      wbTemporal: "all",
      wbTemporalSeason: "winter",
      wbTemporalYear: "2024",
      wbTemporalMonth: "1",
      wbDesignSites: "0",
      wbBasinSel: "none",
      wbBasinCustom: "",
      wbPreferSel: "",
      wbPreferCustom: "",
      wbWeightSel: "",
      wbWeightCustom: "",
      wbSeSel: "",
      wbSeCustom: "",
      wbFormat: "png_200",
      wbPalette: "standard",
      wbOmp: "1",
      wbRVersionSel: "4.4.3",
      wbRVersionCustom: "",
      wbPostDir1: "output/results/det_0117_orn_rsn",
      wbCompareMode: "dirs",
      wbPostDir2: "output/results/det_0180_orn_rsn",
      wbPostLabels: "",
      wbPostOut: "",
      wbSnapCatchmentSel: "all",
      wbSnapCatchmentCustom: "",
      wbSnapOut: "",
      wbPlotDir: "output/results/det_0117_orn_rsn",
      wbPlotBasin: "",
      wbPlotVmin: "",
      wbPlotVmax: "",
      wbPlotFigSize: "",
      wbPlotNorm: "none",
      wbPlotPalette: "auto",
      wbPlotFormat: "png_200"
    };

    function resolveSelectOrCustom(selEl, customEl, fallback) {
      if (!selEl) return fallback || "";
      if (selEl.value === "custom") {
        return (customEl ? customEl.value.trim() : "") || fallback || "";
      }
      return selEl.value;
    }

    function syncCustomVisibility() {
      const pairs = [
        [els.wbTargetSel, els.wbTargetCustom],
        [els.wbSourceSel, els.wbSourceCustom],
        [els.wbNetworkSel, els.wbNetworkCustom],
        [els.wbBasinSel, els.wbBasinCustom],
        [els.wbPreferSel, els.wbPreferCustom],
        [els.wbWeightSel, els.wbWeightCustom],
        [els.wbSeSel, els.wbSeCustom],
        [els.wbRVersionSel, els.wbRVersionCustom],
        [els.wbSnapCatchmentSel, els.wbSnapCatchmentCustom],
        [els.wbVarsSel, els.wbVarsCustom]
      ];
      pairs.forEach(([sel, inp]) => {
        if (!sel || !inp) return;
        const show = sel.value === "custom";
        inp.hidden = !show;
      });
      if (els.wbThresholdSel && thresholdWrap) {
        thresholdWrap.hidden = els.wbThresholdSel.value !== "custom";
      }
      if (els.wbTemporal) {
        if (els.wbTemporalSeason) els.wbTemporalSeason.hidden = els.wbTemporal.value !== "custom_season";
        if (els.wbTemporalYear) els.wbTemporalYear.hidden = els.wbTemporal.value !== "custom_year";
        if (els.wbTemporalMonth) els.wbTemporalMonth.hidden = els.wbTemporal.value !== "custom_month";
      }
      if (postDir2Wrap && els.wbCompareMode) {
        const isDirs = els.wbCompareMode.value === "dirs";
        postDir2Wrap.classList.toggle("is-disabled", !isDirs);
        if (els.wbPostDir2) els.wbPostDir2.disabled = !isDirs;
        if (postDir2AddBtn) postDir2AddBtn.disabled = !isDirs;
      }
    }

    function syncWindowsROptions(mode) {
      const showWinBar = shellMode === "powershell" || shellMode === "cmd";
      if (winRBar) winRBar.hidden = !showWinBar;
      const isPrep = mode === "prep";
      const useFull = Boolean(winFullPathEl && winFullPathEl.checked);
      if (winFullPathEl) winFullPathEl.disabled = !showWinBar || isPrep;
      if (winFullPathLabel) winFullPathLabel.classList.toggle("is-disabled", isPrep);
      const enableVer = showWinBar && !isPrep && useFull;
      if (winRVersionWrap) winRVersionWrap.classList.toggle("is-disabled", !enableVer);
      if (els.wbRVersionSel) els.wbRVersionSel.disabled = !enableVer;
      if (els.wbRVersionCustom) els.wbRVersionCustom.disabled = !enableVer;
    }

    function syncFieldApplicability(mode) {
      document.querySelectorAll(".wb-sections-grid .field[data-applies]").forEach((field) => {
        const applies = (field.getAttribute("data-applies") || "all").split(",").map((s) => s.trim());
        const active = applies.includes("all") || applies.includes(mode);
        field.classList.toggle("is-disabled", !active);
        field.querySelectorAll("select, input").forEach((ctrl) => {
          ctrl.disabled = !active;
        });
      });

      if (postCardEl) {
        postCardEl.classList.toggle("is-mode-active", mode === "stats" || mode === "compare");
      }
      if (snapCardEl) {
        snapCardEl.classList.toggle("is-mode-active", mode === "snap");
      }
      if (plotCardEl) {
        plotCardEl.classList.toggle("is-mode-active", mode === "plots");
      }
      document.querySelectorAll("#wbPostActionBtns [data-post-mode]").forEach((btn) => {
        btn.classList.toggle("active", btn.dataset.postMode === mode);
      });
      if (activateSnapBtn) {
        activateSnapBtn.classList.toggle("active", mode === "snap");
      }
      if (activatePlotBtn) {
        activatePlotBtn.classList.toggle("active", mode === "plots");
      }

      // Keep Trend / Co-Kriging options coherent for active mode
      if (els.wbTrend) {
        Array.from(els.wbTrend.options).forEach((opt) => {
          const isCokrigMethod = opt.value === "lmc" || opt.value === "colocated";
          if (mode === "cokrig") {
            opt.hidden = !isCokrigMethod;
          } else {
            opt.hidden = isCokrigMethod;
          }
        });
        if (mode === "cokrig" && !["lmc", "colocated"].includes(els.wbTrend.value)) {
          els.wbTrend.value = "lmc";
        } else if (mode === "det" && ["lmc", "colocated"].includes(els.wbTrend.value)) {
          els.wbTrend.value = "linear";
        }
      }

      if (els.wbTemporal) {
        Array.from(els.wbTemporal.options).forEach((opt) => {
          const isBatchBy = opt.value.startsWith("by_");
          opt.hidden = mode === "prep" && isBatchBy;
        });
        if (mode === "prep" && els.wbTemporal.value.startsWith("by_")) {
          els.wbTemporal.value = "all";
        }
      }
    }

    function escHtml(val) {
      return String(val == null ? "" : val)
        .replace(/&/g, "&amp;")
        .replace(/</g, "&lt;")
        .replace(/>/g, "&gt;")
        .replace(/"/g, "&quot;");
    }

    function updateShellButtons() {
      document.querySelectorAll("#wbShellBtns [data-shell]").forEach((btn) => {
        btn.classList.toggle("active", btn.dataset.shell === shellMode);
      });
    }

    function formatArg(arg) {
      if (shellMode === "powershell") {
        if (/[,;\s"'$&|<>(){}@`~]/.test(arg)) {
          return `"${arg.replace(/"/g, "`\"")}"`;
        }
        return arg;
      }
      if (shellMode === "cmd") {
        if (/[\s"&|<>^()]/.test(arg)) {
          return `"${arg.replace(/"/g, '\\"')}"`;
        }
        return arg;
      }
      if (/[\s"'$&|;<>()!*?{}~]/.test(arg)) {
        return `'${arg.replace(/'/g, "'\\''")}'`;
      }
      return arg;
    }

    function formatCommand(runner, scriptPath, args) {
      const baseCmd = `${runner} ${scriptPath}`;
      if (!args.length) return baseCmd;
      const formattedArgs = args.map(formatArg);
      if (wrapFlags && formattedArgs.length >= 2) {
        const cont = shellMode === "powershell" ? " `" : shellMode === "cmd" ? " ^" : " \\";
        return baseCmd + cont + "\n  " + formattedArgs.join(cont + "\n  ");
      }
      return baseCmd + " " + formattedArgs.join(" ");
    }

    function getTargetDatasetsCount() {
      const counts = {};
      (window.WORKBENCH_TARGETS || []).forEach((grp) => {
        if (!grp.dataset || grp.dataset === "multi") return;
        const seenInGroup = new Set();
        (grp.items || []).forEach((it) => {
          if (!seenInGroup.has(it.value)) {
            seenInGroup.add(it.value);
            if (!counts[it.value]) counts[it.value] = new Set();
            counts[it.value].add(grp.dataset);
          }
        });
      });
      return counts;
    }

    const targetDatasetsMap = getTargetDatasetsCount();

    function populateTargetSelect(datasetFilter, mode) {
      const sel = els.wbTargetSel;
      const groups = window.WORKBENCH_TARGETS || [];
      if (!sel || !groups.length) return;

      const prevVal = sel.value || "0117";
      const prevOpt = sel.selectedOptions && sel.selectedOptions[0];
      const prevDs = prevOpt ? prevOpt.dataset.ds || "" : "";

      const visibleGroups = groups.filter((g) => {
        if (g.dataset === "multi") return datasetFilter === "auto" || datasetFilter === "wq";
        return datasetFilter === "auto" || g.dataset === datasetFilter;
      });

      sel.innerHTML = "";
      let matchedOpt = null;
      let firstOpt = null;

      visibleGroups.forEach((g) => {
        const og = document.createElement("optgroup");
        og.label = g.label;
        (g.items || []).forEach((it) => {
          const opt = document.createElement("option");
          opt.value = it.value;
          opt.textContent = it.label;
          opt.dataset.ds = it.ds || g.dataset || "";
          if (it.rsnOnly) opt.dataset.rsnOnly = "true";
          og.appendChild(opt);
          if (!firstOpt) firstOpt = opt;
          if (it.value === prevVal && (!prevDs || opt.dataset.ds === prevDs || datasetFilter !== "auto")) {
            if (!matchedOpt) matchedOpt = opt;
          }
        });
        if (og.children.length) sel.appendChild(og);
      });

      const customOpt = document.createElement("option");
      customOpt.value = "custom";
      customOpt.textContent = "Custom determinand code / CAS / index...";
      customOpt.dataset.ds = "";
      sel.appendChild(customOpt);

      if (prevVal === "custom") {
        sel.value = "custom";
      } else if (matchedOpt) {
        matchedOpt.selected = true;
      } else if (firstOpt) {
        firstOpt.selected = true;
      }
    }

    function applyState(state) {
      if (state.wbDataset && els.wbDataset) els.wbDataset.value = state.wbDataset;
      if (state.wbMode && els.wbMode) els.wbMode.value = state.wbMode;
      populateTargetSelect(els.wbDataset ? els.wbDataset.value : "auto", els.wbMode ? els.wbMode.value : "det");
      Object.entries(state).forEach(([k, v]) => {
        if (els[k]) els[k].value = v;
      });
      rebuild();
    }

    function rebuild() {
      const mode = els.wbMode ? els.wbMode.value : "det";
      syncFieldApplicability(mode);
      syncCustomVisibility();
      syncWindowsROptions(mode);

      const target = resolveSelectOrCustom(els.wbTargetSel, els.wbTargetCustom, "0117");
      const dataset = els.wbDataset ? els.wbDataset.value : "auto";
      const selectedOpt = els.wbTargetSel && els.wbTargetSel.selectedOptions ? els.wbTargetSel.selectedOptions[0] : null;
      const selectedTargetDs = selectedOpt ? (selectedOpt.dataset.ds || "") : "";
      const isRsnOnlyTarget = Boolean(selectedOpt && selectedOpt.dataset.rsnOnly === "true");
      const isSharedEcolCol = Boolean(targetDatasetsMap[target] && targetDatasetsMap[target].size > 1);
      const looksLikeLcmsCas = Boolean(
        els.wbTargetSel && els.wbTargetSel.value === "custom" &&
        (/^[0-9]{5,}$/.test(target) || /^[0-9]+-[0-9]+-[0-9]$/.test(target))
      );
      const effectiveDetDataset = dataset !== "auto"
        ? dataset
        : looksLikeLcmsCas
          ? "lcms"
          : (selectedTargetDs === "lcms" || (selectedTargetDs && !["wq", "multi"].includes(selectedTargetDs) && isSharedEcolCol))
            ? selectedTargetDs
            : "auto";
      const source = resolveSelectOrCustom(els.wbSourceSel, els.wbSourceCustom, "rsn");
      const network = resolveSelectOrCustom(els.wbNetworkSel, els.wbNetworkCustom, "orn");
      const transform = els.wbTransform ? els.wbTransform.value : "default";
      const norm = els.wbNorm ? els.wbNorm.value : "none";
      const censor = els.wbCensor ? els.wbCensor.value : "default";
      const trend = els.wbTrend ? els.wbTrend.value : "linear";
      const kernel = els.wbKernel ? els.wbKernel.value : "exponential";
      const afv = els.wbAfv ? els.wbAfv.value : "length";
      const cvMode = els.wbCvMode ? els.wbCvMode.value : "random";
      const thresholdMode = els.wbThresholdSel ? els.wbThresholdSel.value : "none";
      const thresholdVal = els.wbThresholdVal ? els.wbThresholdVal.value.trim() : "10";
      const thresholdDir = els.wbThresholdDir ? els.wbThresholdDir.value : "above";
      const temporal = els.wbTemporal ? els.wbTemporal.value : "all";
      const designSites = els.wbDesignSites ? els.wbDesignSites.value : "0";
      const basin = resolveSelectOrCustom(els.wbBasinSel, els.wbBasinCustom, "");
      const prefer = resolveSelectOrCustom(els.wbPreferSel, els.wbPreferCustom, "");
      const weight = resolveSelectOrCustom(els.wbWeightSel, els.wbWeightCustom, "");
      const sourceSe = resolveSelectOrCustom(els.wbSeSel, els.wbSeCustom, "");
      const fmtVal = els.wbFormat ? els.wbFormat.value : "png_200";
      const palette = els.wbPalette ? els.wbPalette.value : "standard";
      const omp = els.wbOmp ? els.wbOmp.value : "1";

      const normalizeSourceTag = (srcVal) => {
        const s = (srcVal || "rsn").trim().toLowerCase();
        const toks = s.split(/[,+]+/).map((t) => t.trim()).filter(Boolean).sort();
        if (toks.length === 2 && toks[0] === "rsn" && toks[1] === "ssn") return "rss";
        return s;
      };

      const formatVarsFolderTag = (specStr) => {
        const specs = (specStr || "").split(",").map((s) => s.trim()).filter(Boolean);
        if (!specs.length) return "";
        const slugMap = {
          elev: "Elev", slope: "Slope", urban: "Urban", chalk: "Chalk",
          sand: "Sand", clay: "Clay", aquifer: "Aquifer", peat: "Peat",
          easting: "Easting", northing: "Northing", sinuosity: "Sinuosity",
          srcelevation: "SrcElevation", dist2mth: "Dist2Mth",
          length: "Length", shreve: "Shreve", strahler: "Strahler"
        };
        const toSlug = (tok) => {
          const clean = tok.toLowerCase().replace(/[^a-z0-9_]+/g, "");
          return slugMap[clean] || tok.replace(/[^a-zA-Z0-9_]+/g, "");
        };
        if (specs.every((sp) => /^[+-]/.test(sp))) {
          return `_vars_${specs.map((sp) => `${sp[0] === "+" ? "add" : "no"}_${toSlug(sp.slice(1))}`).join("_")}`;
        }
        const terms = specs.map((sp) => toSlug(sp.replace(/^\+/, ""))).filter(Boolean);
        if (!terms.length) return "";
        const capped = terms.length > 5 ? [...terms.slice(0, 4), `plus${terms.length - 4}`] : terms;
        return `_vars_${capped.join("_")}`;
      };

      const args = [];
      const notes = [];
      let runner = "Rscript";
      let scriptRel = "scripts/run_s3n_determinand.R";

      if (mode === "prep") {
        runner = shellMode === "bash" ? "python3" : "python";
        scriptRel = "scripts/01_preprocess_orn.py";
        const prepTarget = (target.includes(",") ? target.split(",")[0].trim() : target) || "0117";
        const prepDataset = effectiveDetDataset !== "auto"
          ? effectiveDetDataset
          : (/^[0-9]{4}$/.test(prepTarget) ? "wq" : (selectedTargetDs || "wq"));
        if (source && source !== "rsn") args.push(`--source=${source}`);
        if (prepDataset && prepDataset !== "wq") args.push(`--dataset=${prepDataset}`);
        if (prepTarget !== "0117" || prepDataset !== "wq") args.push(`--det-code=${prepTarget}`);
        if (network && network !== "orn") args.push(`--network=${network}`);
        if (transform && transform !== "default") {
          args.push(`--transform=${transform}`);
        } else if (["inverts", "diatoms", "macrophytes", "rhs"].includes(prepDataset)) {
          args.push("--transform=none");
        }
        if (afv && afv !== "length") args.push(`--afv-scheme=${afv}`);
        if (censor && censor !== "default") {
          if (censor === "em") {
            notes.push("Note: Spatial Tobit EM kriging runs at model-fit time in <code>run_s3n_determinand.R</code>; Python preprocessing defaults to <code>ros</code>.");
          } else {
            args.push(`--censor-method=${censor}`);
          }
        }
        if (temporal === "custom_season" && els.wbTemporalSeason) args.push(`--season=${els.wbTemporalSeason.value}`);
        else if (temporal === "custom_year" && els.wbTemporalYear && els.wbTemporalYear.value.trim()) args.push(`--year=${els.wbTemporalYear.value.trim()}`);
        else if (temporal === "custom_month" && els.wbTemporalMonth) args.push(`--month=${els.wbTemporalMonth.value}`);
        if (source.toLowerCase().includes("waterblitz")) {
          args.push("--max-snap-dist=350");
          notes.push("Citizen-science WaterBlitz sources automatically benefit from <code>--max-snap-dist=350</code> snapping tolerance.");
        }
        notes.push("Writes preprocessed reach topology, AFVs, covariates, and snapped station tables under <code>data/processed/</code> without R/GDAL spatial dependencies.");
      } else if (mode === "stats") {
        runner = "Rscript";
        scriptRel = "scripts/05_descriptive_stats.R";
        const dir1 = (els.wbPostDir1 ? els.wbPostDir1.value.trim() : "") || (postDir1Sel ? postDir1Sel.value : "") || "output/results/det_0117_orn_rsn";
        const postLabel = els.wbPostLabels ? els.wbPostLabels.value.trim() : "";
        const postOut = els.wbPostOut ? els.wbPostOut.value.trim() : "";
        args.push(dir1);
        if (postLabel) args.push(`--label=${postLabel}`);
        if (postOut) args.push(`--out=${postOut}`);
        if (fmtVal === "svg") args.push("--fig-format=svg");
        else if (fmtVal === "pdf") args.push("--fig-format=pdf");
        else if (fmtVal === "png_300") args.push("--dpi=300");
        notes.push("Computes detailed descriptive statistics, representativeness tests, and stratified environmental/topological tables for a single saved model run.");
        notes.push(`Writes outputs to <code>${escHtml(postOut || `${dir1}/descriptives/`)}</code>.`);
      } else if (mode === "compare") {
        runner = "Rscript";
        scriptRel = "scripts/04_compare_model_runs.R";
        const dir1 = (els.wbPostDir1 ? els.wbPostDir1.value.trim() : "") || (postDir1Sel ? postDir1Sel.value : "") || "output/results/det_0117_orn_rsn";
        const cmpMode = els.wbCompareMode ? els.wbCompareMode.value : "dirs";
        const dir2Raw = (els.wbPostDir2 ? els.wbPostDir2.value.trim() : "") || (postDir2Sel ? postDir2Sel.value : "");
        const postLabels = els.wbPostLabels ? els.wbPostLabels.value.trim() : "";
        const postOut = els.wbPostOut ? els.wbPostOut.value.trim() : "";
        args.push(dir1);
        if (cmpMode === "dirs") {
          const extraDirs = (dir2Raw || "output/results/det_0180_orn_rsn")
            .split(/[\s,]+/)
            .map((s) => s.trim())
            .filter(Boolean);
          extraDirs.forEach((d) => args.push(d));
        } else {
          args.push(`--by=${cmpMode}`);
        }
        if (postLabels) args.push(`--labels=${postLabels}`);
        if (postOut) args.push(`--out=${postOut}`);
        if (fmtVal === "svg") args.push("--fig-format=svg");
        else if (fmtVal === "pdf") args.push("--fig-format=pdf");
        else if (fmtVal === "png_300") args.push("--dpi=300");
        notes.push("Compares CV performance, parameter CIs, reach-level predictions, and stratified diagnostics across runs or temporal slices, rendering <code>comparison_report.md</code> and <code>fig_comp1</code>–<code>fig_comp5</code>.");
        notes.push(`Writes outputs to <code>${escHtml(postOut || "output/comparisons/")}</code>.`);
      } else if (mode === "snap") {
        runner = "Rscript";
        scriptRel = "scripts/06_station_snapping_diagnostics.R";
        const snapCatchment = resolveSelectOrCustom(els.wbSnapCatchmentSel, els.wbSnapCatchmentCustom, "all");
        const snapOut = els.wbSnapOut ? els.wbSnapOut.value.trim() : "";
        if (source && source !== "rsn") args.push(`--source=${source}`);
        if (network && network !== "orn") args.push(`--network=${network}`);
        if (snapCatchment && snapCatchment !== "all") args.push(`--catchment=${snapCatchment}`);
        if (snapOut) args.push(`--out=${snapOut}`);
        if (fmtVal === "svg") args.push("--fig-format=svg");
        else if (fmtVal === "pdf") args.push("--fig-format=pdf");
        else if (fmtVal === "png_300") args.push("--dpi=300");
        notes.push("Audits station-to-reach snapping distances, generates <code>station_snapping_all.csv</code>, <code>catchment_snapping_summary.csv</code>, <code>station_snapping_report.md</code>, and metre-scale catchment zoom maps.");
        notes.push(`Writes outputs to <code>${escHtml(snapOut || `output/station_snapping/${network || "orn"}/`)}</code>.`);
      } else if (mode === "plots") {
        runner = "Rscript";
        scriptRel = "scripts/03_generate_plots.R";
        const plotDir = (els.wbPlotDir ? els.wbPlotDir.value.trim() : "") || (plotDirSel ? plotDirSel.value : "") || "output/results/det_0117_orn_rsn";
        const plotBasin = els.wbPlotBasin ? els.wbPlotBasin.value.trim() : "";
        const plotVmin = els.wbPlotVmin ? els.wbPlotVmin.value.trim() : "";
        const plotVmax = els.wbPlotVmax ? els.wbPlotVmax.value.trim() : "";
        const plotFigSize = els.wbPlotFigSize ? els.wbPlotFigSize.value.trim() : "";
        const plotNorm = els.wbPlotNorm ? els.wbPlotNorm.value : "none";
        const plotPalette = els.wbPlotPalette ? els.wbPlotPalette.value : "auto";
        const plotFmt = els.wbPlotFormat ? els.wbPlotFormat.value : "png_200";
        args.push(`--dir=${plotDir}`);
        if (plotBasin) args.push(`--basin=${plotBasin}`);
        if (plotVmin !== "") args.push(`--vmin=${plotVmin}`);
        if (plotVmax !== "") args.push(`--vmax=${plotVmax}`);
        if (plotFigSize) args.push(`--fig-size=${plotFigSize}`);
        if (plotNorm && plotNorm !== "none") args.push(`--norm=${plotNorm}`);
        if (plotPalette === "standard") args.push("--normal-color");
        else if (plotPalette === "reverse") args.push("--reverse-color");
        if (plotFmt === "svg") args.push("--fig-format=svg");
        else if (plotFmt === "pdf") args.push("--fig-format=pdf");
        else if (plotFmt === "png_300") args.push("--dpi=300");
        notes.push("Re-renders Figures 1–4 (plus custom basin maps and synchronized temporal color bars) from saved prediction CSVs without re-fitting models.");
        notes.push(`Reads and writes figures inside <code>${escHtml(plotDir)}</code>.`);
      } else if (mode === "cokrig") {
        runner = "Rscript";
        scriptRel = "scripts/run_s3n_cokriging.R";
        const singleCoSrcDs = effectiveDetDataset !== "auto"
          ? effectiveDetDataset
          : (/^[0-9]{4}$/.test(target) ? "wq" : (selectedTargetDs || "wq"));
        const coTargets = target.includes(",")
          ? target
          : (target === "0180"
            ? "0117,0180"
            : (singleCoSrcDs !== "wq" ? `${singleCoSrcDs}:${target},wq:0180` : `${target},0180`));
        args.push(`--targets=${coTargets}`);
        if (dataset && dataset !== "auto" && dataset !== "wq" && !coTargets.includes(":")) args.push(`--dataset=${dataset}`);
        if (source && source !== "rsn") args.push(`--source=${source}`);
        if (network && network !== "orn") args.push(`--network=${network}`);
        if (trend === "colocated") args.push("--method=colocated");
        else args.push("--method=lmc");
        if (kernel && kernel !== "exponential") args.push(`--kernel=${kernel}`);
        if (cvMode && cvMode !== "random") args.push(`--cv-mode=${cvMode}`);
        if (thresholdMode === "uk") {
          args.push("--threshold=uk");
        } else if (thresholdMode === "custom" && thresholdVal !== "") {
          args.push(`--threshold=${thresholdVal}`);
          if (thresholdDir === "below") args.push("--threshold-dir=below");
        }
        if (prefer) args.push(`--prefer-source=${prefer}`);
        if (weight) args.push(`--source-weight=${weight}`);
        if (sourceSe) args.push(`--source-se=${sourceSe}`);
        if (fmtVal === "svg") args.push("--fig-format=svg");
        else if (fmtVal === "pdf") args.push("--fig-format=pdf");
        else if (fmtVal === "png_300") args.push("--dpi=300");
        if (omp && omp !== "1") args.push(`--omp=${omp}`);

        const coTokenParts = coTargets.split(",").map((s) => s.trim()).filter(Boolean);
        const coSlugs = coTokenParts.map((s) => (s.includes(":") ? s.split(":").slice(1).join(":") : s).replace(/[^a-zA-Z0-9_]+/g, "_")).join("_");
        const coTag = coSlugs.length > 64 ? `${coSlugs.slice(0, 52)}_p${coTokenParts.length}` : coSlugs;
        const coMethodTag = trend === "colocated" ? "_colocated" : (coTokenParts.length > 3 ? "_k3" : "");
        const coSourceTag = normalizeSourceTag(source);
        const coOutDir = `output/results/cokrig_${coTag}_${network || "orn"}_${coSourceTag}${coMethodTag}/`;
        notes.push("Fits Multivariate Factor-LMC or Co-Located Co-Kriging with Dual CV (cross-target univariate vs. joint gain).");
        notes.push(`Expected output folder: <code>${escHtml(coOutDir)}</code>.`);
        if (transform !== "default") {
          notes.push("Note: Co-Kriging applies each target's configured transform in <code>configs/determinands.json</code> automatically.");
        }
      } else {
        runner = "Rscript";
        scriptRel = "scripts/run_s3n_determinand.R";
        const singleTarget = (target.includes(",") ? target.split(",")[0].trim() : target) || "0117";
        args.push(singleTarget);
        if (effectiveDetDataset && effectiveDetDataset !== "auto") args.push(`--dataset=${effectiveDetDataset}`);
        if (source && source !== "rsn") args.push(`--source=${source}`);
        if (network && network !== "orn") args.push(`--network=${network}`);
        if (transform && transform !== "default") args.push(`--transform=${transform}`);
        if (norm && norm !== "none") args.push(`--norm=${norm}`);
        if (censor && censor !== "default") args.push(`--censor=${censor}`);
        if (trend && trend !== "linear") args.push(`--trend=${trend}`);
        if (kernel && kernel !== "exponential") args.push(`--kernel=${kernel}`);
        if (afv && afv !== "length") args.push(`--afv=${afv}`);
        const combineMode = els.wbCombine ? els.wbCombine.value : "joint";
        if (combineMode && combineMode !== "joint") args.push(`--combine=${combineMode}`);
        const varsSpec = resolveSelectOrCustom(els.wbVarsSel, els.wbVarsCustom, "");
        if (varsSpec) args.push(`--vars=${varsSpec}`);
        if (cvMode && cvMode !== "random") args.push(`--cv-mode=${cvMode}`);
        const bootVal = els.wbBoot ? els.wbBoot.value : "50";
        if (bootVal && bootVal !== "50") args.push(`--boot=${bootVal}`);

        if (thresholdMode === "uk") {
          args.push("--threshold=uk");
        } else if (thresholdMode === "custom" && thresholdVal !== "") {
          args.push(`--threshold=${thresholdVal}`);
          if (thresholdDir === "below") args.push("--threshold-dir=below");
        }

        if (temporal === "by_season") args.push("--by=season", "--st-borrow");
        else if (temporal === "by_year") args.push("--by=year", "--st-borrow");
        else if (temporal === "by_month") args.push("--by=month", "--st-borrow");
        else if (temporal === "custom_season" && els.wbTemporalSeason) args.push(`--season=${els.wbTemporalSeason.value}`);
        else if (temporal === "custom_year" && els.wbTemporalYear && els.wbTemporalYear.value.trim()) args.push(`--year=${els.wbTemporalYear.value.trim()}`);
        else if (temporal === "custom_month" && els.wbTemporalMonth) args.push(`--month=${els.wbTemporalMonth.value}`);

        if (designSites !== "0") args.push(`--design-sites=${designSites}`);
        if (basin && basin !== "none") args.push(`--basin=${basin}`);
        if (prefer) args.push(`--prefer-source=${prefer}`);
        if (weight) args.push(`--source-weight=${weight}`);
        if (sourceSe) args.push(`--source-se=${sourceSe}`);

        if (fmtVal === "svg") args.push("--fig-format=svg");
        else if (fmtVal === "pdf") args.push("--fig-format=pdf");
        else if (fmtVal === "png_300") args.push("--dpi=300");

        if (palette === "reverse") args.push("--reverse-color");

        if (omp && omp !== "1") args.push(`--omp=${omp}`);

        const inferredDs = effectiveDetDataset !== "auto" ? effectiveDetDataset : (/^[0-9]{4}$/.test(singleTarget) ? "wq" : (selectedTargetDs || "ecol"));
        const folderPrefix = inferredDs === "wq" ? `det_${singleTarget}` : `${inferredDs}_${singleTarget.replace(/[^a-zA-Z0-9_]+/g, "_")}`;
        const defaultTrans = (inferredDs === "wq" || inferredDs === "lcms") ? "log" : "none";
        let expectedFolder = `${folderPrefix}_${network || "orn"}_${normalizeSourceTag(source)}`;
        if (transform && transform !== "default" && transform !== defaultTrans) {
          expectedFolder += `_trans_${transform === "none" ? "raw" : transform}`;
        }
        if (varsSpec) expectedFolder += formatVarsFolderTag(varsSpec);
        if (kernel && kernel !== "exponential") expectedFolder += `_k${kernel}`;
        if (combineMode && combineMode !== "joint") expectedFolder += `_c${combineMode}`;
        if (afv && afv !== "length") expectedFolder += `_afv_${afv}`;
        if (trend && trend !== "linear") expectedFolder += `_trend_${trend}`;
        if (cvMode && cvMode !== "random") expectedFolder += `_cv_${cvMode}`;
        if (censor === "em") expectedFolder += "_tobit";
        if (["by_season", "by_year", "by_month"].includes(temporal)) expectedFolder += "_st";
        let subSlice = "";
        if (temporal === "custom_season" && els.wbTemporalSeason) subSlice = `season_${els.wbTemporalSeason.value}/`;
        else if (temporal === "custom_year" && els.wbTemporalYear && els.wbTemporalYear.value.trim()) subSlice = `year_${els.wbTemporalYear.value.trim()}/`;
        else if (temporal === "custom_month" && els.wbTemporalMonth) subSlice = `month_${String(els.wbTemporalMonth.value).padStart(2, "0")}/`;

        notes.push("Outputs full network reach predictions (<code>predictions_*.csv</code>), CV metrics (<code>cv_comparison_*.csv</code>), parameter estimates (<code>parameters_*.csv</code>), and publication maps.");
        notes.push(`Expected output folder: <code>output/results/${escHtml(expectedFolder)}/${escHtml(subSlice)}</code>`);
      }

      function resolveRunnerInfo(baseRunner) {
        if (shellMode !== "bash" && baseRunner === "Rscript" && winFullPathEl && winFullPathEl.checked) {
          const rawVer = resolveSelectOrCustom(els.wbRVersionSel, els.wbRVersionCustom, "4.4.3");
          const looksLikePath = /[\\/]|\.exe$/i.test(rawVer);
          const rExePath = looksLikePath
            ? rawVer
            : `C:\\Program Files\\R\\R-${rawVer.replace(/^R-/i, "")}\\bin\\x64\\Rscript.exe`;
          return {
            runner: shellMode === "powershell" ? `& "${rExePath}"` : `"${rExePath}"`,
            rExePath
          };
        }
        return { runner: baseRunner, rExePath: "" };
      }

      const scriptPath = shellMode === "bash" ? scriptRel : scriptRel.replace(/\//g, "\\");
      const shellLabel = shellMode === "powershell" ? "PowerShell" : shellMode === "cmd" ? "Command Prompt" : unixShellLabel;
      if (targetScriptEl) targetScriptEl.textContent = `${shellLabel} · ${scriptPath}`;

      const resolvedRunner = resolveRunnerInfo(runner);
      runner = resolvedRunner.runner;
      if (resolvedRunner.rExePath) {
        notes.push(
          shellMode === "powershell"
            ? `Uses PowerShell call operator <code>&amp;</code> around <code>${escHtml(resolvedRunner.rExePath)}</code> so <code>Rscript</code> runs without being on Windows <code>PATH</code>.`
            : `Quotes <code>${escHtml(resolvedRunner.rExePath)}</code> for Windows Command Prompt so <code>Rscript</code> runs without being on Windows <code>PATH</code>.`
        );
      }

      if (mode === "det" || mode === "cokrig") {
        if (thresholdMode !== "none") {
          notes.push("Generates 5-tier Gaussian exceedance risk maps and catchment class area breakdowns (<strong>Uncertain</strong> band highlights monitoring priorities).");
        }
        if (designSites !== "0" && mode === "det") {
          notes.push(`Runs OpenMP C++ greedy Schur-complement selection for the top <strong>K = ${designSites}</strong> reaches and renders <code>fig_optimal_design_*</code> maps.`);
        }
        if (trend === "gam" && mode === "det") {
          notes.push("Fits REML penalised splines plus orthogonal out-of-fold spatial residual kriging and generates <code>fig_gam_effects_*</code>.");
        }
        if (censor === "em") {
          notes.push("Activates iterative <strong>Spatial Tobit EM kriging</strong> for left-censored <code>&lt; LOD</code> and interval/range citizen data.");
        }
        if (prefer || weight || sourceSe) {
          notes.push("Applies hierarchical multi-source de-duplication and analytic variance inflation via <code>s3n_fit_ml_weighted</code>.");
        }
        if (source === "ssn" && (isRsnOnlyTarget || effectiveDetDataset === "rhs" || effectiveDetDataset === "lcms")) {
          notes.push("Warning: The selected target/modality is monitored on <code>rsn</code> (River Surveillance Network) only and is unavailable in <code>--source=ssn</code> alone.");
        }
      }

      function updateListButtons() {
        document.querySelectorAll("#wbListBtns [data-list-cmd]").forEach((b) => {
          b.classList.toggle("active", b.dataset.listCmd === activeListCmd);
        });
      }

      if (activeListCmd) {
        const listSpec = buildListSpec(activeListCmd, source, network, effectiveDetDataset);
        const listSp = shellMode === "bash" ? listSpec.script : listSpec.script.replace(/\//g, "\\");
        const listRunnerInfo = resolveRunnerInfo("Rscript");
        if (targetScriptEl) targetScriptEl.textContent = `${shellLabel} · ${listSp}`;
        const listNotes = [
          listSpec.desc,
          "Read-only discovery lookup (writes no model files). Click the highlighted button again or edit any Configure Analysis control to return to the pipeline command."
        ];
        if (listRunnerInfo.rExePath) {
          listNotes.push(
            shellMode === "powershell"
              ? `Uses PowerShell call operator <code>&amp;</code> around <code>${escHtml(listRunnerInfo.rExePath)}</code> so <code>Rscript</code> runs without being on Windows <code>PATH</code>.`
              : `Quotes <code>${escHtml(listRunnerInfo.rExePath)}</code> for Windows Command Prompt so <code>Rscript</code> runs without being on Windows <code>PATH</code>.`
          );
        }
        outEl.textContent = formatCommand(listRunnerInfo.runner, listSp, listSpec.args);
        if (notesEl) {
          notesEl.innerHTML = listNotes.map((n) => `<li>${n}</li>`).join("");
        }
        updateListButtons();
        return;
      }

      outEl.textContent = formatCommand(runner, scriptPath, args);
      if (notesEl) {
        notesEl.innerHTML = notes.map((n) => `<li>${n}</li>`).join("");
      }
      updateListButtons();
      if (typeof syncHashState === "function") syncHashState();
    }

    function buildListSpec(cmdKey, source, network, effectiveDetDataset) {
      const specs = {
        sources: {
          script: "scripts/run_s3n_determinand.R",
          args: ["--list-sources"],
          desc: "Prints all built-in (<code>rsn</code>, <code>ssn</code>, <code>rss</code>) and registered custom monitoring sources (<code>configs/sources/*.json</code>), abbreviations, and data directories."
        },
        networks: {
          script: "scripts/run_s3n_determinand.R",
          args: ["--list-networks"],
          desc: "Lists built-in and custom BYOWN stream network profiles (<code>configs/networks/*.json</code>). Also available via <code>01_preprocess_orn.py --list-networks</code>."
        },
        determinands: {
          script: "scripts/run_s3n_determinand.R",
          args: [
            "--list",
            ...(effectiveDetDataset && effectiveDetDataset !== "auto" && effectiveDetDataset !== "wq" ? [`--dataset=${effectiveDetDataset}`] : []),
            ...(source && source !== "rsn" ? [`--source=${source}`] : [])
          ],
          desc: "Lists available targets, sample counts, and units for the active dataset (<code>wq</code>, <code>inverts</code>, <code>diatoms</code>, <code>macrophytes</code>, <code>rhs</code>, or <code>lcms</code>) and source."
        },
        basins: {
          script: "scripts/run_s3n_determinand.R",
          args: [
            "--list-basins",
            ...(source && source !== "rsn" ? [`--source=${source}`] : []),
            ...(network && network !== "orn" ? [`--network=${network}`] : [])
          ],
          desc: "Lists all 11 Environment Agency River Basin Districts (<code>RBD_ID</code>, <code>RBD_Name</code>) and Top-30 monitored catchments (<code>CatchID</code>, <code>Main_River</code>) for <code>--basin=</code>."
        },
        catchments: {
          script: "scripts/06_station_snapping_diagnostics.R",
          args: [
            "--list-catchments",
            ...(source && source !== "rsn" ? [`--source=${source}`] : []),
            ...(network && network !== "orn" ? [`--network=${network}`] : [])
          ],
          desc: "Lists all river catchments and station-snapping counts available for snapping diagnostics (<code>--catchment=</code> / <code>--basin=</code>)."
        },
        outputs: {
          script: "scripts/05_descriptive_stats.R",
          args: ["--list"],
          desc: "Scans <code>output/results/</code> and lists fitted model output directories available for single-run descriptive statistics (<code>05_descriptive_stats.R</code>)."
        },
        comparisons: {
          script: "scripts/04_compare_model_runs.R",
          args: ["--list"],
          desc: "Scans <code>output/results/</code> and lists model output directories available for multi-run or N:P stoichiometry comparisons (<code>04_compare_model_runs.R</code>)."
        },
        vars: {
          script: "scripts/run_s3n_determinand.R",
          args: ["--list-vars"],
          desc: "Lists all 15 catchment &amp; topological covariates (<code>COVARIATE_CATALOG</code>), short CLI aliases, and default inclusion status (<code>--vars=+...,-...</code>)."
        },
        transforms: {
          script: "scripts/run_s3n_determinand.R",
          args: ["--list-transforms"],
          desc: "Prints all 12 supported response transformations (<code>log</code>, <code>boxcox</code>, <code>yeojohnson</code>, <code>sqrt</code>, <code>logit</code>, etc.), formulas, and valid domains (<code>--transform=</code>)."
        },
        kernels: {
          script: "scripts/run_s3n_determinand.R",
          args: ["--list-kernels"],
          desc: "Lists the 5 spatial stream covariance kernels (<code>exponential</code>, <code>spherical</code>, <code>mariah</code>, <code>linear</code>, <code>auto</code> Profile MLE) accepted by <code>--kernel=</code>."
        },
        trends: {
          script: "scripts/run_s3n_determinand.R",
          args: ["--list-trends"],
          desc: "Lists non-spatial mean/trend engines (<code>linear</code>, <code>gam</code>, <code>spline</code>, <code>rf</code>, <code>gbm</code>) accepted by <code>--trend=</code> before orthogonal S3N residual kriging."
        },
        cv: {
          script: "scripts/run_s3n_determinand.R",
          args: ["--list-cv"],
          desc: "Lists the 5 hydrological cross-validation holdout strategies (<code>random</code>, <code>spatial</code>, <code>catchment</code>, <code>downstream</code>, <code>upstream</code>) for <code>--cv-mode=</code>."
        },
        thresholds: {
          script: "scripts/run_s3n_determinand.R",
          args: ["--list-thresholds"],
          desc: "Displays built-in indicative UK threshold lookup values (<code>--threshold=uk</code>) and direction defaults (<code>--threshold-dir=above|below</code>)."
        }
      };

      return specs[cmdKey] || specs.sources;
    }

    document.querySelectorAll("#wbListBtns [data-list-cmd]").forEach((btn) => {
      btn.addEventListener("click", () => {
        const clicked = btn.dataset.listCmd || "";
        activeListCmd = activeListCmd === clicked ? "" : clicked;
        rebuild();
      });
    });

    const postFieldIds = new Set(["wbPostDir1", "wbCompareMode", "wbPostDir2", "wbPostLabels", "wbPostOut"]);
    const snapFieldIds = new Set(["wbSnapCatchmentSel", "wbSnapCatchmentCustom", "wbSnapOut"]);

    document.querySelectorAll("#wbPostActionBtns [data-post-mode]").forEach((btn) => {
      btn.addEventListener("click", () => {
        activeListCmd = "";
        if (els.wbMode) els.wbMode.value = btn.dataset.postMode || "stats";
        rebuild();
      });
    });

    if (activateSnapBtn) {
      activateSnapBtn.addEventListener("click", () => {
        activeListCmd = "";
        if (els.wbMode) els.wbMode.value = "snap";
        rebuild();
      });
    }

    const plotFieldIds = new Set(["wbPlotDir", "wbPlotBasin", "wbPlotVmin", "wbPlotVmax", "wbPlotFigSize", "wbPlotNorm", "wbPlotPalette", "wbPlotFormat"]);
    if (activatePlotBtn) {
      activatePlotBtn.addEventListener("click", () => {
        activeListCmd = "";
        if (els.wbMode) els.wbMode.value = "plots";
        rebuild();
      });
    }

    ids.forEach((id) => {
      if (!els[id]) return;
      const evt = els[id].tagName === "INPUT" ? "input" : "change";
      els[id].addEventListener(evt, () => {
        activeListCmd = "";
        if (postFieldIds.has(id) && els.wbMode) {
          if (id === "wbCompareMode" || id === "wbPostDir2") {
            els.wbMode.value = "compare";
          } else if (els.wbMode.value !== "stats" && els.wbMode.value !== "compare") {
            els.wbMode.value = "stats";
          }
        } else if (snapFieldIds.has(id) && els.wbMode) {
          els.wbMode.value = "snap";
        } else if ((id === "wbTargetSel" || id === "wbTargetCustom") && els.wbMode && els.wbMode.value === "det") {
          const selOpt = els.wbTargetSel && els.wbTargetSel.selectedOptions ? els.wbTargetSel.selectedOptions[0] : null;
          const curTarget = resolveSelectOrCustom(els.wbTargetSel, els.wbTargetCustom, "0117");
          if ((selOpt && selOpt.dataset.ds === "multi") || curTarget.includes(",")) {
            els.wbMode.value = "cokrig";
          }
        }
        if (id === "wbDataset" || id === "wbMode") {
          populateTargetSelect(els.wbDataset ? els.wbDataset.value : "auto", els.wbMode ? els.wbMode.value : "det");
        }
        rebuild();
      });
    });

    const WIN_R_KEY = "s3n_wb_win_r";
    const HISTORY_KEY = "s3n_wb_cmd_history";
    const postDir2SelRow = document.getElementById("wbPostDir2SelRow");
    const postDir2AddBtn = document.getElementById("wbPostDir2AddBtn");
    const historyBtn = document.getElementById("wbHistoryBtn");
    const historyMenu = document.getElementById("wbHistoryMenu");
    let scannedDirs = [];
    let isHydrating = false;

    function loadWinRPrefs() {
      try {
        const raw = localStorage.getItem(WIN_R_KEY);
        if (!raw) return;
        const prefs = JSON.parse(raw);
        if (winFullPathEl && typeof prefs.fullPath === "boolean") winFullPathEl.checked = prefs.fullPath;
        if (els.wbRVersionSel && prefs.version) els.wbRVersionSel.value = prefs.version;
        if (els.wbRVersionCustom && typeof prefs.custom === "string") els.wbRVersionCustom.value = prefs.custom;
      } catch (e) {}
    }

    function saveWinRPrefs() {
      try {
        localStorage.setItem(WIN_R_KEY, JSON.stringify({
          fullPath: Boolean(winFullPathEl && winFullPathEl.checked),
          version: els.wbRVersionSel ? els.wbRVersionSel.value : "4.4.3",
          custom: els.wbRVersionCustom ? els.wbRVersionCustom.value : ""
        }));
      } catch (e) {}
    }

    function syncHashState() {
      if (isHydrating || typeof history === "undefined" || !history.replaceState) return;
      try {
        const params = new URLSearchParams();
        if (shellMode !== (isWindows ? "powershell" : "bash")) params.set("shell", shellMode);
        if (wrapFlags) params.set("wrap", "1");
        if (activeSectionTab !== "model") params.set("tab", activeSectionTab);
        if (activeListCmd) params.set("list", activeListCmd);
        ids.forEach((id) => {
          if (!els[id] || id === "wbRVersionSel" || id === "wbRVersionCustom") return;
          if (els[id].value !== DEFAULTS[id]) params.set(id, els[id].value);
        });
        const str = params.toString();
        const nextHash = str ? `#${str}` : "";
        if (window.location.hash !== nextHash) {
          history.replaceState(null, "", nextHash || window.location.pathname + window.location.search);
        }
      } catch (e) {}
    }

    function hydrateFromHash() {
      try {
        const h = (window.location.hash || "").replace(/^#/, "");
        if (!h) return false;
        const params = new URLSearchParams(h);
        isHydrating = true;
        if (params.has("shell")) shellMode = params.get("shell") || shellMode;
        if (params.has("wrap")) {
          wrapFlags = params.get("wrap") === "1";
          if (wrapBtn) {
            wrapBtn.textContent = wrapFlags ? "Wrap Flags: On" : "Wrap Flags: Off";
            wrapBtn.classList.toggle("active", wrapFlags);
          }
        }
        if (params.has("wbDataset") && els.wbDataset) els.wbDataset.value = params.get("wbDataset");
        if (params.has("wbMode") && els.wbMode) els.wbMode.value = params.get("wbMode");
        populateTargetSelect(els.wbDataset ? els.wbDataset.value : "auto", els.wbMode ? els.wbMode.value : "det");
        ids.forEach((id) => {
          if (id === "wbRVersionSel" || id === "wbRVersionCustom") return;
          if (params.has(id) && els[id]) els[id].value = params.get(id);
        });
        if (params.has("list")) activeListCmd = params.get("list") || "";
        if (params.has("tab")) activeSectionTab = params.get("tab") || "model";
        isHydrating = false;
        return true;
      } catch (e) {
        isHydrating = false;
        return false;
      }
    }

    function loadHistory() {
      try {
        const arr = JSON.parse(localStorage.getItem(HISTORY_KEY) || "[]");
        return Array.isArray(arr) ? arr.slice(0, 8) : [];
      } catch (e) { return []; }
    }

    function renderHistory() {
      const items = loadHistory();
      if (!historyBtn || !historyMenu) return;
      historyBtn.textContent = `History (${items.length})`;
      historyBtn.hidden = items.length === 0;
      if (!items.length) {
        historyMenu.hidden = true;
        historyBtn.setAttribute("aria-expanded", "false");
        return;
      }
      historyMenu.innerHTML = "";
      items.forEach((cmd) => {
        const b = document.createElement("button");
        b.type = "button";
        b.className = "wb-history-item";
        const preview = cmd.replace(/\s+/g, " ").trim();
        b.textContent = preview.length > 90 ? preview.slice(0, 87) + "..." : preview;
        b.title = cmd;
        b.addEventListener("click", () => {
          const done = () => {
            historyMenu.hidden = true;
            historyBtn.setAttribute("aria-expanded", "false");
            const orig = historyBtn.textContent;
            historyBtn.textContent = "Copied!";
            setTimeout(() => { historyBtn.textContent = orig; }, 1200);
          };
          if (navigator.clipboard && navigator.clipboard.writeText) {
            navigator.clipboard.writeText(cmd).then(done).catch(() => fallbackCopy(cmd, done));
          } else {
            fallbackCopy(cmd, done);
          }
        });
        historyMenu.appendChild(b);
      });
      const clearBtn = document.createElement("button");
      clearBtn.type = "button";
      clearBtn.className = "wb-history-clear";
      clearBtn.textContent = "Clear History";
      clearBtn.addEventListener("click", () => {
        try { localStorage.removeItem(HISTORY_KEY); } catch (e) {}
        renderHistory();
      });
      historyMenu.appendChild(clearBtn);
    }

    window.__wbRecordHistory = function (cmdText) {
      const clean = (cmdText || "").trim();
      if (!clean) return;
      const prev = loadHistory().filter((c) => c !== clean);
      const next = [clean, ...prev].slice(0, 8);
      try { localStorage.setItem(HISTORY_KEY, JSON.stringify(next)); } catch (e) {}
      renderHistory();
    };

    if (historyBtn && historyMenu) {
      historyBtn.addEventListener("click", (e) => {
        e.stopPropagation();
        const willOpen = historyMenu.hidden;
        historyMenu.hidden = !willOpen;
        historyBtn.setAttribute("aria-expanded", willOpen ? "true" : "false");
      });
      document.addEventListener("click", (e) => {
        if (!e.target.closest("#wbHistoryWrap")) {
          historyMenu.hidden = true;
          historyBtn.setAttribute("aria-expanded", "false");
        }
      });
    }

    if (winFullPathEl) {
      winFullPathEl.addEventListener("change", () => {
        saveWinRPrefs();
        rebuild();
      });
    }
    [els.wbRVersionSel, els.wbRVersionCustom].forEach((el) => {
      if (el) el.addEventListener(el.tagName === "INPUT" ? "input" : "change", saveWinRPrefs);
    });

    document.querySelectorAll("#wbShellBtns [data-shell]").forEach((btn) => {
      btn.addEventListener("click", () => {
        shellMode = btn.dataset.shell || "bash";
        updateShellButtons();
        rebuild();
      });
    });

    if (wrapBtn) {
      wrapBtn.addEventListener("click", () => {
        wrapFlags = !wrapFlags;
        wrapBtn.textContent = wrapFlags ? "Wrap Flags: On" : "Wrap Flags: Off";
        wrapBtn.classList.toggle("active", wrapFlags);
        rebuild();
      });
    }

    if (resetBtn) {
      resetBtn.addEventListener("click", () => {
        activeListCmd = "";
        applySectionTab("model");
        const preservedWinFull = winFullPathEl ? winFullPathEl.checked : false;
        const preservedWinVer = els.wbRVersionSel ? els.wbRVersionSel.value : DEFAULTS.wbRVersionSel;
        const preservedWinCustom = els.wbRVersionCustom ? els.wbRVersionCustom.value : DEFAULTS.wbRVersionCustom;
        isHydrating = true;
        if (els.wbPostDir2) delete els.wbPostDir2.dataset.multiAdded;
        applyState(DEFAULTS);
        if (winFullPathEl) winFullPathEl.checked = preservedWinFull;
        if (els.wbRVersionSel) els.wbRVersionSel.value = preservedWinVer;
        if (els.wbRVersionCustom) els.wbRVersionCustom.value = preservedWinCustom;
        if (scannedDirs.length) populateScannedOutputs(scannedDirs);
        isHydrating = false;
        rebuild();
      });
    }

    function applySectionTab(tab) {
      activeSectionTab = tab || "all";
      document.querySelectorAll("#wbSectionTabs [data-wb-tab]").forEach((btn) => {
        const isActive = btn.dataset.wbTab === activeSectionTab;
        btn.classList.toggle("active", isActive);
        btn.setAttribute("aria-selected", isActive ? "true" : "false");
      });
      if (discoveryCardEl) {
        discoveryCardEl.hidden = activeSectionTab !== "all" && activeSectionTab !== "discovery";
      }
      if (configureSectionEl) {
        configureSectionEl.hidden = activeSectionTab === "discovery";
      }
      document.querySelectorAll(".wb-group-card[data-wb-group]").forEach((card) => {
        const grp = card.dataset.wbGroup || "";
        card.hidden = activeSectionTab !== "all" && grp !== activeSectionTab;
      });
    }

    document.querySelectorAll("#wbSectionTabs [data-wb-tab]").forEach((btn) => {
      btn.addEventListener("click", () => {
        applySectionTab(btn.dataset.wbTab || "all");
        syncHashState();
      });
    });

    const dirPickerEl = document.getElementById("wbOutputDirPicker");
    const dirDatalistEl = document.getElementById("wbOutputDirsList");
    const scanPostBtn = document.getElementById("wbScanOutputsPostBtn");
    const scanPlotBtn = document.getElementById("wbScanOutputsPlotBtn");
    const scanStatusPostEl = document.getElementById("wbScanStatusPost");
    const scanStatusPlotEl = document.getElementById("wbScanStatusPlot");
    const postDir1Sel = document.getElementById("wbPostDir1Sel");
    const postDir2Sel = document.getElementById("wbPostDir2Sel");
    const plotDirSel = document.getElementById("wbPlotDirSel");

    function extractOutputRunDirs(fileList) {
      const dirs = new Set();
      Array.from(fileList || []).forEach((file) => {
        const rel = (file.webkitRelativePath || "").replace(/\\/g, "/");
        if (!rel) return;
        const parts = rel.split("/").filter(Boolean);
        if (parts.length < 2) return;
        const rootName = parts[0].toLowerCase();
        const fileName = parts[parts.length - 1].toLowerCase();
        if (!/\.(csv|rds|png|svg|pdf|md)$/i.test(fileName)) return;

        // Normalize to repo-relative path under output/results/<run> where possible
        const resultsIdx = parts.findIndex((p) => p.toLowerCase() === "results");
        if (resultsIdx >= 0 && parts.length > resultsIdx + 2) {
          const runFolder = parts[resultsIdx + 1];
          if (runFolder && !runFolder.startsWith(".")) {
            dirs.add(`output/results/${runFolder}`);
          }
          const subFolder = parts[resultsIdx + 2];
          if (/^(season|year|month)_/i.test(subFolder) && parts.length > resultsIdx + 3) {
            dirs.add(`output/results/${runFolder}/${subFolder}`);
          }
        } else if (rootName === "results" && parts.length >= 3) {
          dirs.add(`output/results/${parts[1]}`);
        } else if (parts.length >= 3 && !parts[1].startsWith(".")) {
          dirs.add(`${parts[0]}/${parts[1]}`);
        } else if (parts.length === 2) {
          dirs.add(parts[0]);
        }
      });
      return Array.from(dirs).sort();
    }

    function populateScannedOutputs(dirs) {
      scannedDirs = Array.isArray(dirs) ? dirs : [];
      if (postDir2SelRow) postDir2SelRow.hidden = !scannedDirs.length;
      if (dirDatalistEl) {
        dirDatalistEl.innerHTML = "";
        dirs.forEach((d) => {
          const opt = document.createElement("option");
          opt.value = d;
          dirDatalistEl.appendChild(opt);
        });
      }

      const selectConfigs = [
        [postDir1Sel, els.wbPostDir1, "Pick primary run directory...", "Or write a custom directory here..."],
        [postDir2Sel, els.wbPostDir2, "Pick additional run directory...", "Or write custom directory(s) here..."],
        [plotDirSel, els.wbPlotDir, "Pick run directory to plot...", "Or write a custom directory here..."]
      ];
      selectConfigs.forEach(([selEl, inpEl, placeholder, customHint]) => {
        if (!selEl) return;
        selEl.innerHTML = "";
        if (!dirs.length) {
          selEl.hidden = true;
          if (inpEl) inpEl.classList.remove("wb-input-muted");
          return;
        }
        const ph = document.createElement("option");
        ph.value = "";
        ph.textContent = placeholder;
        selEl.appendChild(ph);
        dirs.forEach((d) => {
          const opt = document.createElement("option");
          opt.value = d;
          opt.textContent = d;
          selEl.appendChild(opt);
        });
        if (inpEl) {
          const cur = inpEl.value.trim();
          const fallbackDir = inpEl === els.wbPostDir2 && dirs.length > 1 ? dirs[1] : dirs[0];
          const matched = dirs.includes(cur) ? cur : fallbackDir;
          if (dirs.includes(cur) || cur === DEFAULTS[inpEl.id]) {
            selEl.value = matched;
            inpEl.value = "";
            inpEl.placeholder = customHint;
            inpEl.classList.add("wb-input-muted");
          } else if (cur) {
            selEl.value = "";
            inpEl.classList.remove("wb-input-muted");
          } else {
            selEl.value = fallbackDir;
            inpEl.placeholder = customHint;
            inpEl.classList.add("wb-input-muted");
          }
        }
        selEl.hidden = false;
      });

      const statusMsg = dirs.length ? `${dirs.length} folder${dirs.length === 1 ? "" : "s"} found` : "No run folders found";
      if (scanStatusPostEl) scanStatusPostEl.textContent = statusMsg;
      if (scanStatusPlotEl) scanStatusPlotEl.textContent = statusMsg;
    }

    async function scanDirectoryHandle(rootHandle) {
      const dirs = new Set();
      const rootName = (rootHandle.name || "").toLowerCase();

      async function collectSubdirs(parentHandle, prefix, includeTemporalChildren) {
        for await (const [name, handle] of parentHandle.entries()) {
          if (handle.kind !== "directory" || name.startsWith(".")) continue;
          const full = prefix ? `${prefix}/${name}` : name;
          dirs.add(full);
          if (includeTemporalChildren) {
            for await (const [subName, subHandle] of handle.entries()) {
              if (subHandle.kind === "directory" && /^(season|year|month)_/i.test(subName)) {
                dirs.add(`${full}/${subName}`);
              }
            }
          }
        }
      }

      if (rootName === "results") {
        await collectSubdirs(rootHandle, "output/results", true);
      } else {
        let foundResults = false;
        for await (const [name, handle] of rootHandle.entries()) {
          if (handle.kind !== "directory" || name.startsWith(".")) continue;
          if (name.toLowerCase() === "results") {
            foundResults = true;
            await collectSubdirs(handle, "output/results", true);
          } else if (name.toLowerCase() === "output") {
            for await (const [subName, subHandle] of handle.entries()) {
              if (subHandle.kind === "directory" && subName.toLowerCase() === "results") {
                foundResults = true;
                await collectSubdirs(subHandle, "output/results", true);
              }
            }
          }
        }
        if (!foundResults) {
          const basePrefix = rootName === "output" ? "output" : (rootHandle.name || "output");
          await collectSubdirs(rootHandle, basePrefix, true);
        }
      }
      return Array.from(dirs).sort();
    }

    [scanPostBtn, scanPlotBtn].forEach((btn) => {
      if (!btn) return;
      btn.addEventListener("click", async () => {
        if (typeof window.showDirectoryPicker === "function") {
          try {
            const handle = await window.showDirectoryPicker({ mode: "read" });
            const foundDirs = await scanDirectoryHandle(handle);
            populateScannedOutputs(foundDirs);
            return;
          } catch (err) {
            if (err && err.name === "AbortError") return;
          }
        }
        if (dirPickerEl) {
          dirPickerEl.value = "";
          dirPickerEl.click();
        }
      });
    });

    if (dirPickerEl) {
      dirPickerEl.addEventListener("change", () => {
        const foundDirs = extractOutputRunDirs(dirPickerEl.files);
        populateScannedOutputs(foundDirs);
      });
    }

    [
      [postDir1Sel, els.wbPostDir1, "stats"],
      [postDir2Sel, els.wbPostDir2, "compare"],
      [plotDirSel, els.wbPlotDir, null]
    ].forEach(([selEl, inpEl, targetMode]) => {
      if (!selEl || !inpEl) return;
      selEl.addEventListener("change", () => {
        if (!selEl.value) return;
        const multiTokens = inpEl.value.trim().split(/[\s,]+/).filter(Boolean);
        if (inpEl === els.wbPostDir2 && (multiTokens.length > 1 || inpEl.dataset.multiAdded === "1")) {
          return;
        }
        inpEl.value = "";
        inpEl.classList.add("wb-input-muted");
        activeListCmd = "";
        if (targetMode && els.wbMode) {
          if (targetMode === "stats" && els.wbMode.value === "compare") {
            // keep compare mode when picking primary run for comparison
          } else {
            els.wbMode.value = targetMode;
          }
        }
        rebuild();
      });
      inpEl.addEventListener("input", () => {
        const hasVal = Boolean(inpEl.value.trim());
        if (!hasVal) delete inpEl.dataset.multiAdded;
        inpEl.classList.toggle("wb-input-muted", !selEl.hidden && !hasVal);
        if (hasVal && !selEl.hidden) selEl.value = "";
      });
    });

    if (postDir2AddBtn && postDir2Sel && els.wbPostDir2) {
      postDir2AddBtn.addEventListener("click", () => {
        const chosen = postDir2Sel.value;
        if (!chosen) return;
        const currentTokens = els.wbPostDir2.value
          .split(/[\s,]+/)
          .map((s) => s.trim())
          .filter(Boolean);
        if (!currentTokens.includes(chosen)) {
          currentTokens.push(chosen);
        }
        els.wbPostDir2.value = currentTokens.join(" ");
        els.wbPostDir2.dataset.multiAdded = "1";
        els.wbPostDir2.classList.remove("wb-input-muted");
        activeListCmd = "";
        if (els.wbCompareMode) els.wbCompareMode.value = "dirs";
        if (els.wbMode) els.wbMode.value = "compare";
        rebuild();
      });
    }

    loadWinRPrefs();
    renderHistory();
    isHydrating = true;
    applyState(DEFAULTS);
    isHydrating = false;
    hydrateFromHash();
    applySectionTab(activeSectionTab);
    updateShellButtons();
    rebuild();
  }

  initThemeToggle();
  initInfoTooltips();
  initCopyButtons();
  initWorkbench();
})();
