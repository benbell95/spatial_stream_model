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
      "wbCvMode", "wbThresholdSel", "wbThresholdVal", "wbThresholdDir",
      "wbTemporal", "wbDesignSites", "wbBasinSel", "wbBasinCustom",
      "wbPreferSel", "wbPreferCustom", "wbWeightSel", "wbWeightCustom",
      "wbSeSel", "wbSeCustom", "wbFormat", "wbPalette", "wbOmp",
      "wbRVersionSel", "wbRVersionCustom"
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
    if (!outEl) return;

    let wrapFlags = true;
    let shellMode = "bash";
    let activeListCmd = "";

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
      wbCvMode: "random",
      wbThresholdSel: "none",
      wbThresholdVal: "10",
      wbThresholdDir: "above",
      wbTemporal: "all",
      wbDesignSites: "0",
      wbBasinSel: "none",
      wbBasinCustom: "",
      wbPreferSel: "",
      wbPreferCustom: "",
      wbWeightSel: "",
      wbWeightCustom: "",
      wbSeSel: "",
      wbSeCustom: "",
      wbFormat: "png_300",
      wbPalette: "auto",
      wbOmp: "4",
      wbRVersionSel: "4.4.3",
      wbRVersionCustom: ""
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
        [els.wbRVersionSel, els.wbRVersionCustom]
      ];
      pairs.forEach(([sel, inp]) => {
        if (!sel || !inp) return;
        const show = sel.value === "custom";
        inp.hidden = !show;
      });
      if (els.wbThresholdSel && thresholdWrap) {
        thresholdWrap.hidden = els.wbThresholdSel.value !== "custom";
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
        if (g.dataset === "multi") return datasetFilter === "auto" || mode === "cokrig";
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
      const effectiveDetDataset = dataset !== "auto"
        ? dataset
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
      const fmtVal = els.wbFormat ? els.wbFormat.value : "png_300";
      const palette = els.wbPalette ? els.wbPalette.value : "auto";
      const omp = els.wbOmp ? els.wbOmp.value : "4";

      const args = [];
      const notes = [];
      let runner = "Rscript";
      let scriptRel = "scripts/run_s3n_determinand.R";

      if (mode === "prep") {
        runner = shellMode === "bash" ? "python3" : "python";
        scriptRel = "scripts/01_preprocess_orn.py";
        if (source && source !== "rsn") args.push(`--source=${source}`);
        if (dataset && dataset !== "auto" && dataset !== "wq") args.push(`--dataset=${dataset}`);
        if (network && network !== "orn") args.push(`--network=${network}`);
        if (transform && transform !== "default") args.push(`--transform=${transform}`);
        if (censor && censor !== "default") args.push(`--censor=${censor}`);
        if (source.toLowerCase().includes("waterblitz")) {
          args.push("--max-snap-dist=350");
          notes.push("Citizen-science WaterBlitz sources automatically benefit from <code>--max-snap-dist=350</code> snapping tolerance.");
        }
        notes.push("Writes <code>data/processed/orn_snapped_*.rds</code>, topology weights, and sub-catchment AFVs without R/GDAL spatial dependencies.");
      } else if (mode === "cokrig") {
        runner = "Rscript";
        scriptRel = "scripts/run_s3n_cokriging.R";
        const coTargets = target.includes(",") ? target : `${target},0180`;
        args.push(`--targets=${coTargets}`);
        if (dataset && dataset !== "auto" && dataset !== "wq") args.push(`--dataset=${dataset}`);
        if (source && source !== "rsn") args.push(`--source=${source}`);
        if (network && network !== "orn") args.push(`--network=${network}`);
        if (trend === "colocated") args.push("--method=colocated");
        else args.push("--method=lmc");
        if (kernel && kernel !== "exponential") args.push(`--kernel=${kernel}`);
        if (afv && afv !== "length") args.push(`--afv=${afv}`);
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
        else args.push(`--dpi=${fmtVal.replace("png_", "")}`);
        args.push(`--omp=${omp}`);

        notes.push("Fits Multivariate Factor-LMC or Co-Located Co-Kriging with Dual CV (cross-target univariate vs. joint gain).");
        if (transform !== "default") {
          notes.push("Note: Co-Kriging applies each target's configured transform in <code>configs/determinands.json</code> automatically.");
        }
      } else {
        runner = "Rscript";
        scriptRel = "scripts/run_s3n_determinand.R";
        const singleTarget = target.includes(",") ? target.split(",")[0].trim() : target;
        args.push(singleTarget || "0117");
        if (effectiveDetDataset && effectiveDetDataset !== "auto") args.push(`--dataset=${effectiveDetDataset}`);
        if (source && source !== "rsn") args.push(`--source=${source}`);
        if (network && network !== "orn") args.push(`--network=${network}`);
        if (transform && transform !== "default") args.push(`--transform=${transform}`);
        if (norm && norm !== "none") args.push(`--norm=${norm}`);
        if (censor && censor !== "default") args.push(`--censor=${censor}`);
        if (trend && trend !== "linear") args.push(`--trend=${trend}`);
        if (kernel && kernel !== "exponential") args.push(`--kernel=${kernel}`);
        if (afv && afv !== "length") args.push(`--afv=${afv}`);
        if (cvMode && cvMode !== "random") args.push(`--cv-mode=${cvMode}`);

        if (thresholdMode === "uk") {
          args.push("--threshold=uk");
        } else if (thresholdMode === "custom" && thresholdVal !== "") {
          args.push(`--threshold=${thresholdVal}`);
          if (thresholdDir === "below") args.push("--threshold-dir=below");
        }

        if (temporal === "by_season") args.push("--by=season", "--st-borrow");
        else if (temporal === "by_year") args.push("--by=year", "--st-borrow");
        else if (temporal === "by_month") args.push("--by=month", "--st-borrow");
        else if (temporal === "season_winter") args.push("--season=winter");
        else if (temporal === "year_2024") args.push("--year=2024");

        if (designSites !== "0") args.push(`--design-sites=${designSites}`);
        if (basin && basin !== "none") args.push(`--basin=${basin}`);
        if (prefer) args.push(`--prefer-source=${prefer}`);
        if (weight) args.push(`--source-weight=${weight}`);
        if (sourceSe) args.push(`--source-se=${sourceSe}`);

        if (fmtVal === "svg") args.push("--fig-format=svg");
        else if (fmtVal === "pdf") args.push("--fig-format=pdf");
        else args.push(`--dpi=${fmtVal.replace("png_", "")}`);

        if (palette === "standard") args.push("--normal-color");
        else if (palette === "reverse") args.push("--reverse-color");

        args.push(`--omp=${omp}`);

        notes.push("Outputs full <code>90,876</code>-reach predictions (<code>s3n_predictions_*.rds</code>), CV table, parameter estimates, and publication maps.");
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
      const shellLabel = shellMode === "powershell" ? "PowerShell" : shellMode === "cmd" ? "Command Prompt" : "Bash";
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

      if (mode !== "prep") {
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
      }

      if (source === "ssn" && (isRsnOnlyTarget || effectiveDetDataset === "rhs" || effectiveDetDataset === "lcms")) {
        notes.push("Warning: The selected target/modality is monitored on <code>rsn</code> (River Surveillance Network) only and is unavailable in <code>--source=ssn</code> alone.");
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

    ids.forEach((id) => {
      if (!els[id]) return;
      const evt = els[id].tagName === "INPUT" ? "input" : "change";
      els[id].addEventListener(evt, () => {
        activeListCmd = "";
        if (id === "wbDataset" || id === "wbMode") {
          populateTargetSelect(els.wbDataset ? els.wbDataset.value : "auto", els.wbMode ? els.wbMode.value : "det");
        }
        rebuild();
      });
    });

    if (winFullPathEl) {
      winFullPathEl.addEventListener("change", rebuild);
    }

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
        if (winFullPathEl) winFullPathEl.checked = false;
        applyState(DEFAULTS);
      });
    }

    updateShellButtons();
    applyState(DEFAULTS);
  }

  initThemeToggle();
  initInfoTooltips();
  initCopyButtons();
  initWorkbench();
})();
