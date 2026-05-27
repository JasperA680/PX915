"""Cellular-automaton tab: preset selector, model toggle, network preview, run + plots."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QComboBox,
    QPushButton, QLabel, QScrollArea,
)

from python.road_network import PRESETS, NetworkSpec, LayoutSpec
from python.io import load_network_netcdf

from python.gui.param_form import ParamForm
from python.gui.network_widget import NetworkWidget
from python.gui.plot_panel import PlotPanel
from python.gui.road_edit_dialog import RoadEditDialog, JunctionEditDialog
from python.gui.runners import RunnerThread, FDSweepThread
from python.gui._common import set_button_stale, ALWAYS_VISIBLE_SCROLLBAR_QSS
from python.gui._run_tab import RunTab


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_BINARY = REPO_ROOT / "build" / "run_network"
DEFAULT_OUTDIR = REPO_ROOT / "data" / "output" / "gui"

# Available CA update rules.  Only "NS" is implemented physics-side; "TASEP"
# is plumbed through the JSON config but the Fortran driver errors out on it.
CA_MODELS = ("NS", "TASEP")

# Presets where TASEP physics is appropriate.  TASEP is a 1D single-chain
# process — any preset whose roads decompose into independent 1D lanes with
# right-ward hops is fair game.  That includes the junction/network presets,
# whose roads are still individual chains connected by junctions.  Only
# ``two_lane`` is excluded because the whole point of that preset is the
# lane-change feature, which TASEP doesn't model.
TASEP_OK_PRESETS = {"single_lane", "t_junction", "crossroads", "roundabout", "town"}

# Presets where the fundamental-diagram sweep is meaningful.  The sweep
# uses pure-Python TASEP on a single chain, so junctions / multilane / routing
# scenarios don't relate to its J(ρ) curve.
FD_OK_PRESETS = {"single_lane"}


class CATab(RunTab):
    """Cellular-automaton tab: preset selector, model toggle, network preview, run + plots."""

    def __init__(self, binary: Path, output_dir: Path, parent=None):
        super().__init__(parent)
        self._binary = binary
        self._output_dir = output_dir
        self._fd_runner: Optional[FDSweepThread] = None
        self._current_spec: Optional[NetworkSpec] = None
        self._current_layout: Optional[LayoutSpec] = None
        # Per-road parameter overrides set via the click-to-edit dialog.
        # Keyed by road_id; reset when the preset combo changes.
        self._road_overrides: dict = {}
        # Per-junction routing overrides (dict[junction_id -> list[list[float]]]).
        self._junction_overrides: dict = {}
        # FD-sweep tracking is independent of the main run.
        self._has_fd_results: bool = False

        # --- Toolbar row: preset + model + run ---
        toolbar = QWidget()
        toolbar_row = QHBoxLayout(toolbar)
        toolbar_row.addWidget(QLabel("Preset:"))
        self.preset_combo = QComboBox()
        for name in PRESETS:
            self.preset_combo.addItem(name)
        self.preset_combo.setToolTip(
            "Network topology to simulate.  Changing the preset rebuilds the\n"
            "network preview and resets any per-road or per-junction edits.\n"
            "  single_lane — one open chain\n"
            "  two_lane    — two parallel same-direction lanes\n"
            "  t_junction  — three-leg merge/diverge\n"
            "  crossroads  — four-leg intersection with routing\n"
            "  roundabout  — circular ring with arms\n"
            "  town        — small grid of crossroads"
        )
        toolbar_row.addWidget(self.preset_combo)
        toolbar_row.addSpacing(12)
        toolbar_row.addWidget(QLabel("Model:"))
        self.model_combo = QComboBox()
        for m in CA_MODELS:
            self.model_combo.addItem(m)
        self.model_combo.setToolTip(
            "Cellular-automaton update rule.\n"
            "  NS    — Nagel–Schreckenberg (multi-speed, dawdling probability)\n"
            "  TASEP — strict 1D single-lane chain (greyed out where it isn't\n"
            "          physically meaningful, e.g. multi-lane presets)."
        )
        toolbar_row.addWidget(self.model_combo)
        toolbar_row.addSpacing(12)
        self.run_button = QPushButton("Run")
        self.run_button.setMinimumWidth(120)
        toolbar_row.addWidget(self.run_button)
        toolbar_row.addSpacing(12)
        self.fd_button = QPushButton("Run FD sweep")
        self.fd_button.setMinimumWidth(140)
        self.fd_button.setToolTip(
            "Run a fundamental-diagram sweep using the currently selected\n"
            "Model (NS or TASEP). Populates the Fundamental Diagram plot tab.\n\n"
            "  • TASEP: open-chain α/β sweep (β=1 then α=1 branches).\n"
            "  • NS:    periodic-ring density sweep at the form's v_max / p_slow.\n\n"
            "Only meaningful for the single_lane preset — both sweeps run a\n"
            "homogeneous 1D system, so junctions / multilane / routing\n"
            "scenarios don't relate to the J(ρ) curve and the button is\n"
            "greyed out there."
        )
        toolbar_row.addWidget(self.fd_button)
        toolbar_row.addStretch(1)

        # --- Left: params (scrollable, capped height) + network preview (fills rest) ---
        # The form has ~10 controls (~450 px tall natural).  Wrapping it in a
        # QScrollArea with maximumHeight=340 lets the network widget claim the
        # rest of the vertical space without the form ever pushing it out.
        self.param_form = ParamForm()
        form_scroll = QScrollArea()
        form_scroll.setWidget(self.param_form)
        form_scroll.setWidgetResizable(True)
        form_scroll.setMaximumHeight(340)
        form_scroll.setFrameShape(QScrollArea.NoFrame)
        # Keep the vertical scrollbar visible at all times so users get a
        # static visual cue that the form has more content below the fold.
        # macOS auto-hides native scrollbars regardless of ScrollBarAlwaysOn;
        # applying *any* QScrollBar stylesheet switches Qt off the native
        # style and onto its own renderer, which honours AlwaysOn properly.
        form_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        form_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        form_scroll.setStyleSheet(ALWAYS_VISIBLE_SCROLLBAR_QSS)

        self.network_widget = NetworkWidget()
        self.network_widget.setMinimumSize(360, 360)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.addWidget(form_scroll)
        left_layout.addWidget(self.network_widget, stretch=1)

        # --- Right: plot panel ---
        self.plot_panel = PlotPanel()

        # --- Central layout ---
        body = QSplitter(Qt.Horizontal)
        body.addWidget(left)
        body.addWidget(self.plot_panel)
        body.setStretchFactor(0, 0)
        body.setStretchFactor(1, 1)
        body.setSizes([520, 800])

        outer = QVBoxLayout(self)
        outer.addWidget(toolbar)
        outer.addWidget(body, stretch=1)

        # --- Signals ---
        self.preset_combo.currentTextChanged.connect(self._on_preset_changed_ca)
        self.model_combo.currentTextChanged.connect(self._on_model_changed)
        self.param_form.params_changed.connect(self._rebuild_spec)
        self.param_form.params_changed.connect(self._mark_stale)
        self.run_button.clicked.connect(self._on_run)
        self.fd_button.clicked.connect(self._on_run_fd)
        self.network_widget.road_clicked.connect(self._on_road_clicked)
        self.network_widget.junction_clicked.connect(self._on_junction_clicked)

        # Bootstrap with the first preset (also greys out TASEP / FD if needed).
        self._sync_model_availability()
        self._sync_fd_availability()
        self._rebuild_spec()

    def _on_preset_changed_ca(self, name: str):
        """Reset per-road and per-junction overrides on preset change, then rebuild."""
        self._road_overrides = {}
        self._junction_overrides = {}
        self._sync_model_availability()
        self._sync_fd_availability()
        self._rebuild_spec()
        # Preset switch invalidates any displayed results.
        self._mark_stale()

    def _on_model_changed(self, _name: str):
        # Switching the model means the next Run would produce different
        # results; any currently-shown plots are stale from the new model's POV.
        self._mark_stale()

    def _mark_stale(self):
        """Mark the plot panel + buttons as out-of-date.

        Drives three indicators in concert:
          - corner banner on the tab bar
          - Run / Run-FD button styling (orange dot when their output is stale)
          - FD tab title (so the user notices it inside the panel too)
        """
        super()._mark_stale()
        set_button_stale(self.fd_button, "Run FD sweep", self._has_fd_results)
        self.plot_panel.set_fd_stale(self._has_fd_results)

    def _mark_fd_fresh(self):
        """Clear the stale indicator on the FD button + tab after a sweep."""
        self._has_fd_results = True
        set_button_stale(self.fd_button, "Run FD sweep", False)
        self.plot_panel.set_fd_stale(False)

    def _sync_fd_availability(self):
        """Grey out the Run-FD button + FD tab for presets where it isn't meaningful.

        The FD sweep is a pure-Python single-lane TASEP exercise, so junctions
        / multilane / routing presets don't relate to its J(ρ) curve.
        """
        name = self.preset_combo.currentText()
        ok = name in FD_OK_PRESETS
        self.fd_button.setEnabled(ok)
        if not ok:
            # Tab back to placeholder + disabled; any prior sweep result is wiped.
            self.plot_panel.reset_fd()
            self._has_fd_results = False
            set_button_stale(self.fd_button, "Run FD sweep", False)

    def _sync_model_availability(self):
        """Grey out TASEP for any preset where it is not physically meaningful.

        TASEP is a strict 1D single-lane chain — anything with multiple lanes
        or junctions falls outside its definition.
        """
        name = self.preset_combo.currentText()
        tasep_ok = name in TASEP_OK_PRESETS
        # Find the TASEP entry in the combo and toggle its enabled flag.
        idx = self.model_combo.findText("TASEP")
        if idx >= 0:
            model = self.model_combo.model()
            item = model.item(idx)
            if item is not None:
                # Qt.ItemIsEnabled = 32; remove or add it to control greying.
                if tasep_ok:
                    item.setFlags(item.flags() | Qt.ItemIsEnabled)
                else:
                    item.setFlags(item.flags() & ~Qt.ItemIsEnabled)
        # If TASEP is now disabled but currently selected, switch to NS.
        if not tasep_ok and self.model_combo.currentText() == "TASEP":
            ns_idx = self.model_combo.findText("NS")
            if ns_idx >= 0:
                self.model_combo.setCurrentIndex(ns_idx)

    def _rebuild_spec(self, *_args):
        name = self.preset_combo.currentText()
        if name not in PRESETS:
            return
        kwargs = self.param_form.preset_kwargs(name)
        try:
            spec, layout = PRESETS[name](**kwargs)
        except TypeError as exc:
            self.log.emit(f"preset error: {exc}")
            return
        # Apply any per-road overrides set via the click-to-edit dialog.
        for rid, ov in self._road_overrides.items():
            road = next((r for r in spec.roads if r.id == rid), None)
            if road is None:
                continue
            for ln in road.lanes:
                if ov.get("alpha") is not None and getattr(ln, "open_in", False):
                    ln.alpha = ov["alpha"]
                if ov.get("beta") is not None and getattr(ln, "open_out", False):
                    ln.beta = ov["beta"]
                if ov.get("length") is not None:
                    ln.length = int(ov["length"])
        # Apply any per-junction routing overrides.
        for jid, routes in self._junction_overrides.items():
            j = next((j for j in spec.junctions if j.id == jid), None)
            if j is None:
                continue
            j.routes = [list(row) for row in routes]
        self._current_spec = spec
        self._current_layout = layout
        self.network_widget.set_network(spec, layout)
        # Enable lane-change UI only when the spec has at least one road with
        # multiple same-direction lanes (otherwise the parameters do nothing).
        has_parallel = any(
            len(r.lanes) >= 2
            and len({ln.flow_direction for ln in r.lanes}) == 1
            for r in spec.roads
        )
        self.param_form.set_lane_change_enabled(has_parallel)
        self.status.emit(
            f"preset {name}: {len(spec.roads)} roads, {len(spec.junctions)} junctions"
        )

    def _on_road_clicked(self, road_id: int):
        """Open the per-road editor dialog and apply any changes."""
        if self._current_spec is None:
            return
        road = next((r for r in self._current_spec.roads if r.id == road_id), None)
        if road is None:
            return

        # Read current α/β/length from the road's lanes.  Length lives on
        # ``LaneSpec`` (RoadSpec has no length attribute), so we read it
        # from the first lane — the click-to-edit dialog applies the same
        # length to every lane on the road in _rebuild_spec.
        alpha = next((ln.alpha for ln in road.lanes if getattr(ln, "open_in", False)), 0.5)
        beta  = next((ln.beta  for ln in road.lanes if getattr(ln, "open_out", False)), 0.5)
        has_alpha = any(getattr(ln, "open_in",  False) for ln in road.lanes)
        has_beta  = any(getattr(ln, "open_out", False) for ln in road.lanes)
        length = road.lanes[0].length if road.lanes else 20

        dlg = RoadEditDialog(road_id, alpha, beta, length,
                             has_alpha, has_beta, parent=self)
        if dlg.exec_() == RoadEditDialog.Accepted:
            vals = dlg.values()
            self._road_overrides[road_id] = vals
            self._rebuild_spec()
            self._mark_stale()
            self.log.emit(
                f"R{road_id} overridden: α={vals['alpha']}  β={vals['beta']}  L={vals['length']}"
            )

    def _on_junction_clicked(self, junction_id: int):
        """Open the junction routing editor and apply any changes."""
        if self._current_spec is None:
            return
        j = next((j for j in self._current_spec.junctions if j.id == junction_id), None)
        if j is None or not j.routes:
            self.log.emit(f"J{junction_id}: no routing matrix to edit")
            return

        in_road_ids  = [int(leg.road) for leg in j.in_legs]
        out_road_ids = [int(leg.road) for leg in j.out_legs]
        dlg = JunctionEditDialog(
            junction_id, j.routes,
            in_road_ids=in_road_ids, out_road_ids=out_road_ids,
            parent=self,
        )
        if dlg.exec_() == JunctionEditDialog.Accepted:
            new_routes = dlg.routes()
            self._junction_overrides[junction_id] = new_routes
            self._rebuild_spec()
            self._mark_stale()
            summary = "  ".join(
                "[" + ", ".join(f"{v:.2f}" for v in row) + "]"
                for row in new_routes
            )
            self.log.emit(f"J{junction_id} routes overridden: {summary}")

    def _on_run(self):
        if self._current_spec is None or self._current_layout is None:
            self.log.emit("no spec to run")
            return
        if not self._binary.exists():
            self.log.emit(f"binary not found: {self._binary} — try `make run_network`")
            return
        if self._runner is not None and self._runner.isRunning():
            return

        params = self.param_form.sim_params()
        params.model = self.model_combo.currentText()
        preset_name = self.preset_combo.currentText()
        out_dir = self._output_dir / preset_name
        self.run_button.setEnabled(False)
        self.progress.emit(0)
        self.status.emit(f"running {preset_name} model={params.model} ({params.n_steps} steps)...")
        self.log.emit(f"--- run {preset_name} model={params.model} ---")

        self._runner = RunnerThread(
            self._current_spec, params, self._current_layout,
            output_dir=out_dir, binary=self._binary, parent=self,
        )
        self._runner.progress.connect(self._on_runner_progress)
        self._runner.finished_ok.connect(self._on_finished_ok)
        self._runner.finished_error.connect(self._on_finished_error)
        self._mark_running(f"{preset_name} / {params.model}")
        self._runner.start()

    # ------------------------------------------------------------------
    # Fundamental-diagram sweep (pure-Python TASEP)
    # ------------------------------------------------------------------

    def _on_run_fd(self):
        if self._fd_runner is not None and self._fd_runner.isRunning():
            return
        model = self.model_combo.currentText()   # "NS" or "TASEP"
        # Pull the road length from the current spec — that way the sweep
        # uses whatever the user configured (preset default or per-road
        # click-to-edit override).  FD button is only enabled for single_lane
        # so the first road's first lane is the canonical length.
        if (self._current_spec is not None
                and self._current_spec.roads
                and self._current_spec.roads[0].lanes):
            L = int(self._current_spec.roads[0].lanes[0].length)
        else:
            L = 100 if model == "NS" else 50
        # TASEP runs an alpha- and a beta-branch, so n_points = 30 yields 60
        # output rows; bump NS to 60 so both diagrams have the same point
        # density. NS rounds each target rho to an integer vehicle count, so
        # cap at L to avoid two adjacent points collapsing onto the same
        # n_vehicles on short rings.
        n_points = min(60, L) if model == "NS" else 30
        v_max = int(self.param_form.v_max.value())
        p_slow = float(self.param_form.p_slow.value())
        self.fd_button.setEnabled(False)
        if model == "NS":
            self.status.emit(
                f"running NS FD sweep L={L} v_max={v_max} p_slow={p_slow} "
                f"({n_points} points)…"
            )
            self.log.emit(
                f"--- NS FD sweep L={L} v_max={v_max} p_slow={p_slow} "
                f"points={n_points} ---"
            )
        else:
            self.status.emit(f"running TASEP FD sweep L={L} ({2 * n_points} points)…")
            self.log.emit(f"--- TASEP FD sweep L={L} points={n_points} ---")
        self._fd_runner = FDSweepThread(
            model=model, L=L, n_points=n_points,
            v_max=v_max, p_slow=p_slow, parent=self,
        )
        self._fd_runner.finished_ok.connect(self._on_fd_ok)
        self._fd_runner.finished_error.connect(self._on_fd_error)
        self._fd_runner.start()

    def _on_fd_ok(self, rho, J, model: str, info: dict):
        self.fd_button.setEnabled(True)
        self.status.emit(f"{model} FD sweep done — {len(rho)} points")
        self.log.emit(
            f"{model} FD sweep complete: L={info['L']}, {len(rho)} points  "
            f"(ρ range {float(rho.min()):.2f}–{float(rho.max()):.2f}, "
            f"J range {float(J.min()):.3f}–{float(J.max()):.3f})"
        )
        title = (
            f"Fundamental diagram ({model}, L={info['L']}, "
            f"v_max={info['v_max']})" if model == "NS"
            else f"Fundamental diagram ({model}, L={info['L']})"
        )
        self.plot_panel.set_fd_sweep(
            rho, J, title=title,
            model=model, v_max=int(info['v_max']),
        )
        self._mark_fd_fresh()

    def _on_fd_error(self, err: str):
        self.fd_button.setEnabled(True)
        self.status.emit("FD sweep error")
        self.log.emit(f"FD ERROR: {err}")

    # ------------------------------------------------------------------

    def _on_finished_ok(self, nc_path: str):
        self.progress.emit(100)
        self.status.emit(f"done — {nc_path}")
        self.run_button.setEnabled(True)
        self.log.emit(f"loaded {nc_path}")
        try:
            result = load_network_netcdf(nc_path)
            self.plot_panel.show_result(result)
            self._log_diagnostics(result)
            preset = self.preset_combo.currentText()
            model = self.model_combo.currentText()
            params = self.param_form.sim_params()
            self._mark_fresh(f"{preset} / {model} / {params.n_steps} steps")
        except Exception as exc:
            self.log.emit(f"plot error: {exc}")
            self._mark_stale()

    def _log_diagnostics(self, result):
        """Print per-road steady-state ρ̄ and J̄ to the run log."""
        import numpy as np
        n_steps, n_roads = result.road_density.shape
        # Average over the second half of the run to skip the transient.
        half = max(1, n_steps // 2)
        lane_road = np.asarray(result.lane_road_id)

        self.log.emit("--- diagnostics (steady-state averages, last half of run) ---")
        for r in range(n_roads):
            rho_bar = float(result.road_density[half:, r].mean())
            n_lanes = max(1, int(np.sum(lane_road == r + 1)))
            J_bar = float(result.road_exits[half:, r].mean()) / n_lanes
            total_in  = int(result.road_entries[:, r].sum())
            total_out = int(result.road_exits[:, r].sum())
            self.log.emit(
                f"  R{r + 1}: <ρ>={rho_bar:.3f}  <J>={J_bar:.3f}  "
                f"(in={total_in}, out={total_out}, lanes={n_lanes})"
            )

    def stop_running(self):
        super().stop_running()
        if self._fd_runner is not None and self._fd_runner.isRunning():
            self._fd_runner.requestInterruption()
            self._fd_runner.wait(2000)
