"""Run and visualise multi-lane LWR PDE simulations.

Five scenarios:

  A) Independent lanes (k=0)    — staggered IC, lanes stay different forever
  B) Coupled lanes (k=0.5)      — same IC, lanes equilibrate via lane changing
  C) Conservation check         — periodic sine IC, total mass must stay flat
  D) Fast/slow lanes (k=0.5)    — lane 1 v_max=1.0, lane 2 v_max=1.5;
                                   slow lane ends up denser, fast lane has higher flow
  E) Speed limit + fast/slow    — v_limit=0.6 caps both lanes; free-flow Δv vanishes
                                   so lane-changing is suppressed until congestion builds

Usage
-----
    python scripts/run_multilane_pde.py           # build, run all five, display
    python scripts/run_multilane_pde.py --save    # also save to plots/
    python scripts/run_multilane_pde.py --no-run  # re-plot existing .nc files
    python scripts/run_multilane_pde.py --scenario E  # single scenario
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "src" / "python"))

import matplotlib.pyplot as plt
from pde_runner import run_pde, load_pde_netcdf
from visualisation import plot_lane_densities, plot_total_mass
from analysis import compute_total_mass

PLOTS_DIR = ROOT / "plots"

SCENARIOS = {
    # Lane 1 starts dense (ρ=0.7, slow), lane 2 starts sparse (ρ=0.2, fast).
    # With k=0 there is no coupling so they stay different forever.
    "A_independent": dict(
        n_lanes=2, lane_change_rate=0.0,
        v_max=1.0, rho_max=1.0,
        rho_left_bc=0.7, rho_right_bc=0.2,
        ic_type="staggered", flux_type="godunov", bc_type="periodic",
        M=200, n_steps=800,
        label="Independent lanes  (k=0)  —  stay different",
    ),
    # Same starting point as A, but k=0.5 couples the lanes.
    # Vehicles move from the slow dense lane into the fast sparse lane
    # until the lanes equilibrate toward equal density.
    "B_coupled": dict(
        n_lanes=2, lane_change_rate=0.5,
        v_max=1.0, rho_max=1.0,
        rho_left_bc=0.7, rho_right_bc=0.2,
        ic_type="staggered", flux_type="godunov", bc_type="periodic",
        M=200, n_steps=800,
        label="Lane changing  (k=0.5)  —  lanes equilibrate",
    ),
    # Periodic sine wave: both lanes start identically, k=0.5.
    # Total mass must stay flat — the conservation diagnostic.
    "C_conservation": dict(
        n_lanes=2, lane_change_rate=0.5,
        v_max=1.0, rho_max=1.0,
        rho_left_bc=0.3, rho_right_bc=0.7,
        ic_type="sine", flux_type="godunov", bc_type="periodic",
        M=200, n_steps=800,
        label="Conservation check  (k=0.5, periodic sine IC)",
    ),
    # Fast/slow lanes: lane 1 v_max=1.0 (slow), lane 2 v_max=1.5 (fast).
    # Both lanes start at the same uniform density (rho=0.45).
    # Under periodic BCs, total mass is conserved.
    #
    # Physics: at equal density, lane 2 is always faster, so the source
    # term moves vehicles FROM lane 1 INTO lane 2.  Equilibrium is reached
    # when velocities equalise: v1(rho1) = v2(rho2).
    # With total mass 0.9:  rho1_eq ~ 0.34,  rho2_eq ~ 0.56.
    #
    # Expected outcome (confirmed by Phase 5 simulation):
    #   - Fast lane (lane 2) ends up DENSER  (rho2 > rho1)
    #   - Fast lane carries HIGHER flow      (q2 > q1)
    #   - Slow lane ends up SPARSER         (vehicles drained into fast lane)
    "D_fast_slow": dict(
        n_lanes=2, lane_change_rate=0.5,
        v_max=1.0, rho_max=1.0,
        rho_left_bc=0.45, rho_right_bc=0.45,
        v_max_lanes=[1.0, 1.5],
        ic_type="constant", flux_type="godunov", bc_type="periodic",
        M=200, n_steps=1200,
        label="Fast/slow lanes  (v=[1.0,1.5], k=0.5)  —  fast lane fills up",
    ),
    # Speed limit + fast/slow lanes.
    # v_limit=0.6, v_max_lane2=1.5.  rho* = rho_max*(1 - v_limit/v_max) = 0.4.
    # Below rho*: both lanes travel at v_limit → Δv = 0 → no lane-change incentive.
    # Above rho*: speed limit inactive → normal asymmetric exchange resumes.
    # Staggered IC (lane1 dense, lane2 sparse) so the contrast is visible early.
    "E_speed_limit": dict(
        n_lanes=2, lane_change_rate=0.5,
        v_max=1.0, rho_max=1.0,
        v_limit=0.6,
        v_max_lanes=[1.0, 1.5],
        rho_left_bc=0.2, rho_right_bc=0.7,
        ic_type="staggered", flux_type="godunov", bc_type="periodic",
        M=200, n_steps=1200,
        label="Speed limit  (v_limit=0.6, v=[1.0,1.5], k=0.5)\n"
              "—  limit suppresses switching in free-flow zone (ρ < ρ*=0.4)",
    ),
}


def run_scenario(key, params, nc_path):
    p = {k: v for k, v in params.items() if k != "label"}
    run_pde(p, output_path=nc_path, exe=ROOT / "build" / "pde_solver")
    print(f"  [{key}] done -> {nc_path.name}")


def plot_scenario(key, label, data, save_dir=None):
    n_lanes = data["n_lanes"]

    fig = plt.figure(figsize=(16, 10))
    fig.suptitle(label, fontsize=13, fontweight="bold")

    gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.38)

    # Top row: per-lane heatmaps (up to 2)
    ax_l1 = fig.add_subplot(gs[0, 0])
    ax_l2 = fig.add_subplot(gs[0, 1])
    ax_tot = fig.add_subplot(gs[0, 2])

    # Bottom row: per-lane density at x=0.5, total mass, total space-time
    ax_xt  = fig.add_subplot(gs[1, 0])
    ax_mass = fig.add_subplot(gs[1, 1])
    ax_snap = fig.add_subplot(gs[1, 2])

    x    = data["x"]
    time = data["time"]
    rho_pl = data["density_per_lane"]    # (time, lane, x)
    attrs  = data["attrs"]
    rho_max = float(attrs.get("rho_max", 1.0))

    # Per-lane heatmaps
    v_limit = attrs.get("v_limit", None)
    v_max_g = float(attrs.get("v_max", 1.0))
    rho_star = rho_max * (1.0 - float(v_limit) / v_max_g) if v_limit is not None else None

    for lane_idx, ax in enumerate([ax_l1, ax_l2]):
        if lane_idx < n_lanes:
            im = ax.imshow(
                rho_pl[:, lane_idx, :],
                aspect="auto", origin="lower", cmap="viridis",
                vmin=0, vmax=rho_max,
                extent=[float(x[0]), float(x[-1]), float(time[0]), float(time[-1])],
            )
            plt.colorbar(im, ax=ax, label="ρ")
            ax.set_xlabel("x"); ax.set_ylabel("t")
            title = f"Lane {lane_idx + 1}  ρ(x, t)"
            if rho_star is not None:
                title += f"\n(v_limit={float(v_limit):.2f}, ρ*={rho_star:.2f})"
            ax.set_title(title, fontsize=9)
        else:
            ax.set_visible(False)

    # Total space-time
    rho_tot = rho_pl.sum(axis=1)
    im2 = ax_tot.imshow(
        rho_tot, aspect="auto", origin="lower", cmap="viridis",
        vmin=0, vmax=rho_max * n_lanes,
        extent=[float(x[0]), float(x[-1]), float(time[0]), float(time[-1])],
    )
    plt.colorbar(im2, ax=ax_tot, label="ρ_tot")
    ax_tot.set_xlabel("x"); ax_tot.set_ylabel("t")
    ax_tot.set_title("Total  ρ_tot(x, t)")

    # Per-lane density vs time at x=0.5
    plot_lane_densities(data, x_pos=0.5, ax=ax_xt)

    # Total mass vs time
    plot_total_mass(data, ax=ax_mass)

    # Final snapshot
    colors = plt.cm.tab10([0.0, 0.1])
    for lane_idx in range(n_lanes):
        ax_snap.plot(x, rho_pl[-1, lane_idx, :],
                     color=colors[lane_idx], linewidth=1.4,
                     label=f"Lane {lane_idx + 1} (final)")
    ax_snap.set_xlabel("x"); ax_snap.set_ylabel("ρ")
    ax_snap.set_title("Final density snapshot")
    ax_snap.set_ylim(-0.02, rho_max * 1.1)
    ax_snap.legend(fontsize=9)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    if save_dir:
        save_dir.mkdir(parents=True, exist_ok=True)
        path = save_dir / f"multilane_{key}.png"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"  Saved {path}")

    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--no-run", action="store_true",
                        help="Skip build and run; replot existing .nc files")
    parser.add_argument("--save", action="store_true",
                        help="Save figures to plots/")
    parser.add_argument("--scenario", metavar="KEY",
                        help="Run only one scenario, e.g. D or D_fast_slow")
    args = parser.parse_args()

    if not args.no_run:
        print("Building PDE solver...")
        subprocess.run(["make", "pde"], check=True, cwd=ROOT, capture_output=True)

    save_dir = PLOTS_DIR if args.save else None

    # Filter to a single scenario if requested
    scenarios = SCENARIOS
    if args.scenario:
        key_filter = args.scenario.upper()
        scenarios = {k: v for k, v in SCENARIOS.items()
                     if k.upper().startswith(key_filter)}
        if not scenarios:
            print(f"No scenario matching '{args.scenario}'. Available: {list(SCENARIOS)}")
            return

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        nc_files = {key: tmp / f"{key}.nc" for key in scenarios}

        if not args.no_run:
            print("Running scenarios...")
            for key, params in scenarios.items():
                run_scenario(key, params, nc_files[key])

        print("Plotting...")
        for key, params in scenarios.items():
            data = load_pde_netcdf(nc_files[key])
            plot_scenario(key, params["label"], data, save_dir=save_dir)

        plt.show()


if __name__ == "__main__":
    main()
