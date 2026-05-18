"""End-to-end Python pipeline test.

For each preset:
  * build & validate the spec
  * serialise to JSON, invoke ./build/run_network
  * load the NetCDF result
  * assert shape, schema, mass conservation, and config roundtrip

Run from repo root:
    .venv/bin/python tests/test_run_pipeline.py
"""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from python.road_network import PRESETS, SimParams, validate, spec_to_dict
from python.run_simulation import run_simulation
from python.io import load_network_netcdf


BINARY = ROOT / "build" / "run_network"


def assert_pipeline(preset_name: str, tmp: Path):
    print(f"--- {preset_name} ---")
    spec, layout = PRESETS[preset_name]()
    validate(spec)

    # Config roundtrip: serialise -> parse.
    cfg = spec_to_dict(spec, SimParams(n_steps=200, rng_seed=7), layout)
    cfg_round = json.loads(json.dumps(cfg, sort_keys=True))
    assert cfg == cfg_round, "JSON roundtrip failed"

    params = SimParams(n_steps=200, rng_seed=7)
    nc = run_simulation(spec, params, layout, output_dir=tmp / preset_name, binary=BINARY)
    result = load_network_netcdf(nc)

    n_steps = result.occupancy.shape[0]
    assert n_steps == params.n_steps, f"n_steps mismatch: {n_steps} vs {params.n_steps}"
    assert result.occupancy.shape[1] == len(result.lane_road_id), \
        "lane_index dim mismatches lane_road_id length"
    assert result.occupancy.shape[2] >= max(result.lane_length), \
        "cell dim smaller than max lane length"

    # final_occupancy mirrors occupancy[-1].
    assert (result.final_occupancy == result.occupancy[-1]).all(), \
        "final_occupancy mismatch"

    # config_json roundtrips through the global attribute.
    cfg_from_nc = result.config
    assert cfg_from_nc["schema_version"] == 1
    assert len(cfg_from_nc["roads"]) == len(spec.roads)
    assert len(cfg_from_nc["junctions"]) == len(spec.junctions)
    assert cfg_from_nc["params"]["model"] == "NS", \
        f"default model should roundtrip as 'NS', got {cfg_from_nc['params']['model']!r}"

    # Per-step mass conservation: N(t+1) - N(t) == entries(t+1) - exits(t+1).
    n_total = result.occupancy.reshape(n_steps, -1).sum(axis=1)
    delta = n_total[1:] - n_total[:-1]
    entries_sum = result.road_entries[1:].sum(axis=1)
    exits_sum = result.road_exits[1:].sum(axis=1)
    diff = delta - (entries_sum - exits_sum)
    bad = (diff != 0).sum()
    assert bad == 0, f"{preset_name}: {bad} steps violate mass conservation; first diffs {diff[diff != 0][:5]}"

    print(f"  PASS: shape={result.occupancy.shape}  mean_density={result.road_density.mean():.3f}")


def assert_tasep_rejected(tmp: Path):
    """The Fortran driver should refuse to run with model='TASEP' (not implemented)."""
    print("--- TASEP rejection ---")
    spec, layout = PRESETS["crossroads"]()
    params = SimParams(n_steps=50, rng_seed=1, model="TASEP")
    try:
        run_simulation(spec, params, layout, output_dir=tmp / "tasep_reject", binary=BINARY)
    except RuntimeError as exc:
        print(f"  PASS: {exc}")
        return
    raise AssertionError("expected RuntimeError when model='TASEP'")


def main():
    if not BINARY.exists():
        sys.exit(f"missing binary {BINARY} — run `make run_network` first")

    tmp = Path(tempfile.mkdtemp(prefix="px915_pipeline_"))
    try:
        for name in PRESETS:
            assert_pipeline(name, tmp)
        assert_tasep_rejected(tmp)
        print("ALL OK")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
