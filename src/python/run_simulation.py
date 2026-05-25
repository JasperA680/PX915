"""Launcher that runs the Fortran network driver from a Python NetworkSpec.

Writes a NetCDF config, invokes ./build/run_network, parses the 'step k/n'
lines on stdout to drive a progress callback, and returns the NC path.
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path
from typing import Callable, Optional

from python.road_network import NetworkSpec, SimParams, LayoutSpec, validate
from python.io import write_config_netcdf


_STEP_RE = re.compile(r"^step (\d+)/(\d+)\s*$")


def run_simulation(
    spec: NetworkSpec,
    params: SimParams,
    layout: LayoutSpec,
    output_dir: os.PathLike,
    binary: os.PathLike = "./build/run_network",
    progress: Optional[Callable[[float, str], None]] = None,
    config_name: str = "config.nc",
    output_name: str = "result.nc",
) -> Path:
    """Run the Fortran simulator and return the output NetCDF path.

    `progress(frac, message)` is called with `frac` in [0, 1] as the binary
    emits 'step k/n' lines, and with frac=1.0 + a 'wrote: ...' message at end.
    """
    validate(spec)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = output_dir / config_name
    nc_path  = output_dir / output_name

    write_config_netcdf(spec, params, layout, cfg_path)

    proc = subprocess.Popen(
        [str(binary), str(cfg_path), str(nc_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    # Keep a tail of non-progress lines so we can include the binary's error
    # message in the RuntimeError if it exits non-zero.
    tail: list[str] = []
    assert proc.stdout is not None
    for line in proc.stdout:
        line = line.rstrip()
        m = _STEP_RE.match(line)
        if m and progress is not None:
            k, n = int(m.group(1)), int(m.group(2))
            progress(k / n, line)
        else:
            tail.append(line)
            if len(tail) > 20:
                tail.pop(0)
            if progress is not None:
                progress(-1.0, line)

    rc = proc.wait()
    if rc != 0:
        msg = f"run_network exited with code {rc}"
        if tail:
            msg += ": " + " | ".join(tail[-3:])
        raise RuntimeError(msg)
    return nc_path


def write_config(
    spec: NetworkSpec,
    params: SimParams,
    layout: LayoutSpec,
    path: os.PathLike,
) -> Path:
    """Serialise a spec + params + layout to a NetCDF config file."""
    validate(spec)
    return write_config_netcdf(spec, params, layout, path)
