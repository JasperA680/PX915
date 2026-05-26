"""Python interface to the Fortran fundamental-diagram sweep driver.

The Fortran ``build/fd_sweep`` binary runs both TASEP (open-boundary
chain, α-then-β sweep with the unswept boundary held at 1 so the bulk
sees a deterministic reservoir and the J(ρ) curve recovers the
textbook ``min(ρ, 1−ρ)`` shape) and NS (periodic ring, density sweep)
sweeps; ``run_fd_sweep`` invokes it and ``load_fd_netcdf`` reads the
result.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Union

import numpy as np
import netCDF4 as nc


_DEFAULT_EXE = Path(__file__).parents[2] / "build" / "fd_sweep"
_DEFAULT_OUT = Path(__file__).parents[2] / "data" / "output" / "fundamental_diagram.nc"

# Per-model default measurement length.  These match the values the
# previous pure-Python implementation in ``analysis.py`` used so that
# saved sweeps come out at the same noise level as before.
_DEFAULT_N_STEPS = {"TASEP": 3000, "NS": 1500}


def run_fd_sweep(model: str, L: int, n_points: int,
                 v_max: int = 5, p_slow: float = 0.2,
                 n_steps: int = None, seed: int = None,
                 output_path: Union[str, Path] = _DEFAULT_OUT,
                 exe: Union[str, Path] = _DEFAULT_EXE) -> Path:
    """Invoke the Fortran FD sweep binary.

    Parameters
    ----------
    model:
        ``"TASEP"`` (open-boundary 1D chain, α/β sweep) or ``"NS"``
        (periodic ring, density sweep).
    L:
        Lattice length (sites).
    n_points:
        Number of values per swept branch.  TASEP produces ``2 * n_points``
        output rows (α-sweep then β-sweep); NS produces ``n_points``.
    v_max, p_slow:
        Only used by the NS sweep; TASEP ignores them.
    n_steps:
        Measurement window in time steps.  Defaults to 3000 for TASEP and
        1500 for NS — matches the pure-Python defaults the sweep replaces.
    seed:
        Optional integer to reseed the Fortran RNG so the sweep is
        reproducible across runs.
    output_path:
        NetCDF file to write.
    exe:
        Path to the compiled ``fd_sweep`` binary.
    """
    if model not in ("TASEP", "NS"):
        raise ValueError(f"unknown FD sweep model: {model!r}")
    if n_steps is None:
        n_steps = _DEFAULT_N_STEPS[model]

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(exe), model,
        str(int(L)), str(int(n_points)), str(int(n_steps)),
        str(int(v_max)), str(float(p_slow)),
        str(output_path),
    ]
    if seed is not None:
        cmd.append(str(int(seed)))

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"fd_sweep failed (exit {result.returncode}):\n{result.stderr}"
        )
    return output_path


def load_fd_netcdf(path: Union[str, Path]) -> dict:
    """Read a fundamental-diagram NetCDF written by ``fd_sweep``.

    Returns
    -------
    dict with keys ``rho``, ``J`` (each a ``numpy.ndarray``) plus
    ``model``, ``L``, ``n_burnin``, ``n_measure``, ``v_max``, ``p_slow``
    from the file's global attributes.
    """
    with nc.Dataset(path, "r") as ds:
        rho = np.array(ds.variables["rho"][:], dtype=float)
        J   = np.array(ds.variables["J"][:],   dtype=float)
        attrs = {k: ds.getncattr(k) for k in ds.ncattrs()}
    return dict(rho=rho, J=J, **attrs)
