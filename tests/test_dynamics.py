"""Dynamics: the correlation maths, and the MD command EVODYN would issue."""
import numpy as np
import pandas as pd
import pytest

from evodyn.dynamics import dccm_from_covariance
from evodyn.dynamics.equilibrator import build_command, run_equilibrator


def test_covariance_to_correlation(tmp_path):
    rng = np.random.default_rng(0)
    n = 4
    a = rng.normal(size=(3 * n, 3 * n))
    cov = a @ a.T
    p = tmp_path / "covar.ascii"
    np.savetxt(p, cov)
    dccm = dccm_from_covariance(p)
    assert dccm.shape == (n, n)
    assert np.allclose(np.diag(dccm), 1.0)
    assert np.allclose(dccm.to_numpy(), dccm.to_numpy().T)
    assert dccm.to_numpy().min() >= -1.0001 and dccm.to_numpy().max() <= 1.0001


def test_equilibrator_command_matches_the_published_protocol():
    cmd = build_command("prot.pdb", force_field="amber99sb-ildn", water_model="tip3p",
                        box_type="cubic", box_distance=1.2)
    assert cmd[:3] == ["EquilibraTor", "-p", "prot.pdb"]
    for flag, value in (("-gff", "amber99sb-ildn"), ("-gwm", "tip3p"),
                        ("-gbt", "cubic"), ("-gd", "1.2")):
        assert cmd[cmd.index(flag) + 1] == value


def test_replicates_get_separate_working_directories(tmp_path):
    runs = [run_equilibrator("prot.pdb", tmp_path, replicate=r, dry_run=True)
            for r in (1, 2)]
    assert runs[0].workdir != runs[1].workdir
    assert all(r.workdir.exists() for r in runs)
    # EquilibraTor names its project directory after the input structure, so
    # replicates sharing a working directory would overwrite each other
    assert runs[0].project_dir.name == runs[1].project_dir.name
