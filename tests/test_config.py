import pathlib

import pytest
import yaml

from evodyn.config import load_config
from evodyn.exceptions import ConfigError


def _write(tmp_path, data):
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump(data))
    return p


def test_minimal_config_loads(tmp_path):
    cfg = load_config(_write(tmp_path, {
        "target": {"name": "X", "sequence": "ACDEF", "reference_start": 1},
        "sequences": {"source": "none"},
        "coevolution": {"method": "none", "couplings": "c.csv",
                        "couplings_format": "translated"},
        "dynamics": {"mode": "none"},
    }))
    assert cfg.target.name == "X"
    assert cfg.dynamics.mode == "none"


def test_equilibrate_without_a_structure_is_rejected(tmp_path):
    with pytest.raises(ConfigError, match="starting structure"):
        load_config(_write(tmp_path, {
            "target": {"name": "X"},
            "dynamics": {"mode": "equilibrate"},
        }))


def test_unknown_key_is_rejected(tmp_path):
    with pytest.raises(ConfigError, match="unknown key"):
        load_config(_write(tmp_path, {
            "target": {"name": "X"},
            "network": {"threshhold": "percolation"},
        }))


def test_input_paths_resolve_relative_to_the_config(tmp_path):
    (tmp_path / "sub").mkdir()
    cfg = load_config(_write(tmp_path, {
        "target": {"name": "X", "structure": "sub/s.pdb"},
        "dynamics": {"mode": "none"},
    }))
    assert pathlib.Path(cfg.target.structure).is_absolute()
    assert pathlib.Path(cfg.target.structure).parent.name == "sub"
