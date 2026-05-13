"""Shared pytest fixtures for odp.

These paths replace the ad-hoc sys.path.insert dance in test files. Any
test that needs to import a module from `source/`, `scripts/`, or
`dependencies/fasta-parser/` should depend on the matching fixture (which
also prepends to sys.path the first time it is resolved).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return Path(__file__).resolve().parent


@pytest.fixture(scope="session")
def source_dir(repo_root: Path) -> Path:
    path = repo_root / "source"
    sys.path.insert(0, str(path))
    return path


@pytest.fixture(scope="session")
def scripts_dir(repo_root: Path) -> Path:
    path = repo_root / "scripts"
    sys.path.insert(0, str(path))
    return path


@pytest.fixture(scope="session")
def fasta_parser_dir(repo_root: Path) -> Path:
    path = repo_root / "dependencies" / "fasta-parser"
    sys.path.insert(0, str(path))
    return path


@pytest.fixture(scope="session")
def tests_dir(repo_root: Path) -> Path:
    return repo_root / "tests"


@pytest.fixture(scope="session")
def testdb_dir(tests_dir: Path) -> Path:
    return tests_dir / "testdb"


@pytest.fixture(scope="session")
def mini_hydra(testdb_dir: Path) -> Path:
    path = testdb_dir / "mini_hydra"
    if not path.is_dir():
        pytest.skip(f"mini_hydra fixture missing at {path}")
    return path


@pytest.fixture(scope="session")
def mini_urchin(testdb_dir: Path) -> Path:
    path = testdb_dir / "mini_urchin"
    if not path.is_dir():
        pytest.skip(f"mini_urchin fixture missing at {path}")
    return path


@pytest.fixture(scope="session")
def mini_lg_db(testdb_dir: Path) -> Path:
    path = testdb_dir / "mini_LG_db"
    if not path.is_dir():
        pytest.skip(f"mini_LG_db fixture missing at {path}")
    return path
