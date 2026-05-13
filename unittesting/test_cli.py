"""Tests for the bin/odp CLI wrapper.

Exercises the argument-parsing layer and the command-construction logic
without actually invoking snakemake.
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
CLI_PATH = REPO_ROOT / "bin" / "odp"


# Load bin/odp as a regular module so we can unit-test internal helpers.
# importlib won't auto-detect Python source because the file has no .py
# extension; we have to hand it an explicit SourceFileLoader.
def _load_cli_module():
    from importlib.machinery import SourceFileLoader
    loader = SourceFileLoader("odp_cli", str(CLI_PATH))
    spec = importlib.util.spec_from_loader("odp_cli", loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def cli():
    return _load_cli_module()


# ---- Module-level helpers ----------------------------------------------


def test_repo_root_resolves(cli):
    assert cli.repo_root() == REPO_ROOT


def test_snakefiles_map_all_resolve(cli):
    """Every Snakefile name the CLI dispatches to must actually exist on disk."""
    for name in cli.SNAKEFILES:
        sf = cli.snakefile_for(name)
        assert sf.is_file(), f"missing Snakefile for subcommand {name!r}: {sf}"


def test_default_cores_returns_positive_int(cli):
    n = cli.default_cores()
    assert isinstance(n, int)
    assert n >= 1


# ---- argparse plumbing -------------------------------------------------


def test_parser_lists_every_subcommand(cli):
    parser = cli.build_parser()
    # Hook to inspect subparsers without invoking them
    subparsers_action = next(
        a for a in parser._actions if a.dest == "subcommand"
    )
    names = set(subparsers_action.choices.keys())
    expected = set(cli.SNAKEFILES) | {"init", "validate", "version"}
    assert names == expected


def test_no_subcommand_exits_nonzero(cli):
    """`odp` with no subcommand should be an error (subparser is required)."""
    parser = cli.build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args([])
    assert exc.value.code != 0


def test_subcommand_help_executes(cli, capsys):
    """`odp run --help` should print help and exit 0."""
    parser = cli.build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["run", "--help"])
    assert exc.value.code == 0
    out = capsys.readouterr().out
    assert "snakemake" in out.lower() or "cores" in out.lower()


@pytest.mark.parametrize("subcommand", [
    "run", "only-db", "filecheck", "nway-rbh",
    "rbh-to-alignments", "rbh-to-hmm", "rbh-to-ribbon", "plot-mixing",
])
def test_each_run_subcommand_parses(cli, subcommand):
    """All snakefile-backed subcommands should accept the common arg shape."""
    parser = cli.build_parser()
    args = parser.parse_args([
        subcommand,
        "--cores", "4",
        "--config", "config.yaml",
        "--dry-run",
    ])
    assert args.subcommand == subcommand
    assert args.cores == 4
    assert args.config == "config.yaml"
    assert args.dry_run is True


def test_rerun_incomplete_default_on(cli):
    parser = cli.build_parser()
    args = parser.parse_args(["run"])
    assert args.rerun_incomplete is True


def test_no_rerun_incomplete_flag(cli):
    parser = cli.build_parser()
    args = parser.parse_args(["run", "--no-rerun-incomplete"])
    assert args.rerun_incomplete is False


def test_snakemake_arg_passthrough_repeatable(cli):
    # argparse needs `--flag=--value` form when the value starts with '--'
    # so it doesn't confuse it for another flag.
    parser = cli.build_parser()
    args = parser.parse_args([
        "run",
        "--snakemake-arg=--quiet",
        "--snakemake-arg=--keep-going",
    ])
    assert args.snakemake_arg == ["--quiet", "--keep-going"]


# ---- build_snakemake_cmd -----------------------------------------------


def test_build_cmd_minimal(cli):
    parser = cli.build_parser()
    args = parser.parse_args(["run"])
    cmd = cli.build_snakemake_cmd(
        snakefile=cli.snakefile_for("run"),
        cores=args.cores,
        workdir=None,
        configfile=None,
        dry_run=False,
        rerun_incomplete=False,
        printshellcmds=False,
        extra=[],
    )
    assert cmd[0].endswith("snakemake")
    assert "--snakefile" in cmd
    assert str(cli.snakefile_for("run")) in cmd
    assert "--cores" in cmd
    assert str(args.cores) in cmd


def test_build_cmd_all_flags(cli, tmp_path):
    cfg = tmp_path / "config.yaml"
    cfg.write_text("species: {}\n")
    cmd = cli.build_snakemake_cmd(
        snakefile=cli.snakefile_for("only-db"),
        cores=12,
        workdir=tmp_path,
        configfile=cfg,
        dry_run=True,
        rerun_incomplete=True,
        printshellcmds=True,
        extra=["--quiet"],
    )
    assert "--directory" in cmd
    assert str(tmp_path.resolve()) in cmd
    assert "--configfile" in cmd
    assert str(cfg.resolve()) in cmd
    assert "--dry-run" in cmd
    assert "--rerun-incomplete" in cmd
    assert "--printshellcmds" in cmd
    assert "--quiet" in cmd
    assert "12" in cmd  # cores value


# ---- init subcommand ---------------------------------------------------


def test_init_writes_starter_config(cli, tmp_path, monkeypatch, capsys):
    out = tmp_path / "config.yaml"
    parser = cli.build_parser()
    args = parser.parse_args(["init", "--output", str(out)])
    rc = cli.cmd_init(args)
    assert rc == 0
    assert out.is_file()
    body = out.read_text()
    assert "species:" in body
    assert "diamond_or_blastp" in body


def test_init_refuses_overwrite_without_force(cli, tmp_path):
    out = tmp_path / "config.yaml"
    out.write_text("existing content")
    parser = cli.build_parser()
    args = parser.parse_args(["init", "--output", str(out)])
    rc = cli.cmd_init(args)
    assert rc == 1
    assert out.read_text() == "existing content"


def test_init_overwrites_with_force(cli, tmp_path):
    out = tmp_path / "config.yaml"
    out.write_text("existing content")
    parser = cli.build_parser()
    args = parser.parse_args(["init", "--output", str(out), "--force"])
    rc = cli.cmd_init(args)
    assert rc == 0
    assert "existing content" not in out.read_text()
    assert "species:" in out.read_text()


# ---- subprocess-level smoke check (no snakemake invoked) ---------------


def test_cli_top_level_help_via_subprocess(cli):
    """Invoking the CLI as a subprocess should print help cleanly."""
    result = subprocess.run(
        [sys.executable, str(CLI_PATH), "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "SUBCOMMAND" in result.stdout
    assert "rbh-to-ribbon" in result.stdout


def test_cli_version_via_subprocess(cli):
    result = subprocess.run(
        [sys.executable, str(CLI_PATH), "version"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "odp checkout" in result.stdout
