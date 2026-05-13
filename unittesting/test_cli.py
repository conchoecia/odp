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
    expected = set(cli.SNAKEFILES) | {
        "init", "validate", "version", "orthologs-to-rbh", "orthofinder-import",
    }
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
        retries=0,
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
        retries=0,
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
    # `odp: ...` line carries the project version (git-describe output or
    # "unknown" if .git is absent / git is not on PATH).
    assert "odp:" in result.stdout


def test_odp_version_helper_returns_string(cli):
    v = cli.odp_version()
    assert isinstance(v, str)
    assert v != ""


# ---- UX improvements: --version, typo-suggest, pre-flight ---------------


def test_top_level_version_flag(cli):
    """`odp --version` should exit 0 and print 'odp <something>'."""
    result = subprocess.run(
        [sys.executable, str(CLI_PATH), "--version"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert result.stdout.startswith("odp ")


def test_unknown_subcommand_suggests_correction(cli):
    result = subprocess.run(
        [sys.executable, str(CLI_PATH), "ru"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "unknown subcommand 'ru'" in result.stderr
    assert "did you mean 'run'" in result.stderr


def test_unknown_subcommand_far_off_no_guess(cli):
    """Garbage subcommand far from any valid name: no false suggestion."""
    result = subprocess.run(
        [sys.executable, str(CLI_PATH), "xyzqwerty"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "unknown subcommand 'xyzqwerty'" in result.stderr
    assert "did you mean" not in result.stderr


def test_run_with_missing_config_dies_cleanly(cli):
    """No snakemake traceback for a missing config — one-line CLI error."""
    result = subprocess.run(
        [sys.executable, str(CLI_PATH), "run", "--config", "/nonexistent/x.yaml"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "config file not found" in result.stderr
    # No traceback lines from snakemake.
    assert "Traceback" not in result.stderr


def test_run_with_missing_workdir_dies_cleanly(cli):
    result = subprocess.run(
        [sys.executable, str(CLI_PATH), "run", "--workdir", "/nonexistent/dir"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "workdir not found" in result.stderr


def test_subcommand_order_in_help(cli):
    """Help text should list pipeline subcommands first, utilities last."""
    parser = cli.build_parser()
    help_text = parser.format_help()
    # `run` should appear before `init`; `version` should be last.
    assert help_text.index("\n    run") < help_text.index("\n    init")
    assert help_text.index("\n    init") < help_text.index("\n    version")


def test_validate_auto_detects_local_config(cli, tmp_path, monkeypatch):
    """`odp validate` with no --config should pick up ./config.yaml."""
    cfg = tmp_path / "config.yaml"
    cfg.write_text("species: {}\n")
    monkeypatch.chdir(tmp_path)
    parser = cli.build_parser()
    args = parser.parse_args(["validate"])
    # We don't want to actually exec snakemake in a unit test, so monkeypatch
    # subprocess.call inside the cli module to no-op + remember the cmd.
    seen = {}
    def fake_call(cmd):
        seen["cmd"] = cmd
        return 0
    monkeypatch.setattr(cli.subprocess, "call", fake_call)
    rc = cli.cmd_validate(args)
    assert rc == 0
    # The CLI should have resolved ./config.yaml as the config.
    assert any("config.yaml" in part for part in seen["cmd"])


def test_validate_no_config_no_local_fails(cli, tmp_path, monkeypatch):
    """No --config and no ./config.yaml → error, not snakemake invoked."""
    monkeypatch.chdir(tmp_path)
    parser = cli.build_parser()
    args = parser.parse_args(["validate"])
    rc = cli.cmd_validate(args)
    assert rc == 2


# ---- Cores over-provisioning warning ------------------------------------


def test_warn_if_over_provisioned_emits_when_too_high(cli, capsys, monkeypatch):
    monkeypatch.setattr(cli, "available_cores", lambda: 8)
    cli.warn_if_over_provisioned(9999)
    captured = capsys.readouterr()
    assert "exceeds detected logical cores" in captured.err
    assert "8" in captured.err


def test_warn_if_over_provisioned_silent_when_under(cli, capsys, monkeypatch):
    monkeypatch.setattr(cli, "available_cores", lambda: 8)
    cli.warn_if_over_provisioned(4)
    captured = capsys.readouterr()
    assert captured.err == ""


def test_warn_if_over_provisioned_silent_at_exact_limit(cli, capsys, monkeypatch):
    monkeypatch.setattr(cli, "available_cores", lambda: 8)
    cli.warn_if_over_provisioned(8)
    captured = capsys.readouterr()
    assert captured.err == ""


def test_run_with_wild_cores_still_proceeds_to_config_check(cli):
    """--cores 9999 should warn but still continue to the config check."""
    result = subprocess.run(
        [sys.executable, str(CLI_PATH),
         "run", "--cores", "9999", "--config", "/nonexistent.yaml"],
        capture_output=True,
        text=True,
    )
    # exit 2 from config not found, not from cores
    assert result.returncode == 2
    assert "exceeds detected logical cores" in result.stderr
    assert "config file not found" in result.stderr


# ---- argcomplete hook ---------------------------------------------------


def test_cli_has_argcomplete_marker(cli):
    """`# PYTHON_ARGCOMPLETE_OK` must be in the first 1024 bytes of bin/odp
    for argcomplete's `register-python-argcomplete` to opt the file in."""
    head = CLI_PATH.read_text()[:1024]
    assert "PYTHON_ARGCOMPLETE_OK" in head


# ---- Auto-workdir from config -------------------------------------------


def test_auto_workdir_from_config(cli, tmp_path, monkeypatch):
    """If --config is given and --workdir is omitted, the CLI should run
    snakemake with --directory set to the config's parent directory.
    Configs frequently contain relative paths, so this is the
    do-what-you-mean default."""
    cfg = tmp_path / "config.yaml"
    cfg.write_text("species: {}\n")

    parser = cli.build_parser()
    args = parser.parse_args(["run", "--config", str(cfg), "--dry-run"])

    seen = {}
    def fake_call(c):
        seen["cmd"] = c
        return 0
    monkeypatch.setattr(cli.subprocess, "call", fake_call)
    rc = cli.run_subcommand("run", args)
    assert rc == 0
    cmd = seen["cmd"]
    assert "--directory" in cmd
    idx = cmd.index("--directory")
    assert cmd[idx + 1] == str(tmp_path.resolve())


def test_explicit_workdir_wins_over_auto(cli, tmp_path, monkeypatch):
    """If --workdir is explicit, the auto-set must not clobber it."""
    cfg = tmp_path / "config.yaml"
    cfg.write_text("species: {}\n")
    other = tmp_path / "elsewhere"
    other.mkdir()

    parser = cli.build_parser()
    args = parser.parse_args([
        "run", "--config", str(cfg),
        "--workdir", str(other),
        "--dry-run",
    ])

    seen = {}
    def fake_call(c):
        seen["cmd"] = c
        return 0
    monkeypatch.setattr(cli.subprocess, "call", fake_call)
    cli.run_subcommand("run", args)
    idx = seen["cmd"].index("--directory")
    assert seen["cmd"][idx + 1] == str(other.resolve())


# ---- --retries flag -----------------------------------------------------


def test_retries_default_zero_no_config_kwarg(cli):
    """With --retries 0 (the default) no '--config retries=...' is added."""
    cmd = cli.build_snakemake_cmd(
        snakefile=cli.snakefile_for("run"),
        cores=2,
        workdir=None,
        configfile=None,
        dry_run=False,
        rerun_incomplete=False,
        printshellcmds=False,
        retries=0,
        extra=[],
    )
    joined = " ".join(cmd)
    assert "retries=" not in joined


def test_retries_nonzero_injects_config_kwarg(cli):
    cmd = cli.build_snakemake_cmd(
        snakefile=cli.snakefile_for("run"),
        cores=2,
        workdir=None,
        configfile=None,
        dry_run=False,
        rerun_incomplete=False,
        printshellcmds=False,
        retries=3,
        extra=[],
    )
    assert "--config" in cmd
    assert "retries=3" in cmd


def test_retries_flag_parsed(cli):
    parser = cli.build_parser()
    args = parser.parse_args(["run", "--retries", "5"])
    assert args.retries == 5


def test_retries_default_is_zero(cli):
    parser = cli.build_parser()
    args = parser.parse_args(["run"])
    assert args.retries == 0
