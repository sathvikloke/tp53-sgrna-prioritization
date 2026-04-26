"""End-to-end smoke test for the CLI argparser (does not run the network)."""

from __future__ import annotations

import pytest

from src.cli import build_parser


def test_cli_parser_help_exits_zero(capsys):
    parser = build_parser()
    with pytest.raises(SystemExit) as ei:
        parser.parse_args(["--help"])
    assert ei.value.code == 0


@pytest.mark.parametrize("cmd", [
    "fetch-cds", "fetch-cohorts", "enumerate", "score-on-target",
    "rank", "tables",
])
def test_cli_subcommand_known(cmd):
    parser = build_parser()
    args = parser.parse_args([cmd])
    assert args.cmd == cmd
    assert callable(args.func)


def test_cli_score_off_target_device_default():
    parser = build_parser()
    args = parser.parse_args(["score-off-target"])
    assert args.device == "C"


def test_cli_unknown_command_errors():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["does-not-exist"])
