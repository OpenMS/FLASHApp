"""The FLASHIda live view (src/render/ida_live.py) refreshes while a run is
acquired: its elements stay drawn as they are during a refresh instead of
greying out, and a run's acquisition decisions are read again only once its
logs changed."""

from pathlib import Path

import pytest

pytest.importorskip('openms_insight')

from src.render import ida_live  # noqa: E402


def test_live_container_keeps_elements_drawn(monkeypatch):
    shown = []
    monkeypatch.setattr(ida_live.st, 'html', shown.append)
    monkeypatch.setattr(ida_live.st, 'container', lambda key: key)
    assert ida_live.live_container('live-view') == 'live-view'
    css = shown[0]
    assert '.st-key-live-view [data-stale="true"]' in css
    assert 'opacity:1 !important' in css


def test_decisions_are_read_again_only_when_the_logs_changed(monkeypatch, tmp_path):
    reads = []
    monkeypatch.setattr(ida_live.ida_log, 'read_run', lambda run: reads.append(run) or run)
    monkeypatch.setattr(ida_live.ida_log, 'decisions', lambda logs: {'run': str(logs), 'n': len(reads)})
    monkeypatch.setattr(ida_live, '_DECISIONS', {})
    run = Path(tmp_path, 'R1')
    first = ida_live.run_decisions(run, 'stamp-1')
    assert ida_live.run_decisions(run, 'stamp-1') is first
    assert len(reads) == 1
    assert ida_live.run_decisions(run, 'stamp-2')['n'] == 2
    assert len(reads) == 2


def test_decisions_keep_few_runs(monkeypatch, tmp_path):
    monkeypatch.setattr(ida_live.ida_log, 'read_run', lambda run: run)
    monkeypatch.setattr(ida_live.ida_log, 'decisions', lambda logs: {})
    monkeypatch.setattr(ida_live, '_DECISIONS', {})
    for i in range(6):
        ida_live.run_decisions(Path(tmp_path, f'R{i}'), 's')
    assert len(ida_live._DECISIONS) == 4
    assert str(Path(tmp_path, 'R5')) in ida_live._DECISIONS
