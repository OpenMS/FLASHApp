"""Viewer cache folders stay short enough for Windows' path length limit."""

from pathlib import Path

from src.render.cache_dirs import KEEP, insight_cache_dir, short_name

LONG_RUN = ('FLASHIda_methodLCMStqs0p9MS3_20261103_Caco2Histone_RealTimeIDMS3_60min_2ul_MS2_'
            '20261109161635_2026-11-09-16-16-36')


def test_short_names_are_kept():
    assert short_name('run_20261109-161635') == 'run_20261109-161635'
    assert short_name('x' * (KEEP + 9)) == 'x' * (KEEP + 9)


def test_long_names_are_cut_and_hashed():
    name = short_name(LONG_RUN)
    assert len(name) <= KEEP + 9
    assert name.startswith(LONG_RUN[:KEEP])
    # Names that differ only after the kept part get folders of their own
    assert short_name(LONG_RUN + 'b') != name
    assert short_name(LONG_RUN) == name


def test_no_trailing_dot_or_space():
    name = short_name('a' * (KEEP - 1) + '.' + 'b' * 40)
    assert '.-' not in name


def test_windows_cache_path_fits():
    root = r'C:\Users\Administrator\AppData\Roaming\FLASHApp\workspaces\default\flashida\cache'
    folder = Path(root, 'insight-live', short_name(LONG_RUN),
                  'decision_map-883152f55527-fcb39173', 'preprocessed', 'level_0_categories.parquet')
    assert len(str(folder)) < 260
    assert insight_cache_dir('cache', LONG_RUN) == Path('cache', 'insight', short_name(LONG_RUN))


def test_one_build_lock_per_cache():
    from src.render.cache_dirs import build_lock
    assert build_lock('cache', 'a-1') is build_lock(Path('cache'), 'a-1')
    assert build_lock('cache', 'a-1') is not build_lock('cache', 'a-2')
