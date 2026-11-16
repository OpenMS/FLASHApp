"""Draw a small table as an openms-insight ``SeriesPlot`` in the dataset's cache:
the template's ``src/viz/series.py`` with FLASHApp's cache folder and the
explanation of the component being drawn.
"""

from pathlib import Path

from src.viz import series as library
from src.viz.series import BOX_POINTS, thin  # noqa: F401  pylint: disable=unused-import  (tests use them)


def show_series(data, name: str, cache_path, key: str, height=None, state_manager=None, help=None, **plot_kwargs):  # pylint: disable=redefined-builtin
    """Draw ``data`` (polars or pandas) as a ``SeriesPlot``; ``plot_kwargs`` go to it."""
    if help is None:
        from src.render import teaching
        help = teaching.help()
    return library.show_series(data, name, cache_path=Path(cache_path, 'insight'), key=key,
                               state_manager=state_manager, height=height, help=help, **plot_kwargs)
