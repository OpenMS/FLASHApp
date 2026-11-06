"""Helpers the Layout Manager dialogs share.

The Viewer arranges its components with the view bar (``src/layout_manager``);
the comparison grid still opens this package's editor in a dialog
(``src/render/quant_view.py``), which uses the help text and the FileManager
helper below.
"""

from src.workflow.FileManager import FileManager

_HELP = (
    "Drag components into the grid (up to 3 per row), add experiments to "
    "compare several files, then press **Save**."
)


def _own_file_manager(file_manager):
    """A FileManager on the same cache with its own SQLite connection.

    A dialog reruns as a fragment, possibly on another script thread, and an
    SQLite connection may only be used by the thread that opened it.
    """
    return FileManager(file_manager.workflow_dir, file_manager.cache_path)
