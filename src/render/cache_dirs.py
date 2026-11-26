"""Short folder names for viewer caches.

A cache folder named after a run or dataset (a raw file name plus a time
stamp) can push the cache files past Windows' 260-character path limit:
openms-insight writes them three folders deeper. `short_name` keeps the
start of the name, so the folder is still recognisable, and replaces the
rest with a hash, so different names never share a folder.
Streamlit-free.
"""

import hashlib
import threading
from pathlib import Path

# Characters of the name kept before the hash
KEEP = 32


def short_name(name: str, keep: int = KEEP) -> str:
    """`name` itself when it is short, else its first `keep` characters and
    eight characters of its hash."""
    name = str(name)
    if len(name) <= keep + 9:
        return name
    digest = hashlib.md5(name.encode('utf-8')).hexdigest()[:8]
    # Windows drops a trailing dot or space from a folder name
    return f"{name[:keep].rstrip('. ')}-{digest}"


def insight_cache_dir(cache_path, dataset_id):
    """Folder of a dataset's viewer component caches below a workspace cache
    (`file_manager.cache_path`)."""
    return Path(cache_path, 'insight', short_name(dataset_id))


def _viewer_caches(cache_path, dataset_id):
    """The viewer caches removed with a dataset's results, or all of them
    (`dataset_id` None) when the whole store is cleared."""
    if dataset_id is None:
        return Path(cache_path, 'insight')
    return insight_cache_dir(cache_path, dataset_id)


def register_cache_hook():
    """Has FileManager.remove_results / clear_cache delete the viewer caches."""
    from src.workflow.FileManager import FileManager
    if _viewer_caches not in FileManager.cache_dir_hooks:
        FileManager.cache_dir_hooks.append(_viewer_caches)


register_cache_hook()


_BUILD_LOCKS: dict[str, threading.Lock] = {}
_BUILD_LOCKS_GUARD = threading.Lock()


def build_lock(cache_dir, cache_id) -> threading.Lock:
    """Lock held while one component cache is built.

    Two browser tabs that see the same new results would otherwise build the
    same cache at once, and on Windows one of them cannot replace a file the
    other has open ("Could not replace cache file"). The second tab waits and
    then loads what the first one built. One lock per cache, within this
    server process (the desktop app runs one).
    """
    key = str(Path(cache_dir, cache_id))
    with _BUILD_LOCKS_GUARD:
        return _BUILD_LOCKS.setdefault(key, threading.Lock())
