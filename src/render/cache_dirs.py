"""Short folder names for viewer caches.

A cache folder named after a run or dataset (a raw file name plus a time
stamp) can push the cache files past Windows' 260-character path limit:
openms-insight writes them three folders deeper. `short_name` keeps the
start of the name, so the folder is still recognisable, and replaces the
rest with a hash, so different names never share a folder.
Streamlit-free.
"""

import hashlib
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
