import pickle
import hashlib

def hash_complex(d):
    serialized = pickle.dumps(d)
    return hashlib.sha256(serialized).hexdigest()


def payload_hash(data, tracker_id):
    """Hash of a grid cell's payload as sent to the Vue component.

    The frontend skips re-parsing a render whose hash equals the previous one. The
    hash therefore has to change whenever the cell must rebuild, and that includes a
    new experiment: two runs of the same input produce byte-identical tables, so the
    data alone is not enough. The StateTracker id is new for every experiment and
    constant across reruns within one, so it keeps the skip for unchanged data and
    forces a rebuild (and the default row selection) when the experiment changes.
    """
    return hash_complex((data, tracker_id))
