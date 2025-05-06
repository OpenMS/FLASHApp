import pickle
import hashlib

def hash_complex(d):
    serialized = pickle.dumps(d)
    return hashlib.sha256(serialized).hexdigest()