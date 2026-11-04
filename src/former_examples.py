"""Example runs of earlier FLASHApp releases in an existing workspace.

The local app copies the demo workspace only when the workspace folder is
missing, so a workspace from an earlier release keeps that release's example
runs. The upload counts every file that is not one of today's examples as the
user's own and then leaves the examples out, so such a workspace never shows
the current examples. ``retire`` takes the earlier copies out once, matched by
content (size, then sha256), never by name alone; the upload then puts today's
examples in their place unless the user has runs of their own.
"""

import hashlib
from pathlib import Path

# (size in bytes, sha256) of every example MS run an earlier release shipped
# in example-data/workspaces/default, from the repository history. Runs that
# were only renamed since are recognised by their content as today's examples.
FORMER_EXAMPLES = {
    (41747, '26ed8150399a2e25a51689088cd3c648a3df7e6a3a6ec238fb2ea0c38ce59da2'),
    (43872, '69f8b735174e6e4f8525caf2f0d8d1b804ea43f27b45c69675d0ef53305de2ad'),
    (39657, 'fbd6aa79d854ccd3b1a3f56dae99b494350c7ca03705c32263ac73fb36e78181'),
    (138626, '58a4f209198951117a2bf7f13a50e503a5276d6b1978071190582d79ecf409f5'),
    (875547, 'bce1715c2ea16ca106c52b4b57890701e44527a493f3d7ffbf6f0b321f333554'),
    (177838, '8ecc4154e2fb73e293a13fd98b25a7d80c7af33b814a79c83b472182b66b2afc'),
    (44670, '5cf589c82df3c6a0b8e56ada8ff3214bd325c07771f5819dab7772fc6b83a2e4'),
    (1357278, 'ac66399db2cea37e9a8653fd3b2d54ef2c37b6d74217bf25a2dce0c4b573c984'),
    (251813, '548aad45408972f6e887382c8eb9592a4d4630366155b90400d4adddf8fc63e6'),
    (11975352, '9d72dbd5140b613a3036d202571f8abdaa84f120797dbcf63c115d7b0f3990b2'),
    (12246895, '3f673807ec6a4ecdaa89992d6266a9521e2746c6aa817f86c471cc46f2bf2f2c'),
    (15234292, 'a0142c4dcdff7e5ee41bf28ce7a3d8e8e24173aa88c5944ad8b6806a8a0363b2'),
    (16017341, 'aecaddf712ae230328e612a59ab534a839eb3da2c57a9d618778581542fd0b18'),
    (7593480, 'cb42c0f05e95fa5d29cfb305e6ed90e322c141adc117cfa875a4fe9f9b5e58ec'),
    (7615270, '6392d27d8549a69e5940719f855f3387def1c67bf61f3fb9f9a0a76d9ef40aa9'),
    (9551199, '7663d38565aa3022013058550d10d38960e0705eec443250439abcf3641c8c9c'),
    (9873880, '22394da3a25d039901f2b01446bcedfaa9b544b84c7388698498e68abbf4a602'),
}

MARKER = '.former-examples-retired'
EXTERNAL = 'external_files.txt'


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


class _Matcher:
    """Tells whether a file has the content of an example, past or present.
    Only a file whose size matches one is hashed."""

    def __init__(self, examples):
        self.known = set(FORMER_EXAMPLES)
        self.current = [Path(e) for e in examples if Path(e).is_file()]
        self.sizes = {s for s, _ in self.known} | {e.stat().st_size for e in self.current}
        self._current_hashed = False

    def __call__(self, path: Path) -> bool:
        try:
            size = path.stat().st_size
        except OSError:
            return False
        if size not in self.sizes:
            return False
        if not self._current_hashed:
            self.known |= {(e.stat().st_size, _sha256(e)) for e in self.current}
            self._current_hashed = True
        return (size, _sha256(path)) in self.known


def is_example(path, examples) -> bool:
    """Whether ``path`` holds the content of an example run (any release)."""
    return _Matcher(examples)(Path(path))


def retire(folder, examples, marker) -> list:
    """Remove the copies of earlier releases' example runs from the MS data
    folder ``folder`` and its references in external_files.txt, once: the
    file ``marker`` (outside ``folder``, which lists every file as an input)
    records that it ran. A file named like one of today's ``examples`` is left
    to the upload. Returns the names removed."""
    folder = Path(folder)
    marker = Path(marker)
    if marker.exists() or not folder.is_dir():
        return []
    current = {Path(e).name for e in examples}
    matches = _Matcher(examples)
    removed = []
    for p in folder.iterdir():
        if p.is_file() and p.name != EXTERNAL and p.name not in current and matches(p):
            p.unlink()
            removed.append(p.name)
    external = Path(folder, EXTERNAL)
    if external.exists():
        lines = [line for line in external.read_text().splitlines() if line.strip()]
        keep = [line for line in lines if not matches(Path(line))]
        removed += [Path(line).name for line in lines if line not in keep]
        if len(keep) < len(lines):
            if keep:
                external.write_text(''.join(f'{line}\n' for line in keep))
            else:
                external.unlink()
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text('\n'.join(removed))
    return removed
