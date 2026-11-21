"""Validate the manifest and compare both certificate runs with the release receipt."""

import hashlib
from pathlib import Path
import subprocess
import sys


def main():
    root = Path(__file__).resolve().parent
    for line in (root / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split("  ", 1)
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != digest:
            raise SystemExit(f"Checksum mismatch in {name}")
    output = []
    for script in ("verify.py", "independent.py"):
        run = subprocess.run([sys.executable, "-B", script], cwd=root,
                             check=True, stdout=subprocess.PIPE)
        output.append(run.stdout)
        sys.stdout.buffer.write(run.stdout)
        sys.stdout.buffer.flush()
    if b"".join(output) != (root / "evidence/check.txt").read_bytes():
        raise SystemExit("Certificate output differs from the release receipt")


if __name__ == "__main__":
    main()
