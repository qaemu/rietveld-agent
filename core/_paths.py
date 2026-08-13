"""Repository-root resolution for the run-from-clone layout.

The engine is designed to run from a clone of the repository: data,
governance schemas and policies live in the repo tree, not inside the
installed package.  ``repo_root()`` finds that tree from the current
working directory first (covers installed console scripts started inside
a clone), then from the package location.
"""

import os
from pathlib import Path


def repo_root() -> Path:
    """Return the repository root (dir containing ``data/`` and ``README.md``)."""
    here = Path(__file__).resolve().parent
    for cand in [Path(os.getcwd()).resolve(), here, *here.parents]:
        if (cand / "data").is_dir() and (cand / "README.md").is_file():
            return cand
    return here.parent