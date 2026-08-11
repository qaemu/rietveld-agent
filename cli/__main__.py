"""Enable `python -m cli ...`."""
import sys

from .analyze import main

if __name__ == "__main__":
    sys.exit(main())