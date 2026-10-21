"""Run the inversion-only multi-species study using the existing Modules 1-8."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ibp.production.run_study import main  # noqa: I001 - must follow sys.path setup


if __name__ == "__main__":
    sys.exit(main())
