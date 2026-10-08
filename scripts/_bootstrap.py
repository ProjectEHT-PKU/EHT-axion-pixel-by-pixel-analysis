"""Make local reusable modules available to direct scripts without installation.

Python puts scripts/ on sys.path when a file is executed directly. Adding its
repository parent works from any working directory; no editable package required.
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
