"""Shared fixtures and path setup for BeeGees tests."""
import sys
from pathlib import Path
from beegees.utils.configs import get_package_dir

# Resolve scripts from the installed beegees package (single source of truth)
# so tests validate exactly what ships, not a copy in the source tree.
SCRIPTS_DIR = get_package_dir() / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
