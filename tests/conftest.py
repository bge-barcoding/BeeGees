"""Shared fixtures and path setup for BeeGees tests."""
import sys
from pathlib import Path

# Make workflow scripts importable without installing them as a package
SCRIPTS_DIR = Path(__file__).parent.parent / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
