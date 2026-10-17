import os
import sys

# Make the repo root importable so `web.backend.control` is a valid path.
# conftest.py lives in web/backend/tests/, so up three levels = repo root.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
