import importlib.util
import os
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)


@pytest.fixture(scope='session')
def be3d_local():
    """examples/be3d_local.py as a module (its __main__ block, which imports the pipeline, does not run)."""
    spec = importlib.util.spec_from_file_location('be3d_local', os.path.join(REPO, 'examples', 'be3d_local.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
