import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # the repo root, so `import solidwright` works
sys.path.insert(0, HERE)                   # so tests can `import shapes, rhino_mock`

import rhino_mock  # noqa: E402

rhino_mock.install()
