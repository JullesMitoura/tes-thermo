"""
    Returns the absolute path to a resource based on a relative path,
    relative to the current working directory.

    This is useful for locating resources such as external binaries,
    configuration files, or solvers that are distributed with the library
    or expected to be in a known location during runtime.

    Parameters:
        relative_path (str): The relative path to the resource.

    Returns:
        str: The absolute path to the resource.
"""

import sys
import os

def path(path: str) -> str:
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, path)