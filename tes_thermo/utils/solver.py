import os
import sys
import platform
import pyomo.environ as pyo

def path(path: str) -> str:
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, path)

def get_solver(solver_path: str = None):
    if platform.system() == "Windows":
        solver = pyo.SolverFactory('ipopt', executable=path(solver_path))
    else:
        solver = pyo.SolverFactory('ipopt')  # Assumes 'ipopt' is in PATH
    return solver