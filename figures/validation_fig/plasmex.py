"""
Where the validation scripts find PlasmeX. Nothing from PlasmeX is copied into this repository:
scripts and data (bohm.py, almeida_arc_bc.py, the transport tables) are always loaded from
the live source tree, so later changes to PlasmeX are picked up automatically.

    PlasmeX root:  $PLASMEX_ROOT, else ~/work/PlasmeX
    scripts/data:  <root>/exe

    from plasmex import exe_module
    bohm = exe_module("bohm")
"""
import importlib.util
import os
import sys

ROOT = os.path.abspath(os.path.expanduser(os.environ.get("PLASMEX_ROOT", "~/work/PlasmeX")))
EXE = os.path.join(ROOT, "exe")


def exe_module(name):
    """Import PlasmeX/exe/<name>.py (once) under its own name."""
    if name in sys.modules:
        return sys.modules[name]
    path = os.path.join(EXE, name + ".py")
    if not os.path.isfile(path):
        raise SystemExit(f"cannot find {path} -- set PLASMEX_ROOT to the PlasmeX source tree")
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod
