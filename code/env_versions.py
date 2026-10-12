# Writes results/environment.txt: interpreter and solver versions of the reproduction (referee item 12).
import sys, platform, subprocess
mods = ['numpy', 'scipy', 'cvxpy', 'clarabel', 'scs', 'mpmath', 'matplotlib']
lines = [f"python {sys.version.split()[0]} ({platform.platform()})"]
for m in mods:
    try: mod = __import__(m); lines.append(f"{m} {getattr(mod, '__version__', 'unknown')}")
    except Exception as e: lines.append(f"{m} not importable ({e})")
try: lines.append("git " + subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd='/home/user/Fable-Research').decode().strip())
except Exception: pass
open('/home/user/Fable-Research/results/environment.txt', 'w').write("\n".join(lines) + "\n"); print("\n".join(lines))
