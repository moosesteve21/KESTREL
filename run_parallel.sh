#!/usr/bin/env bash
# run_parallel.sh - search a target list with several workers sharing one queue.
# The work is done by run_parallel.py (see its header for details); this wrapper
# keeps the familiar command and picks the project's Python.
#
# Usage:
#     ./run_parallel.sh targets/<name>.txt results/<name> [n_workers] [script]
#     (defaults: 8 workers, find_planets.py)
#     MAX_ATTEMPTS (default 3), RETRY_WAIT (seconds before an errored star is
#     retried, default 300) and STAGGER (seconds between worker starts, default 20)
#     can be set in the environment.
#
# Exit status: 0 complete; 3 incomplete (stars in <out>/failed_targets.txt; rerun
# to retry); 4 the folder is busy; 6 stopped (<out>/STOP or interrupted; rerun to
# continue); 1 bad arguments.
here=$(cd "$(dirname "$0")" && pwd)
py="$here/.venv/bin/python"     # the project venv, so it works in a fresh tmux shell
[[ -x "$py" ]] || py=python3
exec "$py" "$here/run_parallel.py" "$@"
