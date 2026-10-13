#!/usr/bin/env bash
# server/start_sessions.sh - recreate the PlanetFinder tmux sessions after a reboot.
# Run at boot by the tess-sessions user service (see server/README.md); safe to run by
# hand too: sessions that already exist are left alone.
#
#   claude   Claude Code in ~/TESS (a fresh session: /resume picks up an old one)
#   db       Datasette, tessdb.py serve (http://127.0.0.1:8001)
#   jupyter  JupyterLab on port 8888
#   tess     a shell for batches; shows any batch a reboot interrupted
#
# Each command is typed into a login shell, as if by hand, so the shell stays when
# the command ends (Ctrl-C a server and press Up, Enter to restart it).
# Batches are never restarted automatically: rerun the command batch_status.sh shows.
set -u
TESS="$HOME/TESS"
tmux="/usr/bin/tmux ${TMUX_SOCKET:+-L $TMUX_SOCKET}"   # TMUX_SOCKET: a separate server, for testing

session() {    # session NAME [COMMAND]
    $tmux has-session -t "=$1" 2>/dev/null && return
    $tmux new-session -d -s "$1" -c "$TESS"
    [[ -n "${2:-}" ]] && $tmux send-keys -t "=$1:" "$2" Enter
}

session db      ".venv/bin/python tessdb.py serve"
session jupyter ".venv/bin/jupyter-lab --no-browser --port 8888"
session tess    "server/batch_status.sh"
session claude  "$HOME/.local/bin/claude"
