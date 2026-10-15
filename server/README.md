# Server setup: sessions that survive a reboot

At boot, the `tess-sessions` systemd user service recreates the tmux sessions this
project uses:

| Session | Runs |
|---|---|
| `claude` | Claude Code in `~/TESS`. It starts a fresh session; `/resume` picks up an earlier one |
| `db` | Datasette (`tessdb.py serve`), http://127.0.0.1:8001 |
| `jupyter` | JupyterLab, port 8888 |
| `tess` | a shell for batches, which first runs `batch_status.sh` |

Each command is typed into a normal login shell, so the shell stays if the command
stops. To restart a server: attach, Ctrl-C, Up, Enter.

**Batches are never restarted automatically.** After a reboot, `tmux attach -t tess`
shows any batch that didn't finish and the exact command to resume it. The search
continues where it stopped (`queue.db`), so nothing is searched twice.

## Setup (once)

```bash
sudo loginctl enable-linger sprice      # user services start at boot, without a login
mkdir -p ~/.config/systemd/user
cp ~/TESS/server/tess-sessions.service ~/.config/systemd/user/
systemctl --user daemon-reload
systemctl --user enable tess-sessions.service
```

After editing `tess-sessions.service`, copy it again and run `daemon-reload`. Changes
to `start_sessions.sh` need nothing more.

## Day to day

```bash
tmux ls                                  # the sessions
tmux attach -t claude                    # (Ctrl-b d to detach)
~/TESS/server/batch_status.sh            # unfinished batches, any time
systemctl --user status tess-sessions    # did the boot start work?
journalctl --user -u tess-sessions       # its log
```

- **Don't stop or restart the service while working.** The tmux server and everything
  in it, including a running batch, belong to the service. Stopping it closes every
  session, just as a shutdown would.
- **Undo:** `systemctl --user disable tess-sessions`, and
  `sudo loginctl disable-linger sprice`.
- **Testing the script** without touching the real sessions:
  `TMUX_SOCKET=test server/start_sessions.sh`, look with `tmux -L test ls`, then
  `tmux -L test kill-server`.

## Discord notifications

`pipeline.sh` posts one line to a Discord channel when a batch ends, mentioning you so
Discord alerts your phone and PC. For example: "b013 complete: 500 stars, 9 signals,
1 clean, 3 for VSX (after 3h 12m)". It also posts if the batch was stopped, is
incomplete, or a step failed.

- **Settings:** in `~/.config/planetfinder/discord` (mode 600, outside git):
  `DISCORD_WEBHOOK=https://discord.com/api/webhooks/...` and `DISCORD_USER_ID=...`. Anyone
  with the webhook URL can post to the channel. If it leaks, delete the webhook in
  Discord (Server settings, Integrations) and make a new one.
- **Send a message by hand:** `server/notify.sh "text"`. With no settings file, or if
  Discord is unreachable, it does nothing, so it never stops a pipeline.
- **Test without sending:** `NOTIFY_DRYRUN=1` prints the message instead.

## Discord channels (2026-10-15)

`server/notify.sh --kind search|batch|prefetch "message"` posts to the webhook for that kind,
from `~/.config/planetfinder/discord` (mode 600, never in git): DISCORD_SEARCH_NOTIFICATION_WEBHOOK
("search done" and search-stage problems, from pipeline.sh), DISCORD_BATCH_NOTIFICATION_WEBHOOK
("complete" / "follow-up FAILED", from followup.sh), DISCORD_PREFETCH_NOTIFICATION_WEBHOOK
(prefetch summaries, from prefetch.py). A kind without its own webhook uses DISCORD_WEBHOOK.
`NOTIFY_DRYRUN=1` prints the message and its kind instead of sending.
