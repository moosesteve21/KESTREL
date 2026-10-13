#!/usr/bin/env bash
# server/notify.sh - post a one-line message to the project's Discord channel.
#
# Usage: server/notify.sh "b013 complete: 7 signals, 1 clean"
#
# The webhook URL and the user to mention are read from ~/.config/planetfinder/discord
# (DISCORD_WEBHOOK=..., DISCORD_USER_ID=...), kept out of git: anyone with the URL can
# post to the channel. Without that file, or if Discord can't be reached, it does
# nothing and still exits 0, so a notification never stops the pipeline.
conf="${PLANETFINDER_DISCORD:-$HOME/.config/planetfinder/discord}"
[[ -r "$conf" && $# -ge 1 ]] || exit 0
# shellcheck disable=SC1090
source "$conf"
[[ -n "${DISCORD_WEBHOOK:-}" ]] || exit 0
mention=""; allowed='{"parse":[]}'
if [[ -n "${DISCORD_USER_ID:-}" ]]; then
    mention="<@$DISCORD_USER_ID> "
    allowed=$(jq -nc --arg u "$DISCORD_USER_ID" '{users:[$u]}')
fi
payload=$(jq -nc --arg c "$mention$*" --argjson a "$allowed" \
    '{content: ($c | .[0:1900]), username: "PlanetFinder", allowed_mentions: $a}')
[[ -n "${NOTIFY_DRYRUN:-}" ]] && { echo "$payload"; exit 0; }    # testing: print, don't send
curl -fsS --max-time 15 -H "Content-Type: application/json" -d "$payload" "$DISCORD_WEBHOOK" \
    >/dev/null 2>&1 || echo "notify: Discord message not sent" >&2
exit 0
