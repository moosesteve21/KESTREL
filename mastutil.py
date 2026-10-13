"""
mastutil.py - a shorter wait for MAST queries in the per-star scripts.

astroquery waits for a MAST query to finish (polling while the server says
"EXECUTING") for up to astroquery.mast.conf.timeout seconds, 600 by default. In
batch b008 two waves of MAST trouble left 24 stars waiting 600 s per attempt (about
200 worker-minutes). A per-star query normally finishes in about 5 s, so these
scripts give up after MAST_TIMEOUT seconds and let their own retries (and the
run's queue) try again.

astroquery copies the setting into its query objects when it is imported, so
changing conf.timeout afterwards has no effect; set_mast_timeout() sets it on the
objects themselves. Not used by tess_target_finder.py: its whole-sky catalogue
queries are large and can legitimately take minutes.
"""
MAST_TIMEOUT = 90       # seconds: polling wait, and silence allowed on one HTTP read
CONNECT_TIMEOUT = 15    # seconds to open a connection


def set_mast_timeout(seconds=MAST_TIMEOUT):
    """Shorten astroquery's polling wait, and give every MAST HTTP request a
    timeout. astroquery sends its requests with no timeout and checks its polling
    limit only between requests, so one request that hangs (as in b008 and on
    2026-10-12 15:00, exactly 600 s each time, apparently dropped by the server)
    isn't cut short by the polling limit alone. Downloads stream, so the read
    timeout applies between chunks, not to a whole large file."""
    import functools
    from astroquery import mast
    from astroquery.mast.discovery_portal import PortalAPI
    PortalAPI.TIMEOUT = seconds
    sessions = set()
    for name in ("Observations", "Catalogs", "Tesscut", "Zcut", "MastMissions"):
        obj = getattr(mast, name, None)
        if obj is None:
            continue
        if getattr(obj, "_session", None) is not None:
            sessions.add(obj._session)
        for attr in ("_portal_api_connection", "_service_api_connection"):
            conn = getattr(obj, attr, None)
            if conn is not None:
                conn.TIMEOUT = seconds
                if getattr(conn, "_session", None) is not None:
                    sessions.add(conn._session)
    for sess in sessions:
        if getattr(sess.request, "_mast_timeout", False):
            continue
        orig = sess.request

        @functools.wraps(orig)
        def request(method, url, *args, _orig=orig, **kw):
            if kw.get("timeout") is None:
                kw["timeout"] = (CONNECT_TIMEOUT, seconds)
            return _orig(method, url, *args, **kw)
        request._mast_timeout = True
        sess.request = request
    return seconds


def retry(func, what="MAST query", tries=3, wait=60):
    """Call func(), retrying transient failures (timeouts, dropped connections)
    up to tries times, wait seconds apart. With the shorter timeout a slow MAST
    spell raises errors sooner, so every per-star MAST call should go through
    this (or its own retry loop)."""
    import logging
    import time
    for attempt in range(1, tries + 1):
        try:
            return func()
        except Exception as e:
            if attempt == tries:
                raise
            logging.warning(f"  {what} failed ({type(e).__name__}); retrying in {wait} s")
            time.sleep(wait)

