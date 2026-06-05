import json
import logging
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from matplotlib import pyplot as plt
from spacetrack import SpaceTrackClient
import spacetrack.operators as op
from plotting import plot_raan_hist, plot_inclination_hist
from dotenv import load_dotenv

# CONFIG
load_dotenv(Path(__file__).parent.parent / '.env')
ST_USER = os.environ.get("ST_USER")
ST_PASS = os.environ.get("ST_PASS")
ST_BASE = "https://space-track.org/"
DISCOS_TOKEN    = os.environ.get("DISCOS_TOKEN")   # https://discosweb.esac.esa.int

# Mission parameters
ALT_MIN_KM      = 35_700          # km
ALT_MAX_KM      = 36_500          # km  (includes GEO + graveyard)
INC_MIN_DEG     = 0.0
INC_MAX_DEG     = 15.0
MASS_MIN_KG     = 1_000/0.7
MASS_MAX_KG     = 2_000/0.7

# Cache config
# Cache lives inside the repo so it can be shared across users via git.
# Commit the contents of CACHE_DIR to share refreshed data with teammates.
# Note on concurrency: if two people refresh and push at the same time, the
# second push will need to be rebased. With a 30-day TTL this should be rare.
CACHE_DIR           = Path(__file__).parent.parent / "cache" / "satellite_data"
CACHE_MAX_AGE_HOURS = 24*30        # Re-fetch if cached file is older than this

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)


def _cache_path(name: str) -> Path:
    """Return the on-disk path for a named cache entry."""
    return CACHE_DIR / f"{name}.pkl.gz"


def _cache_is_fresh(name: str) -> bool:
    """Return True if the cache file for *name* exists and is within CACHE_MAX_AGE_HOURS."""
    path = _cache_path(name)
    if not path.exists():
        return False
    age_hours = (time.time() - path.stat().st_mtime) / 3600
    if age_hours > CACHE_MAX_AGE_HOURS:
        log.info("Cache for '%s' is %.1f h old (max %s h) — will re-fetch.",
                 name, age_hours, CACHE_MAX_AGE_HOURS)
        return False
    log.info("Cache for '%s' is %.1f h old — using cached data.", name, age_hours)
    return True


def _read_cache(name: str) -> pd.DataFrame:
    return pd.read_pickle(_cache_path(name), compression="gzip")


def _write_cache(name: str, df: pd.DataFrame) -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    df.to_pickle(_cache_path(name), compression="gzip")


def clear_cache():
    """Delete all cached data files."""
    if CACHE_DIR.exists():
        for f in CACHE_DIR.glob("*.pkl.gz"):
            f.unlink(missing_ok=True)
    log.info("Cache cleared: %s", CACHE_DIR)


def fetch_data_startrack() -> pd.DataFrame:
    log.info("Fetching SpaceTrack GP data ...")
    st = SpaceTrackClient(ST_USER, ST_PASS)

    st_data_raw = st.gp(
        orderby="NORAD_CAT_ID",
        inclination=op.inclusive_range(INC_MIN_DEG, INC_MAX_DEG),
        periapsis=op.greater_than(ALT_MIN_KM),
        apoapsis=op.less_than(ALT_MAX_KM),
        format="json",
    )
    st_data_json = json.loads(st_data_raw)

    df = pd.DataFrame(st_data_json)
    df = df.apply(pd.to_numeric, errors="ignore")
    log.info("SpaceTrack returned %d records.", len(df))
    return df


def fetch_data_discos(norad_ids) -> pd.DataFrame:
    URL = "https://discosweb.esoc.esa.int/api/objects"
    BATCH = 30
    records = []

    log.info("Fetching DISCOSweb data for %d NORAD IDs ...", len(norad_ids))
    for i in range(0, len(norad_ids), BATCH):
        batch_ids = norad_ids[i: i + BATCH]
        headers = {
            "Authorization": f"Bearer {DISCOS_TOKEN}",
            "Accept": "application/json",
        }
        satnostr = ",".join(map(str, batch_ids))
        params = {"filter": f"in(satno,({satnostr}))"}

        r = requests.get(URL, headers=headers, params=params)
        r.raise_for_status()
        for obj in r.json().get("data", []):
            attrs = obj.get("attributes", {})
            records.append({
                "SATNO": attrs.get("satno", ""),
                "DISCOS_NAME": attrs.get("name", ""),
                "MASS_KG": attrs.get("mass"),
                "OBJECT_CLASS": attrs.get("objectClass", ""),
                "SHAPE": attrs.get("shape", ""),
                "WIDTH": attrs.get("width", ""),
                "HEIGHT": attrs.get("height", ""),
                "DEPTH": attrs.get("depth", ""),
                "DIAMETER": attrs.get("diameter", ""),
                "xSectMax": attrs.get("xSectMax", ""),
                "xSectMin": attrs.get("xSectMin", ""),
                "xSectAvg": attrs.get("xSectAvg", ""),
                "SPAN": attrs.get("span", ""),
                "ACTIVE": attrs.get("active", ""),
            })

    df = pd.DataFrame(records)
    df = df.apply(pd.to_numeric, errors="ignore")
    log.info("DISCOSweb returned %d records.", len(df))
    return df


def get_startrack_data(force_refresh: bool = False) -> pd.DataFrame:
    """Return SpaceTrack data, using cache unless stale or *force_refresh*."""
    if not force_refresh and _cache_is_fresh("startrack"):
        return _read_cache("startrack")
    df = fetch_data_startrack()
    _write_cache("startrack", df)
    return df


def get_discos_data(norad_ids, force_refresh: bool = False) -> pd.DataFrame:
    """Return DISCOSweb data, using cache unless stale or *force_refresh*."""
    if not force_refresh and _cache_is_fresh("discos"):
        return _read_cache("discos")
    df = fetch_data_discos(list(norad_ids))
    _write_cache("discos", df)
    return df


def get_merged_data(force_refresh: bool = False) -> pd.DataFrame:

    data_st = get_startrack_data(force_refresh=force_refresh)
    data_discos = get_discos_data(data_st["NORAD_CAT_ID"], force_refresh=force_refresh)

    data = data_st.merge(data_discos, left_on="NORAD_CAT_ID", right_on="SATNO", how="inner")

    # Filter data to be in the mass range
    mask = data["MASS_KG"].between(MASS_MIN_KG, MASS_MAX_KG, inclusive='both')
    data = data[mask]

    return data


if __name__ == "__main__":
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    data = get_merged_data(force_refresh=False)

    # cluster(data)

    plot_raan_hist(data)
    plt.show()
    plot_inclination_hist(data)
    plt.show()

    # missing = set(data_st["NORAD_CAT_ID"]) - set(data_discos["SATNO"])
    # print(len(missing), list(missing)[:10])
    #TODO: Find out why so many are missing from discosweb. Seems like discosweb returns same object for multiple satno