
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

GEO_PERI_MIN_KM = 35700
GEO_PERI_MAX_KM = 36500
GEO_APO_MIN_KM = 35700
GEO_APO_MAX_KM = 36500
INC_MIN_DEG = 0
INC_MAX_DEG = 20


def find_default_spacetrack_csv() -> Path:
    candidates = sorted(Path.cwd().glob("st_*.csv"))
    if not candidates:
        candidates = sorted(Path.cwd().glob("*.csv"))
    if not candidates:
        raise FileNotFoundError(
            "No CSV file found. Pass it explicitly with --spacetrack path/to/file.csv"
        )
    return candidates[0]


def clean_satcat_id(series: pd.Series) -> pd.Series:
    """Normalise NORAD/Satcat IDs by removing spaces and leading zeroes."""
    return (
        series.astype(str)
        .str.strip()
        .str.replace(r"\.0$", "", regex=True)
        .str.lstrip("0")
        .replace("", pd.NA)
    )


def read_spacetrack_csv(path: Path) -> pd.DataFrame:
    st = pd.read_csv(path, dtype=str)
    st.columns = st.columns.str.strip()

    required = {"NORAD_CAT_ID", "OBJECT_NAME"}
    missing = required - set(st.columns)
    if missing:
        raise KeyError(f"Space-Track CSV is missing required columns: {missing}")

    st["NORAD_KEY"] = clean_satcat_id(st["NORAD_CAT_ID"])

    # Convert useful numeric columns when present.
    for col in [
        "PERIAPSIS",
        "APOAPSIS",
        "INCLINATION",
        "RA_OF_ASC_NODE",
        "ECCENTRICITY",
        "ARG_OF_PERICENTER",
        "MEAN_ANOMALY",
        "MEAN_MOTION",
        "SEMIMAJOR_AXIS",
        "PERIOD",
    ]:
        if col in st.columns:
            st[col] = pd.to_numeric(st[col], errors="coerce")

    # Re-apply REAVER orbital filter if the fields are present.
    if {"PERIAPSIS", "APOAPSIS", "INCLINATION"}.issubset(st.columns):
        st = st[
            st["PERIAPSIS"].between(GEO_PERI_MIN_KM, GEO_PERI_MAX_KM)
            & st["APOAPSIS"].between(GEO_APO_MIN_KM, GEO_APO_MAX_KM)
            & st["INCLINATION"].between(INC_MIN_DEG, INC_MAX_DEG)
        ].copy()

    if "OBJECT_TYPE" in st.columns:
        st = st[st["OBJECT_TYPE"].astype(str).str.upper().str.strip() == "PAYLOAD"].copy()

    return st


def read_gcat_currentcat(path: Path) -> pd.DataFrame:
    """
    Read GCAT currentcat.tsv correctly.

    IMPORTANT:
    The first line is the header but starts with '#JCAT'.
    The second line is an update note, e.g. '# Updated ...'.
    Therefore, do NOT use comment='#', because that skips the real header.
    """
    gcat = pd.read_csv(path, sep="\t", skiprows=[1], dtype=str)
    gcat.columns = [c.strip().lstrip("#") for c in gcat.columns]

    if "Satcat" not in gcat.columns:
        raise KeyError(
            "Could not find 'Satcat' in GCAT currentcat.tsv. "
            f"Columns found were: {list(gcat.columns)}"
        )

    gcat["SATCAT_KEY"] = clean_satcat_id(gcat["Satcat"])

    # Normalise text fields used in filters.
    for col in ["Active", "DDate", "ExpandedStatus", "OpOrbit", "Type", "Name"]:
        if col in gcat.columns:
            gcat[col] = gcat[col].fillna("").astype(str).str.strip()

    return gcat


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spacetrack", type=Path, default=None, help="Space-Track CSV file")
    parser.add_argument("--gcat", type=Path, default=Path("currentcat.tsv"), help="Local GCAT currentcat.tsv file")
    parser.add_argument("--outdir", type=Path, default=Path.cwd(), help="Output directory")
    args = parser.parse_args()

    spacetrack_path = args.spacetrack or find_default_spacetrack_csv()
    gcat_path = args.gcat
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    if not spacetrack_path.exists():
        raise FileNotFoundError(f"Space-Track file not found: {spacetrack_path}")
    if not gcat_path.exists():
        raise FileNotFoundError(f"GCAT currentcat.tsv not found: {gcat_path}")

    st = read_spacetrack_csv(spacetrack_path)
    gcat = read_gcat_currentcat(gcat_path)

    enriched = st.merge(
        gcat,
        left_on="NORAD_KEY",
        right_on="SATCAT_KEY",
        how="left",
        suffixes=("_ST", "_GCAT"),
    )

    # GCAT screening logic:
    # Active == 'P' is used here as inactive payload / non-operational payload proxy.
    # DDate '-' or blank means it has not decayed/re-entered according to currentcat.
    # This is a screening filter, not final mission validation.
    ddate_ok = enriched.get("DDate", pd.Series("", index=enriched.index)).fillna("").astype(str).str.strip().isin(["", "-"])
    active_p = enriched.get("Active", pd.Series("", index=enriched.index)).fillna("").astype(str).str.strip().eq("P")
    in_earth_orbit = enriched.get("ExpandedStatus", pd.Series("", index=enriched.index)).fillna("").astype(str).str.contains("Earth orbit", case=False, na=False)

    defunct_candidates = enriched[active_p & ddate_ok & in_earth_orbit].copy()

    missing = enriched[enriched["SATCAT_KEY"].isna()].copy()

    # Save outputs.
    st.to_csv(outdir / "spacetrack_cleaned_geo_inc20.csv", index=False)
    enriched.to_csv(outdir / "spacetrack_gcat_enriched.csv", index=False)
    defunct_candidates.to_csv(outdir / "reaver_defunct_candidates_gcat.csv", index=False)
    missing.to_csv(outdir / "reaver_missing_gcat_matches.csv", index=False)

    print("Done.")
    print(f"Space-Track input file: {spacetrack_path}")
    print(f"GCAT input file:        {gcat_path}")
    print(f"Cleaned Space-Track candidates: {len(st)}")
    print(f"GCAT enriched rows:              {len(enriched)}")
    print(f"Likely defunct payload candidates: {len(defunct_candidates)}")
    print(f"Missing GCAT matches: {len(missing)}")
    print()
    print("Created files:")
    print(f"- {outdir / 'spacetrack_cleaned_geo_inc20.csv'}")
    print(f"- {outdir / 'spacetrack_gcat_enriched.csv'}")
    print(f"- {outdir / 'reaver_defunct_candidates_gcat.csv'}")
    print(f"- {outdir / 'reaver_missing_gcat_matches.csv'}")

    preview_cols = [
        "OBJECT_NAME",
        "NORAD_CAT_ID",
        "INCLINATION",
        "RA_OF_ASC_NODE",
        "PERIAPSIS",
        "APOAPSIS",
        "Active",
        "Name",
        "ExpandedStatus",
        "DDate",
        "OpOrbit",
    ]
    preview_cols = [c for c in preview_cols if c in defunct_candidates.columns]
    if preview_cols and len(defunct_candidates) > 0:
        print()
        print("Preview:")
        print(defunct_candidates[preview_cols].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
