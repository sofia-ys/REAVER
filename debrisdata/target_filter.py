from data import get_merged_data
import pandas as pd


data = get_merged_data()

print(f"Before country filter: {len(data)}")

######################
### country filter ###
######################

# Allowed ownership/operator codes
ALLOWED_COUNTRIES = {
    "GER","NETH","UK","SWED","NOR","EUTE","EUME","US","SES","ITSO"
}
# Name of the dataframe column containing country/operator codes
COUNTRY_CODE = "COUNTRY_CODE"
#Engine compliance with capture mechanism

#Engine compliance with capture mechanism
# remove_rows = [504,520,525,758,960,1002,1216] #Engines are not compatible
# data = data.drop(remove_rows)


# Check column exists
if COUNTRY_CODE in data.columns:

    # Clean formatting
    data[COUNTRY_CODE] = (data[COUNTRY_CODE].astype(str).str.upper().str.strip())

    # Apply filter
    data = data[(data[COUNTRY_CODE].isin(ALLOWED_COUNTRIES))]

    print(f"After country filter: {len(data)}")

else:
    print(f"{COUNTRY_CODE} column not found")

# Allowed RAAN range in degrees
RAAN_MIN_DEG = 60
RAAN_MAX_DEG = 90
# Name of the dataframe column containing RAAN values
RAAN_COL = "RA_OF_ASC_NODE"

# Check column exists
if RAAN_COL in data.columns:

    # Ensure numeric
    data[RAAN_COL] = pd.to_numeric(data[RAAN_COL], errors="coerce")

    # Apply filter
    data = data[data[RAAN_COL].between(RAAN_MIN_DEG, RAAN_MAX_DEG, inclusive="both")]

    print(f"After RAAN filter: {len(data)}")

else:
    print(f"{RAAN_COL} column not found")


print(data.to_string())
#Generate table of desired parameters
table = [
    (
        int(idx),
        row["OBJECT_NAME"],
        float(row["MASS_KG"]),
        float(row["SEMIMAJOR_AXIS"]),
        float(row["INCLINATION"]),
        float(row["RA_OF_ASC_NODE"]),
        float(row["ECCENTRICITY"]),
        float(row["ARG_OF_PERICENTER"]),
        float(row["MEAN_ANOMALY"]),
        float(row["REV_AT_EPOCH"]),
        float(row["PERIOD"]),
        float(row["APOAPSIS"]),
        float(row["PERIAPSIS"]),
        row["EPOCH"],
        row["LAUNCH_DATE"]
    )
    for idx, row in data.iterrows()
]

print(table)

#epoch, mean motion, eccentricity, arg_of_pericenter, mean_anomaly, REV_AT_EPOCh, Period, APOAPSIS, PERIAPSIS, LAUNCH_DATE