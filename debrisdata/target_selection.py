from data import get_merged_data

data = get_merged_data()

# Mission constraints
MIN_MASS_KG = 1500
MAX_MASS_KG = 2900

MAX_ECC = 0.001

MIN_INC = 0
MAX_INC = 20

MIN_ALT = 35700
MAX_ALT = 36500

data = data[(data['MASS_KG'] >= 1500) & (data['MASS_KG'] <= 2900)]
data = data[(data['ECCENTRICITY']<=0.001)]

######################
### country filter ###
######################

# Allowed ownership/operator codes
ALLOWED_COUNTRIES = {
    "US","JPN","GER","NETH","UK","SWED","NOR","EUTE","SES","ITSO",
}

# Replace OWNER with actual column if different
COUNTRY_COLUMN = "COUNTRY"

if COUNTRY_COLUMN in data.columns:

    data[COUNTRY_COLUMN] = (data[COUNTRY_COLUMN].astype(str).str.upper().str.strip())

    data = data[data[COUNTRY_COLUMN].isin(ALLOWED_COUNTRIES)]

    print(f"After country filter: {len(data)}")

else:
    print(f"{COUNTRY_COLUMN} column not found")

######################
### mass filter ###
######################

data = data[
    (data["MASS_KG"] >= MIN_MASS_KG) &
    (data["MASS_KG"] <= MAX_MASS_KG)
]

print(f"After mass filter: {len(data)}")