from data import get_merged_data

data = get_merged_data()

# Mission parameters
ALT_MIN_KM      = 35_700          # km
ALT_MAX_KM      = 36_500          # km  (includes GEO + graveyard)
INC_MIN_DEG     = 0.0
INC_MAX_DEG     = 20.0
MASS_MIN_KG     = 1_000/0.7
MASS_MAX_KG     = 2_000/0.7
MAX_ECC = 0.001

data = data[(data['MASS_KG'] >= 1500) & (data['MASS_KG'] <= 2900)]
data = data[(data['ECCENTRICITY']<=0.001)]

######################
### country filter ###
######################

# Allowed ownership/operator codes
ALLOWED_COUNTRIES = {
    "US","JPN","GER","NETH","UK","SWED","NOR","EUTE","SES","ITSO",
}
# Name of the dataframe column containing country/operator codes
COUNTRY_CODE = "COUNTRY_CODE"

# Check column exists
if COUNTRY_CODE in data.columns:

    # Clean formatting
    data[COUNTRY_CODE] = (
        data[COUNTRY_CODE]
        .astype(str)
        .str.upper()
        .str.strip()
    )

    # Apply filter
    data = data[(data[COUNTRY_CODE].isin(ALLOWED_COUNTRIES))]

    print(f"After country filter: {len(data)}")

else:
    print(f"{COUNTRY_CODE} column not found")

print(data.columns.tolist())