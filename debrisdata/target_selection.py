from data import get_merged_data

data = get_merged_data()

print(f"Before country filter: {len(data)}")

######################
### country filter ###
######################

# Allowed ownership/operator codes
ALLOWED_COUNTRIES = {
    "GER","NETH","UK","SWED","NOR","EUTE","EUME"
}
# Name of the dataframe column containing country/operator codes
COUNTRY_CODE = "COUNTRY_CODE"

# Check column exists
if COUNTRY_CODE in data.columns:

    # Clean formatting
    data[COUNTRY_CODE] = (data[COUNTRY_CODE].astype(str).str.upper().str.strip())

    # Apply filter
    data = data[(data[COUNTRY_CODE].isin(ALLOWED_COUNTRIES))]

    print(f"After country filter: {len(data)}")

else:
    print(f"{COUNTRY_CODE} column not found")


print(data["OBJECT_ID"])
print(data.loc[628])

#Engine compliance with capture mechanism


######################
#### Optimization ####
######################


