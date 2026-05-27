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

ALLOWED_ENGINES = {
    "GER","NETH","UK","SWED","NOR","EUTE","EUME"
}

#Engine compliance with capture mechanism
remove_rows = [504,520,525,758,960,1002,1216] #Engines are not compatible
data = data.drop(remove_rows)

#Generate table of desired parameters
table = data[["OBJECT_NAME","MASS_KG","PERIOD","SEMIMAJOR_AXIS","INCLINATION","RA_OF_ASC_NODE"]]
print(table.to_string())