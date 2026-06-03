from data import get_merged_data

data = get_merged_data()

print(f"Before country filter: {len(data)}")

######################
### country filter ###
######################

# Allowed ownership/operator codes
ALLOWED_COUNTRIES = {
    "GER","NETH","SWED","NOR","EUTE","EUME","US","SES","ITSO","CA","UK","AUS","JPN","IND"
}
# Name of the dataframe column containing country/operator codes
COUNTRY_CODE = "COUNTRY_CODE"

#Engine compliance with capture mechanism
remove_rows = [504,520,525,758,960,1002,1216] #Engines are not compatible
data = data.drop(remove_rows)

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


#Generate table of desired parameters
table = [
    (
        int(idx),
        row["OBJECT_NAME"],
        float(row["MASS_KG"]),
        float(row["SEMIMAJOR_AXIS"]),
        float(row["INCLINATION"]),
        float(row["RA_OF_ASC_NODE"])
    )
    for idx, row in data.iterrows()
]

print(table)