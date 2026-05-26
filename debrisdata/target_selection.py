from data import get_merged_data

data = get_merged_data()


# Filter for spacecraft within the mass range of 1500-2900 kg
# Assuming 'MASS_KG' is the column name for the spacecraft's mass in kilograms.
# If the column name is different, please let me know.
data = data[(data['MASS_KG'] >= 1500) & (data['MASS_KG'] <= 2900)]
data = data[(data['ECCENTRICITY']<=0.001)]
print(data)
