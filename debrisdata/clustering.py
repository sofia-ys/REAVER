from sklearn.cluster import k_means
import pandas as pd
from plotting import plot_raan_hist, plot_inclination_hist

def cluster(data: pd.DataFrame):
    data = data.set_index('SATNO')
    fit_data = data[['INCLINATION', 'RA_OF_ASC_NODE']].copy()

    centroid, labels, inertia = k_means(fit_data, n_clusters=10, n_init='auto', random_state=0)


    fit_data['CLUSTER'] = labels
    largest_label = fit_data['CLUSTER'].value_counts().idxmax()
    data_largest = fit_data.loc[fit_data['CLUSTER'] == largest_label]

    plot_raan_hist(data_largest, show=True)
    plot_inclination_hist(data_largest, show=True)


    fit_data