# first line: 139
@memory.cache
def _cached_fetch_discos(norad_ids_tuple: tuple) -> pd.DataFrame:
    return fetch_data_discos(list(norad_ids_tuple))
