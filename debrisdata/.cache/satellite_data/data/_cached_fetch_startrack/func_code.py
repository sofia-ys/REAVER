# first line: 123
@memory.cache
def _cached_fetch_startrack() -> pd.DataFrame:
    return fetch_data_startrack()
