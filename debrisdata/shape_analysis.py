import matplotlib.pyplot as plt
import pandas as pd
from data import get_merged_data

SELECTED_PARAMS = ["NORAD_CAT_ID", "DISCOS_NAME", "SHAPE",
                   "WIDTH", "HEIGHT", "DEPTH", "DIAMETER", "SPAN",
                   "xSectMin", "xSectAvg", "xSectMax",
                   "MASS_KG"]


if __name__ == "__main__":
    data = get_merged_data()

    # Drop rows without cross-section data
    cs = data[SELECTED_PARAMS].copy()
    cs["xSectAvg"] = pd.to_numeric(cs["xSectAvg"], errors="coerce")
    cs = cs.dropna(subset=["xSectAvg"])

    # Make a histogram of average cross-section.
    plt.figure()
    plt.hist(cs["xSectAvg"], bins=30, edgecolor="black")
    plt.xlabel("Average Cross-Section (m²)")
    plt.ylabel("Count")
    plt.title("Distribution of Average Cross-Section")
    plt.tight_layout()
    plt.show()

    # Get the average of the average cross sections.
    mean_xsect = cs["xSectAvg"].mean()
    print(f"Mean average cross-section: {mean_xsect:.4f} m²")

    # Retrieve the 5 spacecraft that are closest to the average, but have unique shapes.
    cs = cs.copy()
    cs["dist_from_mean"] = (cs["xSectAvg"] - mean_xsect).abs()
    cs_sorted = cs.sort_values("dist_from_mean")

    selected = []
    seen_shapes = set()
    for _, row in cs_sorted.iterrows():
        shape = str(row["SHAPE"]).strip().lower()
        if shape not in seen_shapes:
            seen_shapes.add(shape)
            selected.append(row)
        if len(selected) == 5:
            break

    result = pd.DataFrame(selected)[SELECTED_PARAMS + ["dist_from_mean"]]
    print("\n5 spacecraft closest to mean cross-section with unique shapes:")
    print(result.to_string(index=False))
