import pandas as pd
import numpy as np
import argparse

def get_ratio_from_csv(csv_file, ke, cos):
    # Load CSV data
    df = pd.read_csv(csv_file)
    
    # Find the closest KE and cos bin
    ke_bins = np.sort(df["pi0_KE"].unique())  # Unique sorted KE bins
    cos_bins = np.sort(df["pi0_cos"].unique())  # Unique sorted cos bins

    ke_bin = ke_bins[np.argmin(np.abs(ke_bins - ke))]  # Find closest KE bin
    cos_bin = cos_bins[np.argmin(np.abs(cos_bins - cos))]  # Find closest cos bin

    # Retrieve the ratio for the closest bins
    row = df[(df["pi0_KE"] == ke_bin) & (df["pi0_cos"] == cos_bin)]
    
    if not row.empty:
        return row["Ratio"].values[0]  # Return the corresponding ratio
    else:
        return None  # Return None if no match is found

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(description="Get ratio from CSV based on KE and cos values.")
    parser.add_argument('csv_file', type=str, help="Path to the CSV file")
    parser.add_argument('--ke', type=float, required=True, help="The KE value")
    parser.add_argument('--cos', type=float, required=True, help="The cos value")

    # Parse arguments
    args = parser.parse_args()

    # Call the function to get the ratio
    ratio = get_ratio_from_csv(args.csv_file, args.ke, args.cos)
    
    if ratio is not None:
        print(f"Ratio for KE={args.ke}, cos={args.cos}: {ratio}")
    else:
        print(f"No data found for KE={args.ke}, cos={args.cos}.")

if __name__ == "__main__":
    main()
