import os
import glob
import pandas as pd

def process_results(
    source_folder="Al_50_fit/Al-parquets",
    output_folder="processed",
    output_filename="Al_50.parquet",
):
    """
    Loads all *_result.parquet files from the source folder,
    concatenates them, and saves the result as a single Parquet in the output folder.
    """
    # Create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Created output directory: {output_folder}")

    # Find all result parquet files
    search_pattern = os.path.join(source_folder, "*_result.parquet")
    files = sorted(glob.glob(search_pattern))

    if not files:
        print(f"No files found matching {search_pattern}")
        return

    print(f"Found {len(files)} result files. Loading...")

    dfs = []
    for f in files:
        try:
            df = pd.read_parquet(f)
            dfs.append(df)
            print(f"Loaded {os.path.basename(f)} - Shape: {df.shape}")
        except Exception as e:
            print(f"Error loading {f}: {e}")

    if not dfs:
        print("No data loaded.")
        return

    # Concatenate all dataframes
    combined_df = pd.concat(dfs, ignore_index=True)
    print(f"Combined data shape: {combined_df.shape}")

    # Save to Parquet
    output_path = os.path.join(output_folder, output_filename)
    combined_df.to_parquet(output_path, index=False, compression="snappy")
    print(f"Saved processed data to {output_path}")

    # Print final dataframe "structure"
    print("\nFinal dataframe summary:")
    print(f"- Shape (rows, cols): {combined_df.shape}")
    print(f"- Columns ({len(combined_df.columns)}): {list(combined_df.columns)}")
    print("- Preview (first 2 rows):")
    try:
        print(combined_df.head(2).to_string(index=False))
    except Exception:
        print(combined_df.head(2))

if __name__ == "__main__":
    process_results()

'''
Example result:
footprint  success       cost  conc_Fe  conc_Al  conc_Mg  conc_Si  conc_Ca  conc_Ti   conc_O    scale  std_dev
 [[74.8803, 67.8275], [73.5856, 68.0755], [73.6073, 71.2391], [74.9037, 71.249]]     True 198.447425 0.000512 0.903656 0.216553 0.040885 0.000226 0.000276 7.103911 0.000020   0.1099
[[74.0747, 67.9827], [72.7776, 68.2078], [72.7983, 71.2353], [74.0969, 71.2439]]     True 184.602053 0.000471 0.890003 0.234446 0.038601 0.000183 0.000214 7.015423 0.000021   0.1099
'''