import pandas as pd
import os

def combine_parquet_files(input_files, output_filename):
    """
    Combines a list of Parquet files into a single Parquet file and prints row counts.
    
    Args:
        input_files (list): A list of file paths to the Parquet files.
        output_filename (str): The name/path of the final combined Parquet file.
    """
    dataframes = []
    expected_total_rows = 0
    
    print("--- Individual File Counts ---")
    # Loop through each file, load it, and count the rows
    for file in input_files:
        if os.path.exists(file):
            # Read the parquet file
            df = pd.read_parquet(file)
            row_count = len(df)
            
            print(f"File: '{file}' -> {row_count} rows")
            
            expected_total_rows += row_count
            dataframes.append(df)
        else:
            print(f"Warning: File '{file}' not found. Skipping.")
            
    if not dataframes:
        print("No valid dataframes found to combine.")
        return

    # Combine all individual dataframes into one
    combined_df = pd.concat(dataframes, ignore_index=True)
    actual_total_rows = len(combined_df)
    
    print("\n--- Final Summary ---")
    print(f"Expected Total Rows: {expected_total_rows}")
    print(f"Combined Total Rows: {actual_total_rows}")
    
    # Verify no rows were lost during concatenation
    if expected_total_rows == actual_total_rows:
        print("Success! Row counts match perfectly.")
    else:
        print("Warning: Row counts do not match. Check your data for index overlaps or data loss.")

    # Save the combined dataframe to a new Parquet file
    combined_df.to_parquet(output_filename, index=False, engine='pyarrow')
    print(f"\nCombined Parquet file successfully saved as: '{output_filename}'")


# ==========================================
# Example Usage:
# ==========================================
if __name__ == "__main__":
    # Replace these with your actual file paths
    files_to_combine = [
        "processed\Al_50.parquet",
        "processed\O_80.parquet", 
        "processed\Si_80.parquet",
        "processed\Ti_80.parquet"
    ]
    
    output_file = "processed/fits_data.parquet"
    
    combine_parquet_files(files_to_combine, output_file)