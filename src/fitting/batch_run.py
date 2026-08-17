import glob
import os

from run_fit import run_all_fits


def main(
    input_folder="../Si_80",
    output_folder="Si_80_fit",
    output_format="parquet",
    max_workers=8,
):
    os.makedirs(output_folder, exist_ok=True)

    pattern = os.path.join(input_folder, "*.parquet")
    parquet_files = sorted(glob.glob(pattern))

    if not parquet_files:
        print(f"No parquet files found in '{input_folder}'. Nothing to do.")
        return

    print(f"Found {len(parquet_files)} parquet files in '{input_folder}'.")
    print(f"Results will be written to '{output_folder}'.")

    for path in parquet_files:
        print(f"\nProcessing file: {path}")
        run_all_fits(
            input_path=path,
            output_format=output_format,
            max_workers=max_workers,
            output_dir=output_folder,
        )

    print("\nAll files processed.")

if __name__ == "__main__":
    main()
