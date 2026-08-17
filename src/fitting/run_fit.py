import numpy as np
import pandas as pd
from multiprocessing import Pool
from tqdm import tqdm
import os
from xrf_fitter import XRFFitter


def fit_single_row(args):
    counts, footprint = args
    try:
        fitter = XRFFitter(counts=counts, std_dev=0.1)
        solution = fitter.fit(method='leastsq')

        result = {
            'footprint': footprint,
            'success': True,
            'cost': fitter.result.cost if hasattr(fitter.result, 'cost') else np.sum(fitter.result.fun**2),
        }

        for i, el in enumerate(fitter.elements):
            result[f'conc_{el}'] = solution[i]

        result['scale'] = solution[-2]
        result['std_dev'] = solution[-1]
        return result

    except Exception:
        return {
            'footprint': footprint,
            'success': False,
            'cost': np.nan,
            'conc_Fe': np.nan, 'conc_Al': np.nan, 'conc_Mg': np.nan,
            'conc_Si': np.nan, 'conc_Ca': np.nan, 'conc_Ti': np.nan,
            'conc_O': np.nan, 'scale': np.nan, 'std_dev': np.nan
        }


def run_all_fits(input_path, output_format='csv', max_workers=8, output_dir=None):
    df = pd.read_parquet(input_path)
    n_rows = len(df)

    base_name = os.path.splitext(os.path.basename(input_path))[0]
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"{base_name}_result.{output_format}")
    else:
        output_path = f"{base_name}_result.{output_format}"

    print(f"Found {n_rows} spectra to fit")
    print(f"Using {max_workers} processes")
    print(f"Results will be saved to: {output_path}")

    args_list = [(df['counts'].iloc[i], df['footprint'].iloc[i]) for i in range(n_rows)]

    with Pool(processes=max_workers) as pool:
        results = list(tqdm(pool.imap(fit_single_row, args_list), total=n_rows, desc="Fitting spectra", unit="spectrum"))

    result_df = pd.DataFrame(results)

    completed = result_df['success'].sum()
    failed = n_rows - completed

    if output_format == 'parquet':
        result_df.to_parquet(output_path, index=False)
    else:
        result_df.to_csv(output_path, index=False)

    print(f"\nDone! Completed: {completed}, Failed: {failed}")
    print(f"Results saved to: {output_path}")


if __name__ == "__main__":
    run_all_fits('catalogue_Al_50_2025_10.parquet', output_format='parquet', max_workers=8)
