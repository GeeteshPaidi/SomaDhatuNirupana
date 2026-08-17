import os
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from skimage.draw import polygon
import ast

class GaussianArray:
    def __init__(self, grid_size=(64, 64), abundance_map_file=None, coverage_file=None, is_map=False):
        if is_map:
            if abundance_map_file is not None and coverage_file is not None:
                with open(coverage_file, 'rb') as f:
                    count_map = pickle.load(f)
                with open(abundance_map_file, 'rb') as f:
                    abundance_map = pickle.load(f)
                self.grid_size = abundance_map.shape
                self.arr = np.zeros((self.grid_size[0], self.grid_size[1], 2))
                self.arr[:, :, 0] = abundance_map
                self.arr[:, :, 1] = count_map
            else:
                raise ValueError("Abundance map and count file must be provided")
        else:
            self.grid_size = grid_size
            self.arr = np.zeros((grid_size[0], grid_size[1], 2))

    def in_block_or_not(self, img_lat, img_lon, block_lat, block_lon):
        return (block_lat[0] <= min(img_lat) <= block_lat[2] and
                block_lat[0] <= max(img_lat) <= block_lat[2] and
                block_lon[0] <= min(img_lon) <= block_lon[2] and
                block_lon[0] <= max(img_lon) <= block_lon[2])

    def convert_coords_to_indices(self, lat, lon, block_lat, block_lon):
        lat_scale = (self.grid_size[0] - 1) / (block_lat[2] - block_lat[0])
        lon_scale = (self.grid_size[1] - 1) / (block_lon[2] - block_lon[0])

        # FIX: Use dynamic length instead of hardcoded range(4)
        lat_indices = [int((lat[i] - block_lat[0]) * lat_scale) for i in range(len(lat))]
        lon_indices = [int((lon[i] - block_lon[0]) * lon_scale) for i in range(len(lon))]

        return lat_indices, lon_indices

    def calculate_diagonal_length(self, lat_indices, lon_indices):
        # FIX: Ensure it doesn't crash if footprint isn't exactly 4 points
        if len(lat_indices) >= 4:
            diag1 = np.sqrt((lat_indices[2] - lat_indices[0]) ** 2 + (lon_indices[2] - lon_indices[0]) ** 2)
            diag2 = np.sqrt((lat_indices[3] - lat_indices[1]) ** 2 + (lon_indices[3] - lon_indices[1]) ** 2)
            return (diag1 + diag2) / 2
        else:
            # Fallback for arbitrary polygons
            return np.sqrt((max(lat_indices) - min(lat_indices))**2 + (max(lon_indices) - min(lon_indices))**2)

    def generate_gaussian_distribution(self, shape, center, sigma):
        x = np.arange(0, shape[0], 1, float)
        y = np.arange(0, shape[1], 1, float)
        
        # FIX: Use 'ij' indexing so x maps to rows (shape[0]) and y maps to cols (shape[1])
        x, y = np.meshgrid(x, y, indexing='ij')
        
        # FIX: Mathematical max is naturally 1.0; dividing by max() risks ZeroDivisionError on small floats
        gauss = np.exp(-((x - center[0]) ** 2 + (y - center[1]) ** 2) / (2 * sigma ** 2))
        return gauss

    # target_diagonal/base_value: empirically fixed constants calibrated to
    # the CLASS footprint geometry at this pipeline's default grid
    # resolution. No further derivation beyond this fixed, reproducible use
    # is documented or claimed (flagged, not re-derived, per project record).
    def fill_up_the_array(self, img_lat, img_lon, block_lat, block_lon, max_value, target_diagonal=17.625, base_value=2.1739):
        if max_value <= 0 or np.isnan(max_value):
            return
            
        if self.in_block_or_not(img_lat, img_lon, block_lat, block_lon):
            img_lat_indices, img_lon_indices = self.convert_coords_to_indices(img_lat, img_lon, block_lat, block_lon)
            poly_points = np.array([img_lat_indices, img_lon_indices]).T
            min_y, min_x = np.min(poly_points, axis=0)
            max_y, max_x = np.max(poly_points, axis=0)

            avg_diagonal = self.calculate_diagonal_length(img_lat_indices, img_lon_indices)
            if avg_diagonal == 0:
                return 
                
            scale_factor = target_diagonal / avg_diagonal
            sigma = base_value * scale_factor

            # FIX: Ensure height and width are integers
            height, width = int(max_y - min_y + 1), int(max_x - min_x + 1)
            
            if height <= 0 or width <= 0:
                return 
                
            # FIX: Mathematically center the distribution rather than using floor division
            center_x, center_y = (height - 1) / 2.0, (width - 1) / 2.0
            
            gaussian_values = self.generate_gaussian_distribution((height, width), (center_x, center_y), sigma) * max_value
            
            if np.isnan(gaussian_values).any():
                return

            rr, cc = polygon(poly_points[:, 0] - min_y, poly_points[:, 1] - min_x, gaussian_values.shape)

            for r, c in zip(rr, cc):
                x, y = r + min_y, c + min_x
                if 0 <= x < self.grid_size[0] and 0 <= y < self.grid_size[1]:
                    if self.arr[x, y, 1] == 0:  
                        self.arr[x, y, 0] = gaussian_values[r, c]
                        self.arr[x, y, 1] = 1
                    else:
                        count = self.arr[x, y, 1]
                        self.arr[x, y, 0] = (self.arr[x, y, 0] * count + gaussian_values[r, c]) / (count + 1)
                        self.arr[x, y, 1] += 1

    def add_gaussian_box(self, img_lat, img_lon, block_lat, block_lon, max_value, target_diagonal=17.625, base_value=2.1739, plot=False):
        self.fill_up_the_array(img_lat, img_lon, block_lat, block_lon, max_value, target_diagonal, base_value)
        if plot:
            self.visualize_heatmap()

    def visualize_heatmap(self):
        heatmap_data = self.arr[:, :, 0]
        plt.figure(figsize=(10, 8))
        sns.heatmap(heatmap_data, cmap='viridis')
        plt.title('Heatmap of Elemental Abundance with Gaussian Distributions')
        plt.xlabel('Longitude (Pixels)')
        plt.ylabel('Latitude (Pixels)')
        plt.show()
        
    def visualize_counts(self):
        heatmap_data = self.arr[:, :, 1]
        plt.figure(figsize=(10, 8))
        sns.heatmap(heatmap_data, cmap='viridis')
        plt.title('Coverage Map (Number of Overlapping Footprints)')
        plt.xlabel('Longitude (Pixels)')
        plt.ylabel('Latitude (Pixels)')
        plt.show()

    def export_map_pkl(self, f_name, plt_show=False):
        arr = self.arr[:, :, 0]
        with open(f_name, 'wb') as f:
            pickle.dump(arr, f)

    def export_coverage(self, f_name, plt_show=False):
        arr = self.arr[:, :, 1]
        with open(f_name, 'wb') as f:
            pickle.dump(arr, f)

    def export_map_png(self, f_name, dpi=192, cmap='YlOrRd', resize_fact=1, plt_show=False, save_dir='./', element_name="Element"):
        from mpl_toolkits.axes_grid1 import make_axes_locatable 

        arr = self.arr[:, :, 0].copy() 
        counts = self.arr[:, :, 1]
        
        arr[counts == 0] = np.nan
        valid_data = arr[~np.isnan(arr)]
        
        if len(valid_data) > 0:
            non_zero_data = valid_data[valid_data > 0]
            if len(non_zero_data) > 0:
                vmin = np.min(non_zero_data)
            else:
                vmin = np.min(valid_data) 
                
            vmax = np.percentile(valid_data, 98)
            
            # FIX: Prevent Matplotlib ValueError if max and min overlap
            if vmax <= vmin:
                vmax = vmin + 1e-6
        else:
            vmin, vmax = 0, 1
            
        print(f"Applying color limits: Min={vmin:.6f}, Max={vmax:.6f}")

        fig, ax = plt.subplots(figsize=(arr.shape[1] / dpi, arr.shape[0] / dpi), facecolor='black', dpi=dpi)
        ax.set_facecolor('black')
        
        current_cmap = plt.get_cmap(cmap).copy()
        current_cmap.set_bad(color='black')
        
        im = ax.imshow(arr, cmap=current_cmap, origin='lower', vmin=vmin, vmax=vmax, extent=[-180, 180, -90, 90])
        
        ax.set_xlabel('Longitude (Degrees)', color='white', fontsize=10)
        ax.set_ylabel('Latitude (Degrees)', color='white', fontsize=10)
        ax.tick_params(colors='white', labelsize=8)
        
        ax.set_xticks(np.arange(-180, 181, 60))
        ax.set_yticks(np.arange(-90, 91, 30))
        
        ax.set_title(f'Moon Surface Mapping - {element_name}', color='white', fontsize=14, pad=15)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.1) 
        cbar = fig.colorbar(im, cax=cax)
        
        cbar.set_label(f'Relative {element_name} Concentration', color='white', fontsize=10)
        cbar.ax.yaxis.set_tick_params(color='white', labelsize=8)
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
        
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        file_path = os.path.join(save_dir, f_name)
        plt.savefig(file_path, facecolor='black', dpi=(dpi * resize_fact), bbox_inches='tight', pad_inches=0.2)
        
        if plt_show:
            plt.show()
        else:
            plt.close()

    def check_coverage(self):
        coverage = np.sum(self.arr[:, :, 1] > 0) / np.prod(self.grid_size)
        return coverage

    def plot_value_distribution(self, element_name="Element", save_dir='./'):
        arr = self.arr[:, :, 0]
        counts = self.arr[:, :, 1]

        valid_data = arr[counts > 0].flatten()

        if len(valid_data) == 0:
            print("No valid data to plot.")
            return

        vmin = np.percentile(valid_data, 2)
        vmax = np.percentile(valid_data, 98)
        median = np.median(valid_data)

        plt.figure(figsize=(12, 6), facecolor='white')
        
        sns.histplot(valid_data, bins=100, kde=True, color='indigo', edgecolor='black')

        plt.axvline(vmin, color='red', linestyle='dashed', linewidth=2, label=f'2nd Pct (vmin): {vmin:.3f}')
        plt.axvline(vmax, color='red', linestyle='dashed', linewidth=2, label=f'98th Pct (vmax): {vmax:.3f}')
        plt.axvline(median, color='dodgerblue', linestyle='dotted', linewidth=2, label=f'Median: {median:.3f}')

        plt.title(f'Frequency Distribution of Mapped {element_name} Values', fontsize=16)
        plt.xlabel(f'{element_name} Concentration Value', fontsize=12)
        plt.ylabel('Frequency (Number of Pixels)', fontsize=12)
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
            
        file_path = os.path.join(save_dir, f"{element_name}_distribution_plot.png")
        plt.savefig(file_path, dpi=300)
        plt.show()

# =========================================================================
# WRAPPER FUNCTION TO PROCESS THE PARQUET FILE
# =========================================================================

def map_parquet_and_calculate_coverage(parquet_file, element="conc_Si", grid_size=(360, 720), num_partitions=20):
    print(f"Loading data from {parquet_file} for {element}...")
    df = pd.read_parquet(parquet_file)
    df = df[df['success'] == True]

    abundance_map = GaussianArray(grid_size=grid_size)
    
    block_lat =[-90, 90, 90, -90]
    block_lon =[-180, -180, 180, 180]
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Mapping Footprints"):
        fp = row['footprint']
        
        if isinstance(fp, str):
            fp = ast.literal_eval(fp)
        
        img_lat = [point[0] for point in fp]  
        img_lon = [point[1] for point in fp]  

        img_lon = [(lon + 180) % 360 - 180 for lon in img_lon]

        if max(img_lon) - min(img_lon) > 180:
            continue

        abundance = row[element]
        
        abundance_map.add_gaussian_box(img_lat, img_lon, block_lat, block_lon, abundance)

    coverage_fraction = abundance_map.check_coverage()
    covered_pct = coverage_fraction * 100
    empty_pct = 100 - covered_pct
    
    print("\n" + "="*40)
    print(" 🌕 COVERAGE STATISTICS 🌕")
    print("="*40)
    print(f"Target Element (ratio) : {element}")
    print(f"Total Area Covered : {covered_pct:.2f} %")
    print(f"Total Area Empty   : {empty_pct:.2f} %")
    print(f"Grid Resolution    : {grid_size[0]}x{grid_size[1]} pixels")
    print("="*40 + "\n")
    
    print("Generating value frequency distribution plot...")
    abundance_map.plot_value_distribution(element_name=element)

    print("Exporting and displaying abundance map...")
    abundance_map.export_map_png(f"{element}_abundance_map.png", plt_show=True, cmap='inferno', save_dir='./')
    
    print("Exporting and displaying footprint overlap heatmap...")
    abundance_map.visualize_counts()
    
    return abundance_map

def _best_partition_shape(num_partitions):
    rows = int(np.sqrt(num_partitions))
    while rows > 1 and num_partitions % rows != 0:
        rows -= 1
    cols = num_partitions // rows
    return rows, cols

def print_patchwise_average_concentration(abundance_map, element_name, num_partitions=20):
    arr = abundance_map.arr[:, :, 0]
    coverage = abundance_map.arr[:, :, 1] > 0

    part_rows, part_cols = _best_partition_shape(num_partitions)
    row_splits = np.array_split(np.arange(arr.shape[0]), part_rows)
    col_splits = np.array_split(np.arange(arr.shape[1]), part_cols)

    rows = []
    patch_id = 1
    for r_idx, r_inds in enumerate(row_splits):
        for c_idx, c_inds in enumerate(col_splits):
            patch_vals = arr[np.ix_(r_inds, c_inds)]
            patch_cov = coverage[np.ix_(r_inds, c_inds)]

            covered_pixels = int(np.sum(patch_cov))
            total_pixels = int(patch_vals.size)

            if covered_pixels > 0:
                avg_conc = float(np.mean(patch_vals[patch_cov]))
            else:
                avg_conc = np.nan

            rows.append({
                "patch_id": patch_id,
                "grid_pos": f"R{r_idx + 1}C{c_idx + 1}",
                "row_range": f"{r_inds[0]}-{r_inds[-1]}",
                "col_range": f"{c_inds[0]}-{c_inds[-1]}",
                "covered_pixels": covered_pixels,
                "total_pixels": total_pixels,
                f"avg_{element_name}": avg_conc
            })
            patch_id += 1

    patch_df = pd.DataFrame(rows)

    print("\n" + "=" * 72)
    print(f"Patch-wise Average Concentration ({num_partitions} partitions)")
    print("=" * 72)
    print(patch_df.to_string(index=False, float_format=lambda x: f"{x:.6f}"))
    print("=" * 72 + "\n")

    return patch_df

if __name__ == "__main__":
    map_obj = map_parquet_and_calculate_coverage(
        parquet_file='processed/Al_50.parquet', 
        element='conc_Al', 
        grid_size=(360, 720),
        num_partitions=20
    )