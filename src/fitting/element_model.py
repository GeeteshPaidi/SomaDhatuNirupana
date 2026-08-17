import numpy as np
import xraylib
from scipy.stats import norm
import matplotlib.pyplot as plt

class ElementModel:
    """
    Model for X-ray fluorescence (XRF) emission lines with Gaussian detector response.
    
    All energies and standard deviations are in keV. The class uses xraylib for
    atomic data including emission line energies, radiative rates, fluorescence yields,
    and cross sections.
    
    Note: For realistic XRF modeling, consider:
    - Using CS_FluorLine for excitation cross sections (implemented in primary_intensity)
    - Including detector efficiency and solid angle factors
    - Modeling self-absorption (matrix effects)
    - Adding continuum background (Compton, Bremsstrahlung)
    """
    def __init__(self, element, conc, std_dev):
        """
        Initialize the ElementModel with the specified element, conc, and standard deviation.
        
        Args:
            element (str): The symbol of the element (e.g., 'Fe').
            conc (float): The concentration. Units depend on application:
                - If representing total counts from element, use counts
                - If representing concentration (wt%, ppm, etc.), use that value
                The Gaussian amplitude will be proportional to conc.
            std_dev (float): The standard deviation of the Gaussian in keV.
                This typically represents detector resolution. Common values:
                - SDD detectors: 0.120-0.150 keV
                - Si(Li) detectors: 0.150-0.200 keV
                Convert from FWHM if needed: std_dev = FWHM / 2.355
        """
        self.element = element
        self.conc = conc
        self.std_dev = std_dev
        self.Z = xraylib.SymbolToAtomicNumber(element)  # Atomic number
        
        # Get K and L emission lines energies and probabilities
        energy_dict = {}
        radrate = {}
        
        lines = []
        self.lines_map = {}
        ka_lines = [
            (xraylib.KA1_LINE, "ka1"),
        ]
        kb_lines = [(xraylib.KB1_LINE, "kb1")]
        
        # Add K lines to means and radrate
        peak_energies = {}
        peak_rad_rates = {}
        for line, label in ka_lines:
            try:
                energy = xraylib.LineEnergy(self.Z, line)
                prob = xraylib.RadRate(self.Z, line)
                if energy > 0 and prob > 0:
                    self.lines_map[label] = line        
                    peak_energies[f"{label}"] = energy
                    peak_rad_rates[f"{label}"] = prob
                    lines.append(label)
            except:
                continue
        
        energy_dict["ka"] = peak_energies
        radrate["ka"] = peak_rad_rates
        
        # Add KB lines
        peak_energies_kb = {}
        peak_rad_rates_kb = {}
        for line, label in kb_lines:
            try:
                energy = xraylib.LineEnergy(self.Z, line)
                prob = xraylib.RadRate(self.Z, line)
                if energy > 0 and prob > 0:
                    self.lines_map[label] = line
                    peak_energies_kb[label] = energy
                    peak_rad_rates_kb[label] = prob
                    lines.append(label)
            except:
                continue
        
        energy_dict["kb"] = peak_energies_kb
        radrate["kb"] = peak_rad_rates_kb
        
        # L lines (you can add more if needed)
        la_lines = [
            # (xraylib.LA1_LINE, "la1"),
        ]
        lb_lines = [  
            #  (xraylib.LB1_LINE, "lb1"),
            # (xraylib.LB2_LINE, "lb2")
            ]
        peak_energies_la = {}
        peak_rad_rates_la = {}
        for line, label in la_lines:
            try:
                energy = xraylib.LineEnergy(self.Z, line)
                prob = xraylib.RadRate(self.Z, line)
                if energy > 0 and prob > 0:
                    self.lines_map[label] = line 
       
                    peak_energies_la[label] = energy
                    peak_rad_rates_la[label] = prob
                    lines.append(label)

            except:
                continue
        
        energy_dict["la"] = peak_energies_la
        radrate["la"] = peak_rad_rates_la
        peak_energies_lb = {}
        peak_rad_rates_lb = {}
        for line, label in lb_lines:
            try:
                energy = xraylib.LineEnergy(self.Z, line)
                prob = xraylib.RadRate(self.Z, line)
                if energy > 0 and prob > 0:
                    self.lines_map[label] = line        
                    peak_energies_lb[label] = energy
                    peak_rad_rates_lb[label] = prob
                    lines.append(label)
            except:
                continue
        
        energy_dict["lb"] = peak_energies_lb
        radrate["lb"] = peak_rad_rates_lb
        
        # Calculate means only for non-empty dictionaries
        means = {}
        for key, levels in energy_dict.items():
            if levels:  # Only calculate mean if dictionary is not empty
                means[key] = np.mean(list(levels.values()))
            else:
                means[key] = 0.0
        
        for key in energy_dict.keys():
            energy_dict[key]["mean"] = means[key]

        for key, levels in radrate.items():
            if levels:  # Only calculate mean if dictionary is not empty
                means[key] = np.mean(list(levels.values()))
            else:
                means[key] = 0.0
        
        for key in radrate.keys():
            radrate[key]["mean"] = means[key]

        self.energy_dict = energy_dict
        
        self.radrates = radrate
        self.lines = lines
        self.line_div = {}
        # line_div={
        #     "ka":["ka1","ka2"],
        #     "kb":["kb1","kb2"],
        #     "la":["la1"],
        #     "lb":["lb1"]
        # }
        for line in lines:
            if line[:2] not in self.line_div.keys():
                self.line_div[line[:2]] = []
            self.line_div[line[:2]].append(line)
        
        # Store per-line-type std_dev, defaulting to the provided std_dev
        self.std_dev_dict = {}
        for line_type in energy_dict.keys():
            self.std_dev_dict[line_type] = std_dev if std_dev else 0.01
        # self.std_devs = np.array(std_dev * len(energy_dict))
        # self.radrate = np.array(list(radrate.values()))

    def calculate_mass_absorption_coefficient(self, energy=None, element=None, line=None):
        """
        Calculate the mass absorption coefficient for the element at the given energy.
        
        Args:
            energy (float, optional): The energy value in keV. If provided, uses this directly.
            element (str, optional): Element symbol (unused, kept for compatibility).
            line (str, optional): Line label (e.g., 'ka1'). If provided, uses line energy.
        
        Returns:
            float: The mass absorption coefficient in cm^2/g.
        """
        if energy is not None:
            # Ensure energy is a float
            energy = float(energy)
            if energy <= 0:
                raise ValueError(f"Energy must be positive, got {energy}")
            return xraylib.CS_Total(self.Z, energy)
        elif line is not None:
            if line not in self.lines:
                return 0.0
            if line[:2] not in self.energy_dict or line not in self.energy_dict[line[:2]]:
                return 0.0
            energy = self.energy_dict[line[:2]][line]
            # Ensure energy is a float
            energy = float(energy)
            if energy <= 0:
                return 0.0
            return xraylib.CS_Total(self.Z, energy)
        else:
            raise ValueError("Either energy or line must be provided")
    
    def calculate_elemental_const(self, line, ey, use_per_line=True):
        """
        Calculate the elemental constant for the given line.
        
        Args:
            line (str): The line label (e.g., 'ka1').
            ey (float): The incident energy value in keV (currently unused but kept for API compatibility).
            use_per_line (bool): If True, uses per-line radiative rate. If False, uses family mean.
        
        Returns:
            float: The elemental constant at the specified energy.
        """
        # Determine line type and shell
        if "ka" in line:
            ltype = "ka"
            shell = 0  # K shell
        elif "kb" in line:
            ltype = "kb"
            shell = 0  # K shell
        elif "la" in line:
            ltype = "la"
            shell = 1  # L shell
        elif "lb" in line:
            ltype = "lb"
            shell = 1  # L shell
        else:
            raise ValueError(f"Unsupported line type: {line}")
        
        # Get jump ratio and check for zero
        try:
            rk = xraylib.JumpFactor(self.Z, shell)
        except (ValueError, OverflowError):
            raise ValueError(f"Could not determine jump factor for {self.element} shell {shell}")
        
        if rk == 0 or rk < 0:
            raise ValueError(f"Invalid jump ratio (rk={rk}) for {self.element} {line}. "
                           f"Cannot compute elemental constant.")
        
        # Calculate jump factor (1 - 1/rk)
        jump_factor = 1 - 1 / rk
        
        # Get fluorescence yield
        try:
            fluor_yield = xraylib.FluorYield(self.Z, shell)
        except (ValueError, OverflowError) as e:
            print(f"Warning: {e} for {self.element} at {line}, returning 0")
            return 0
        
        # Get radiative rate - use per-line if available and requested
        if use_per_line and line in self.lines and line in self.energy_dict[ltype]:
            rad_rate = self.radrates[ltype].get(line, 0)
            if rad_rate == 0:
                # Fall back to mean if per-line not available
                rad_rate = self.radrates[ltype].get("mean", 0)
        else:
            rad_rate = self.radrates[ltype].get("mean", 0)
        
        if rad_rate == 0:
            raise ValueError(f"No radiative rate found for line {line}")
        
        c = jump_factor * fluor_yield * rad_rate
        return c

    def gaussian(self, x, mean, std_dev):
        """
        Calculate the Gaussian function value at the specified energy.
        
        Note: Multiplying by conc makes the area under the peak equal to conc.
        If conc represents total counts, this is appropriate. If conc is a 
        concentration (e.g., wt%), you may want to scale differently based on
        incident flux, detector efficiency, etc.
        
        Args:
            x (float or array): The energy value(s) in keV.
            mean (float): The mean energy in keV.
            std_dev (float): The standard deviation in keV (detector resolution).
        
        Returns:
            float or array: The Gaussian value(s) at the specified energy(s).
            Returns scalar if x is scalar, array if x is array.
        """
        if std_dev == 0 or std_dev < 0:
            # Return zeros of same shape as input
            if isinstance(x, np.ndarray):
                return np.zeros_like(x)
            return 0.0
        
        # Ensure mean and std_dev are in keV (they should be)
        # norm.pdf returns probability density, which is normalized to area = 1
        # Multiplying by conc makes the area under the peak = conc
        return self.conc * norm.pdf(x, mean, std_dev)
    def primary_intensity(self, line, incident_energy=None, use_fluor_line=True):
        """
        Calculate the primary X-ray fluorescence intensity for a given line.
        
        Args:
            line (str): The line label (e.g., 'ka1').
            incident_energy (float, optional): Incident beam energy in keV. 
                If None, uses a value above the absorption edge.
            use_fluor_line (bool): If True, uses CS_FluorLine for more accurate 
                fluorescent cross section. If False, uses CS_Total at emission energy.
        
        Returns:
            float: The primary intensity (proportional to counts).
        """
        if line not in self.lines:
            raise ValueError(f"Line '{line}' not found in available lines")
        
        line_type = line[:2]
        if line_type not in self.std_dev_dict:
            raise ValueError(f"Line type '{line_type}' not found")
        
        # Get per-line energy if available, otherwise use mean
        if line in self.energy_dict[line_type]:
            line_energy = self.energy_dict[line_type][line]
        else:
            line_energy = self.energy_dict[line_type].get("mean", 0)
            if line_energy == 0:
                raise ValueError(f"No energy found for line '{line}'")
        
        if use_fluor_line and line in self.lines_map:
            # Use CS_FluorLine which gives the fluorescent line cross section directly
            try:
                if incident_energy is None:
                    # Auto-determine incident energy (must be above absorption edge)
                    # For K lines, use K edge energy + 10%, for L use L3 edge + 10%
                    if "k" in line.lower():
                        edge_energy = xraylib.EdgeEnergy(self.Z, xraylib.K_SHELL)
                    else:
                        edge_energy = xraylib.EdgeEnergy(self.Z, xraylib.L3_SHELL)
                    incident_energy = edge_energy * 1.1
                
                line_id = self.lines_map[line]
                cs_fluor = xraylib.CS_FluorLine(self.Z, line_id, incident_energy)
                intensity = cs_fluor * self.conc
            except (ValueError, OverflowError) as e:
                # Fall back to CS_Total if CS_FluorLine fails
                mass_absorption_coeff = self.calculate_mass_absorption_coefficient(energy=line_energy)
                intensity = mass_absorption_coeff * self.conc
        else:
            # Fallback: use mass absorption coefficient at emission energy
            # Note: This is simplified and not physically correct for fluorescence
            mass_absorption_coeff = self.calculate_mass_absorption_coefficient(energy=line_energy)
            intensity = mass_absorption_coeff * self.conc
        
        return intensity
    def jump_ratio_factor(self, shell):
        """
        Calculate the jump ratio rk for the given energy using xraylib.
        
        Args:
            energy (float): The energy value in keV.
        
        Returns:
            float: The jump ratio rk at the specified energy.
        """
        rk = xraylib.JumpFactor(self.Z, shell)  # Assuming JumpFactor returns rk
        if rk == 0:
            raise ValueError("rk must not be zero to avoid division by zero.")
        return 1 - 1 / rk
    
    def weighted_gaussian(self, energy, line=None, use_per_line=False):
        """
        Calculate the weighted Gaussian function value at the specified energy.
        
        Args:
            energy (float or array): The energy value(s) in keV.
            line (str, optional): Line label (e.g., 'ka1') or line type (e.g., 'ka'). 
                If None, uses the first available line.
            use_per_line (bool): If True and line is a specific line (e.g., 'ka1'), 
                uses that line's exact energy. If False, uses line type mean energy.
        
        Returns:
            float or array: The weighted Gaussian value(s) at the specified energy.
            Returns scalar if energy is scalar, array if energy is array.
        """
        if line is None:
            # Use the first available line
            if not self.lines:
                raise ValueError("No lines available for this element")
            line = self.lines[0]
            use_per_line = True  # Use specific line energy if we have the label
        
        # Determine if line is a specific line (e.g., 'ka1') or just type (e.g., 'ka')
        if line in self.lines:
            # Specific line label
            line_type = line[:2]
            if use_per_line:
                if line in self.energy_dict[line_type]:
                    line_energy = self.energy_dict[line_type][line]
                else:
                    line_energy = self.energy_dict[line_type].get("mean", 0)
            else:
                line_energy = self.energy_dict[line_type].get("mean", 0)
        else:
            # Line type only
            line_type = line[:2] if len(line) >= 2 else line
            line_energy = self.energy_dict[line_type].get("mean", 0)
        
        if line_type not in self.energy_dict:
            raise ValueError(f"Line type '{line_type}' not found in energy_dict")
        
        if line_energy == 0:
            raise ValueError(f"No energy found for line/type '{line}'")
        
        std_dev = self.std_dev_dict.get(line_type, self.std_dev)
        
        # Ensure std_dev is in keV (it should be, but document the assumption)
        # std_dev is expected to be in keV throughout
        
        mass_absorption_coeff = self.calculate_mass_absorption_coefficient(energy=line_energy)
        
        # Handle both scalar and array inputs
        gaussian_vals = self.gaussian(energy, line_energy, std_dev)
        
        # Multiply by mass absorption coefficient
        # If energy is scalar, gaussian returns scalar; if array, returns array
        result = gaussian_vals * mass_absorption_coeff
        
        return result
    
    def compute_per_line_spectrum(self, energies, incident_energy=None):
        """
        Compute the full spectrum using per-line values (not family means).
        This gives more accurate representation with separate Kα1, Kα2, etc. peaks.
        
        Args:
            energies (array): Energy values in keV to compute spectrum at.
            incident_energy (float, optional): Incident beam energy in keV for CS_FluorLine.
                If None, auto-determines appropriate value.
        
        Returns:
            array: Total intensity at each energy value.
        """
        spectrum = np.zeros_like(energies, dtype=float)
        
        # Auto-determine incident energy if needed
        if incident_energy is None:
            # Use K edge + 10% for K lines, L3 edge + 10% for L lines
            try:
                k_edge = xraylib.EdgeEnergy(self.Z, xraylib.K_SHELL)
                l_edge = xraylib.EdgeEnergy(self.Z, xraylib.L3_SHELL)
                incident_energy = max(k_edge, l_edge) * 1.1
            except (ValueError, OverflowError):
                incident_energy = 20.0  # Default fallback
        
        # Process each individual line
        for line_label in self.lines:
            if line_label not in self.lines_map:
                continue
            
            line_type = line_label[:2]
            if line_type not in self.energy_dict:
                continue
            
            # Get per-line energy
            if line_label in self.energy_dict[line_type]:
                line_energy = self.energy_dict[line_type][line_label]
            else:
                continue  # Skip if line energy not available
            
            # Get per-line radiative rate
            if line_label in self.radrates[line_type]:
                rad_rate = self.radrates[line_type][line_label]
            else:
                rad_rate = self.radrates[line_type].get("mean", 0)
            
            if rad_rate == 0:
                continue
            
            # Get std_dev for this line type
            std_dev = self.std_dev_dict.get(line_type, self.std_dev)
            
            try:
                primary_int = self.primary_intensity(line_label, incident_energy=incident_energy, use_fluor_line=True)
            except (ValueError, ZeroDivisionError):
                # Fallback: calculate using elemental constant
                try:
                    elem_const = self.calculate_elemental_const(line_label, incident_energy, use_per_line=True)
                    primary_int = elem_const * self.conc
                except (ValueError, ZeroDivisionError):
                    continue
            
            gaussian_contrib = norm.pdf(energies, line_energy, std_dev)
            spectrum += gaussian_contrib * primary_int
        
        return spectrum
    
    def plot(self, energy_range=None, line_types=None, show_weighted=True, show_individual_lines=True, figsize=(12, 8), use_per_line=False):
        """
        Plot the Gaussian distributions for emission lines and optionally the weighted Gaussian.
        
        Args:
            energy_range (tuple, optional): (min_energy, max_energy) in keV. If None, auto-calculated.
            line_types (list, optional): List of line types to plot (e.g., ['ka', 'kb']). If None, plots all available.
            show_weighted (bool): Whether to show the weighted Gaussian curve.
            show_individual_lines (bool): Whether to show individual emission line positions as vertical lines.
            figsize (tuple): Figure size (width, height) in inches.
            use_per_line (bool): If True, uses per-line spectrum computation (more accurate).
        
        Returns:
            matplotlib.figure.Figure: The figure object.
        """
        # Determine which line types to plot
        if line_types is None:
            line_types = [lt for lt in self.energy_dict.keys() if self.energy_dict[lt] and self.energy_dict[lt].get("mean", 0) > 0]
        
        if not line_types:
            raise ValueError("No valid line types available for plotting")
        
        # Determine energy range if not provided
        if energy_range is None:
            min_energy = float('inf')
            max_energy = 0
            for lt in line_types:
                if lt in self.energy_dict and self.energy_dict[lt].get("mean", 0) > 0:
                    mean_energy = self.energy_dict[lt]["mean"]
                    std_dev = self.std_dev_dict.get(lt, self.std_dev)
                    min_energy = min(min_energy, mean_energy - 5 * std_dev)
                    max_energy = max(max_energy, mean_energy + 5 * std_dev)
            
            if min_energy == float('inf'):
                raise ValueError("Could not determine energy range")
            
            # Add some padding
            padding = (max_energy - min_energy) * 0.1
            energy_range = (max(0, min_energy - padding), max_energy + padding)
        
        energy_min, energy_max = energy_range
        energies = np.linspace(energy_min, energy_max, 1000)
        
        fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)
        
        # Plot 1: Individual Gaussian distributions and emission lines
        ax1 = axes[0]
        
        colors = {'ka': 'blue', 'kb': 'green', 'la': 'red', 'lb': 'orange'}
        
        for lt in line_types:
            if lt not in self.energy_dict or not self.energy_dict[lt]:
                continue
            
            mean_energy = self.energy_dict[lt].get("mean", 0)
            if mean_energy == 0:
                continue
            
            std_dev = self.std_dev_dict.get(lt, self.std_dev)
            color = colors.get(lt, 'gray')
            
            # Plot Gaussian distribution
            gaussian_vals = self.gaussian(energies, mean_energy, std_dev)
            ax1.plot(energies, gaussian_vals, label=f'{lt.upper()} Gaussian', 
                    color=color, linestyle='-', linewidth=2)
            
            # Plot individual emission line positions
            if show_individual_lines:
                for line_name, line_energy in self.energy_dict[lt].items():
                    if line_name != "mean" and line_energy > 0:
                        # Add a vertical line at the emission line energy
                        max_val = np.max(gaussian_vals) if len(gaussian_vals) > 0 else 1
                        ax1.axvline(x=line_energy, color=color, linestyle='--', 
                                   alpha=0.5, linewidth=1)
                        # Add label
                        if energy_min <= line_energy <= energy_max:
                            ax1.text(line_energy, max_val * 0.9, line_name, 
                                   rotation=90, verticalalignment='bottom',
                                   horizontalalignment='right', fontsize=8, color=color)
        
        ax1.set_ylabel('Intensity', fontsize=12)
        ax1.set_title(f'{self.element} Emission Line Spectra (Concentration: {self.conc})', 
                     fontsize=14, fontweight='bold')
        ax1.legend(loc='best')
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Weighted Gaussian (if requested)
        ax2 = axes[1]
        
        if show_weighted:
            if use_per_line:
                # Use per-line spectrum computation for accuracy
                spectrum_vals = self.compute_per_line_spectrum(energies)
                ax2.plot(energies, spectrum_vals, label='Per-Line Spectrum', 
                        color='black', linestyle='-', linewidth=2)
            else:
                # Use family means (less accurate but faster)
                for lt in line_types:
                    if lt not in self.energy_dict or not self.energy_dict[lt]:
                        continue
                    
                    mean_energy = self.energy_dict[lt].get("mean", 0)
                    if mean_energy == 0:
                        continue
                    
                    color = colors.get(lt, 'gray')
                    
                    # Calculate weighted Gaussian
                    weighted_vals = np.array([self.weighted_gaussian(e, line=lt) for e in energies])
                    ax2.plot(energies, weighted_vals, label=f'{lt.upper()} Weighted', 
                            color=color, linestyle='-', linewidth=2)
        
        ax2.set_xlabel('Energy (keV)', fontsize=12)
        ax2.set_ylabel('Weighted Intensity', fontsize=12)
        ax2.set_title('Weighted Gaussian Distributions', fontsize=14, fontweight='bold')
        ax2.legend(loc='best')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig
    
    def plot_mass_absorption(self, energy_range=None, figsize=(10, 6)):
        """
        Plot the mass absorption coefficient as a function of energy.
        
        Args:
            energy_range (tuple, optional): (min_energy, max_energy) in keV. Defaults to (1, 20).
            figsize (tuple): Figure size (width, height) in inches.
        
        Returns:
            matplotlib.figure.Figure: The figure object.
        """
        if energy_range is None:
            energy_range = (1, 20)
        
        energy_min, energy_max = energy_range
        energies = np.linspace(energy_min, energy_max, 1000)
        
        mac_values = np.array([self.calculate_mass_absorption_coefficient(energy=e) for e in energies])
        
        fig, ax = plt.subplots(figsize=figsize)
        ax.plot(energies, mac_values, 'b-', linewidth=2)
        ax.set_xlabel('Energy (keV)', fontsize=12)
        ax.set_ylabel('Mass Absorption Coefficient (cm²/g)', fontsize=12)
        ax.set_title(f'Mass Absorption Coefficient for {self.element}', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_yscale('log')
        
        plt.tight_layout()
        return fig
        
# Example usage:
if __name__ == "__main__":
    element = 'Fe'  # Iron
    energy = 6.4  
    conc = 1.0
    std_dev = 0.1 
    
    model = ElementModel(element, conc, std_dev)
    result = model.weighted_gaussian(energy)
    print(f"Weighted Gaussian value for {element} at {energy} keV: {result}")
    
    # Plot the spectra
    fig = model.plot()
    fig.savefig('element_spectrum.png')

    # Plot mass absorption coefficient
    fig2 = model.plot_mass_absorption()
    fig2.savefig('mass_absorption_coefficient.png')