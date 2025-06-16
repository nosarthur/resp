"""
Utility functions for RESP/MPFIT comparison examples.
"""
import numpy as np
import os
import sys

# Add the resp module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from resp.vdw_surface import vdw_surface

# Physical constants
bohr_to_angstrom = 0.529177249

def charges_to_esp(molecule, charges, options=None, psi4_esp_file="1_default_grid_esp.dat", verbose=True):
    """Generate ESP grid files from given charges.
    
    Parameters
    ----------
    molecule : psi4.Molecule
        Molecule instance with geometry
    charges : array_like
        Atomic charges to use for ESP calculation
    options : dict, optional
        Options for grid generation (uses same options as resp())
    
    Returns
    -------
    None
       Prints metrics found in `compare_grid_esp` 
    """
    
    # Check options - same as resp()
    options = {k.upper(): v for k, v in sorted(options.items())}    

    # Get coordinates in angstrom - same as resp()
    coordinates = molecule.geometry()
    coordinates = coordinates.np.astype('float')*bohr_to_angstrom
    
    # Read grid points
    points = np.loadtxt(options['GRID'][0])
    if 'Bohr' in str(molecule.units()):
        points *= bohr_to_angstrom

    # Calculate ESP values - following same pattern as invr calculation in resp()
    esp_values = np.zeros(len(points))
    for i in range(len(points)):
        for j in range(len(coordinates)):
            distance = np.linalg.norm(points[i] - coordinates[j])
            # Convert to atomic units: ESP (a.u.) = charge * (1/r_angstrom) * bohr_to_angstrom
            esp_values[i] += charges[j] / distance * bohr_to_angstrom
    
    esp_psi4 = np.loadtxt(psi4_esp_file)
    compare_grid_esp(esp_psi4, esp_values, verbose=verbose) 
    return


def compare_grid_esp(esp1, esp2, verbose=True):
    """Compare two grid_esp.dat files and calculate metrics.
    
    Parameters
    ----------
    file1 : str
        Path to first grid_esp.dat file
    file2 : str
        Path to second grid_esp.dat file
    verbose : bool
        Print detailed comparison metrics
    
    Returns
    -------
    dict
        Dictionary containing comparison metrics
    """
    
    # Check if arrays have same shape
    if esp1.shape != esp2.shape:
        raise ValueError(f"ESP arrays have different shapes: {esp1.shape} vs {esp2.shape}")
    
    # Calculate metrics
    metrics = {}
    
    # L2 norm - this is analogous to ||A @ q - b|| in RESP fitting
    diff = esp1 - esp2
    metrics['l2_norm'] = np.sqrt(np.sum(diff**2))
    
    # RMSE (L2 norm normalized by number of points)
    metrics['rmse'] = np.sqrt(np.mean(diff**2))
    
    # MAE (L1 norm normalized by number of points)
    metrics['mae'] = np.mean(np.abs(diff))
    
    # Max absolute error
    metrics['max_error'] = np.max(np.abs(diff))
    
    # Pearson correlation coefficient
    metrics['correlation'] = np.corrcoef(esp1, esp2)[0, 1]
    
    # R-squared
    ss_res = np.sum(diff**2)
    ss_tot = np.sum((esp1 - np.mean(esp1))**2)
    metrics['r_squared'] = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    if verbose:
        print(f"\nESP Comparison Results:")
        print(f"=======================")
        print(f"Number of grid points: {len(esp1)}")
        print(f"\nError Metrics:")
        print(f"  L2 norm:               {metrics['l2_norm']:.6e}")
        print(f"  RMSE:                  {metrics['rmse']:.6e}")
        print(f"  MAE:                   {metrics['mae']:.6e}")
        print(f"  Max absolute error:    {metrics['max_error']:.6e}")
        print(f"\nCorrelation Metrics:")
        print(f"  Pearson correlation:   {metrics['correlation']:.6f}")
        print(f"  R-squared:             {metrics['r_squared']:.6f}")
    return metrics
