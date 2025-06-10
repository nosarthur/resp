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

def charges_to_esp(molecule, charges, options=None):
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
        Writes grid.dat and grid_esp.dat files
    """
    if options is None:
        options = {}
    
    # Check options - same as resp()
    options = {k.upper(): v for k, v in sorted(options.items())}
    
    # VDW surface options - same defaults as resp()
    if 'VDW_SCALE_FACTORS' not in options:
        options['VDW_SCALE_FACTORS'] = [1.4, 1.6, 1.8, 2.0]
    if 'VDW_POINT_DENSITY' not in options:
        options['VDW_POINT_DENSITY'] = 1.0
    if 'VDW_RADII' not in options:
        options['VDW_RADII'] = {}
    
    # Process VDW radii dict - same as resp()
    radii = {}
    for i in options['VDW_RADII']:
        radii[i.upper()] = options['VDW_RADII'][i]
    options['VDW_RADII'] = radii
    
    # Get symbols - same as resp()
    symbols = []
    for i in range(molecule.natom()):
        symbols.append(molecule.symbol(i))
    
    # Get coordinates in angstrom - same as resp()
    coordinates = molecule.geometry()
    coordinates = coordinates.np.astype('float')*bohr_to_angstrom
    
    # Generate grid points - same as resp()
    points = []
    for scale_factor in options['VDW_SCALE_FACTORS']:
        shell, radii = vdw_surface(coordinates, symbols, scale_factor,
                            options['VDW_POINT_DENSITY'], options['VDW_RADII'])
        points.append(shell)
    points = np.concatenate(points)
    
    # Handle grid filename
    grid_filename = 'grid.dat'
    if 'GRID' in options and options['GRID']:
        grid_filename = options['GRID'][0]
    
    # Save grid file - same format as resp()
    if 'Bohr' in str(molecule.units()):
        points /= bohr_to_angstrom
        np.savetxt(grid_filename, points, fmt='%15.10f')
        points *= bohr_to_angstrom
    else:
        np.savetxt(grid_filename, points, fmt='%15.10f')
    
    # Calculate ESP values - following same pattern as invr calculation in resp()
    esp_values = np.zeros(len(points))
    for i in range(len(points)):
        for j in range(len(coordinates)):
            distance = np.linalg.norm(points[i] - coordinates[j])
            # Convert to atomic units: ESP (a.u.) = charge * (1/r_angstrom) * bohr_to_angstrom
            esp_values[i] += charges[j] / distance * bohr_to_angstrom
    
    # Handle esp filename
    esp_filename = 'grid_esp.dat'
    if 'ESP' in options and options['ESP']:
        esp_filename = options['ESP'][0]
    
    # Save ESP file
    np.savetxt(esp_filename, esp_values, fmt='%15.10f')
    
    # Move files to numbered format if default names were used
    if grid_filename == 'grid.dat':
        os.system("mv grid.dat 1_%s_grid.dat" % molecule.name())
    if esp_filename == 'grid_esp.dat':
        os.system("mv grid_esp.dat 1_%s_grid_esp.dat" % molecule.name())
    
    return


def compare_grid_esp(file1, file2, verbose=True):
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
    # Load ESP values
    esp1 = np.loadtxt(file1)
    esp2 = np.loadtxt(file2)
    
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
    
    # L1 norm (sum of absolute differences)
    metrics['l1_norm'] = np.sum(np.abs(diff))
    
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
        print(f"  L1 norm:               {metrics['l1_norm']:.6e}")
        print(f"  MAE:                   {metrics['mae']:.6e}")
        print(f"  Max absolute error:    {metrics['max_error']:.6e}")
        print(f"\nCorrelation Metrics:")
        print(f"  Pearson correlation:   {metrics['correlation']:.6f}")
        print(f"  R-squared:             {metrics['r_squared']:.6f}")
        
        # ESP value ranges
        print(f"\nESP Value Ranges:")
        print(f"  File 1 {file1}: [{np.min(esp1):.6e}, {np.max(esp1):.6e}]")
        print(f"  File 2 {file2}: [{np.min(esp2):.6e}, {np.max(esp2):.6e}]")
    
    return metrics


def create_difference_esp(file1, file2, output_file='difference_grid_esp.dat'):
    """Create a grid_esp.dat file containing the difference between two ESP files.
    
    Parameters
    ----------
    file1 : str
        Path to first grid_esp.dat file
    file2 : str
        Path to second grid_esp.dat file
    output_file : str
        Path to output difference file (default: 'difference_grid_esp.dat')
    
    Returns
    -------
    np.ndarray
        Array of difference values
    """
    # Load ESP values
    esp1 = np.loadtxt(file1)
    esp2 = np.loadtxt(file2)
    
    # Check if arrays have same shape
    if esp1.shape != esp2.shape:
        raise ValueError(f"ESP arrays have different shapes: {esp1.shape} vs {esp2.shape}")
    
    # Calculate difference
    difference = esp1 - esp2
    print(f"Difference between {file1} - {file2}")
    
    # Save difference to file
    np.savetxt(output_file, difference, fmt='%14.10e')
    
    print(f"Created difference ESP file: {output_file}")
    print(f"Difference range: [{np.min(difference):.6e}, {np.max(difference):.6e}]")
    
    return difference
