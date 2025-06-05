import numpy as np
import os
import sys


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


def grid_to_mol2(grid_file='1_default_grid.dat', 
                 esp_file='1_default_grid_esp.dat', 
                 output_mol2='sample.mol2',
                 normalize=True,
                 save_info=True):
    """Convert grid.dat and grid_esp.dat files to MOL2 format.
    
    This is a Python implementation of the esp.sh script that creates a MOL2 file
    where each grid point is represented as a carbon atom with ESP values.
    
    Parameters
    ----------
    grid_file : str
        Path to grid.dat file containing coordinates
    esp_file : str
        Path to grid_esp.dat file containing ESP values
    output_mol2 : str
        Path to output MOL2 file
    normalize : bool
        Whether to normalize ESP values (default: True)
    save_info : bool
        Whether to save min/max info to Infos.data (default: True)
    
    Returns
    -------
    dict
        Dictionary containing 'min', 'max', and 'normalized_values' (if normalize=True)
    """
    # Load grid coordinates and ESP values
    grid_coords = np.loadtxt(grid_file)
    esp_values = np.loadtxt(esp_file)
    
    # Ensure we have the right shapes
    if len(grid_coords) != len(esp_values):
        raise ValueError(f"Grid and ESP files have different lengths: {len(grid_coords)} vs {len(esp_values)}")
    
    num_points = len(esp_values)
    
    # Process ESP values
    esp_min = np.min(esp_values)
    esp_max = np.max(esp_values)
    
    result_info = {
        'min': esp_min,
        'max': esp_max
    }
    
    # Normalize if requested
    if normalize:
        normalized_values = np.zeros_like(esp_values)
        for i, value in enumerate(esp_values):
            if value < 0:
                # For negative values, normalize using the minimum
                normalized_values[i] = -value / esp_min if esp_min != 0 else 0
            else:
                # For positive values, normalize using the maximum  
                normalized_values[i] = value / esp_max if esp_max != 0 else 0
        
        result_info['normalized_values'] = normalized_values
        charge_values = normalized_values
        
        # Save normalized values
        np.savetxt('normalized_value', normalized_values, fmt='%9.6f')
    else:
        charge_values = esp_values
        # Save unnormalized values
        np.savetxt('unnormalized_value', esp_values, fmt='%9.6f')
    
    # Save info file if requested
    if save_info:
        with open('Infos.data', 'w') as f:
            f.write(f"max={esp_max}\n")
            f.write(f"min={esp_min}\n")
    
    # Create MOL2 file
    with open(output_mol2, 'w') as f:
        # Write header
        f.write("@<TRIPOS>MOLECULE\n")
        f.write("*****\n")
        f.write(f" {num_points} 0 0 0 0\n")
        f.write("SMALL\n")
        f.write("GASTEIGER\n")
        f.write("\n")
        f.write("@<TRIPOS>ATOM\n")
        
        # Write atoms (grid points)
        for i in range(num_points):
            atom_id = i + 1
            x, y, z = grid_coords[i]
            charge = charge_values[i]
            
            # Format matches the shell script output with proper spacing
            f.write(f"{atom_id:7d}     C\t{x:13.10f} {y:13.10f} {z:13.10f}\t  C3      1    UNL1\t{charge:10.6f}\n")
    
    print(f"Created MOL2 file: {output_mol2}")
    print(f"Number of grid points: {num_points}")
    print(f"ESP range: [{esp_min:.6e}, {esp_max:.6e}]")
    if normalize:
        print(f"Normalized range: [{np.min(normalized_values):.6f}, {np.max(normalized_values):.6f}]")
    
    return result_info


def molecule_to_pdb(molecule, filename='molecule.pdb', res_name='MOL'):
    """Write a Psi4 molecule to a PDB file.
    
    Parameters
    ----------
    molecule : psi4.core.Molecule
        The Psi4 molecule object
    filename : str
        Output PDB filename (default: 'molecule.pdb')
    res_name : str
        Residue name for the PDB file (default: 'MOL')
    
    Returns
    -------
    str
        Path to the created PDB file
    """
    # Check if coordinates need to be converted from Bohr to Angstroms
    bohr_to_angstrom = 0.52917721092
    units_str = str(molecule.units())
    convert_units = 'Bohr' in units_str or 'bohr' in units_str
    
    natoms = molecule.natom()
    
    pdb_lines = []
    
    # Write ATOM records
    for i in range(natoms):
        # Get atom information
        symbol = molecule.symbol(i)
        x, y, z = molecule.x(i), molecule.y(i), molecule.z(i)
        
        # Convert from Bohr to Angstroms if needed
        if convert_units:
            x *= bohr_to_angstrom
            y *= bohr_to_angstrom
            z *= bohr_to_angstrom
        
        # Atom name - use symbol + atom number
        atom_name = f"{symbol}{i+1}"
        
        # Format PDB ATOM line (proper PDB format)
        # ATOM serial atom_name res_name chain res_seq x y z occupancy temp_factor element
        pdb_line = (
            f"ATOM  {i+1:5d} {atom_name:<4} {res_name} A{1:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {symbol}"
        )
        pdb_lines.append(pdb_line)
    
    # Write END record
    pdb_lines.append("END")
    
    # Write to file
    with open(filename, 'w') as f:
        f.write('\n'.join(pdb_lines))
    
    print(f"Created PDB file: {filename}")
    return filename
