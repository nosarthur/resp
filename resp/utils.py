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
