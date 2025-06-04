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
        print(f"  File 1: [{np.min(esp1):.6e}, {np.max(esp1):.6e}]")
        print(f"  File 2: [{np.min(esp2):.6e}, {np.max(esp2):.6e}]")
    
    return metrics


if __name__ == "__main__":
    # Get paths to the ESP files
    resp_esp = "../resp/1_default_grid_esp.dat"
    mpfit_esp = "1_default_grid_esp.dat"
    
    # Check if files exist
    if not os.path.exists(resp_esp):
        print(f"Error: RESP ESP file not found at {resp_esp}")
        sys.exit(1)
    
    if not os.path.exists(mpfit_esp):
        print(f"Error: MPFIT ESP file not found at {mpfit_esp}")
        print("Please run generate_esp.py first to create the ESP from MPFIT charges")
        sys.exit(1)
    
    # Compare the ESP files
    print(f"Comparing ESP files:")
    print(f"  RESP:  {resp_esp}")
    print(f"  MPFIT: {mpfit_esp}")
    
    metrics = compare_grid_esp(resp_esp, mpfit_esp, verbose=True)