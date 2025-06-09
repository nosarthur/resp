import psi4
import resp
from resp.utils import compare_grid_esp, create_difference_esp, grid_to_mol2

mol = psi4.geometry(""" 
 N  -1.33084758   0.98032658   0.42627037
 N   3.95377333   0.98032658   0.42627037
 H   3.56713674   4.35218569  -0.94230973
 H  -0.94421098   4.35218569  -0.94230973
 H  -0.94421098  -2.39153253   1.79485047
 H   3.56713674  -2.39153253   1.79485047
 C   2.62292576   2.94097974  -0.42627037
 C   0.00000000   2.94097974  -0.42627037
 C   0.00000000  -0.98032658   1.27881112
 C   2.62292576  -0.98032658   1.27881112
units bohr
""")

mol.update_geometry()

options = {'VDW_SCALE_FACTORS'  : [1.4, 1.6, 1.8, 2.0],
           'VDW_POINT_DENSITY'  : 1.0,
           'RESP_A'             : 0.0005,
           'RESP_B'             : 0.1,
	   'METHOD_ESP'         : 'hf', 
	   'BASIS_ESP'          : 'aug-cc-pvdz', 
           }

# MPFIT-generated charges 
charges = [-0.43756, -0.43999, 0.06191, 0.02255, 0.06130, 0.00309, 0.14618, 0.19927, 0.14510, 0.23815]

# Electrostatic Potential Charges
#charges = [0.54253546, 0.47763582, 0.24555997, 0.01559197, 0.26587405, -0.00413454, -0.91607898, 0.14309417, -1.02731905, 0.25724113]

# RESP Charges
# charges = [-0.33184326, 0.31842386, 0.31842386, 0.31842386, 0.23053808, 0.03031718, -0.83229597, 0.06616556, 0.0636741, -0.18182728]

# Generate ESP grid files from charges
resp.charges_to_esp(mol, charges, options)

print(f"Generated grid files:")
print(f"  1_{mol.name()}_grid.dat")
print(f"  1_{mol.name()}_grid_esp.dat")

# Compare with RESP-generated ESP
print("\nComparing MPFIT-generated ESP with RESP-generated ESP...")
resp_esp = "../resp/1_default_grid_esp.dat"
mpfit_esp = "1_default_grid_esp.dat"

try:
    metrics = compare_grid_esp(resp_esp, mpfit_esp, verbose=True)
    
    # Create difference ESP file
    print("\nCreating difference ESP file...")
    create_difference_esp(resp_esp, mpfit_esp, output_file='difference_grid_esp.dat')
    
    # Convert to MOL2 format for visualization
    print("\nConverting to MOL2 format...")
    grid_to_mol2(grid_file='1_default_grid.dat', 
                 esp_file='1_default_grid_esp.dat',
                 output_mol2='mpfit_esp.mol2',
                 normalize=False)
    
    # Also create MOL2 for the difference
    print("\nCreating MOL2 for difference ESP...")
    grid_to_mol2(grid_file='1_default_grid.dat',
                 esp_file='difference_grid_esp.dat', 
                 output_mol2='difference_esp.mol2',
                 normalize=False)
    
except Exception as e:
    print(f"Error during comparison: {e}")
