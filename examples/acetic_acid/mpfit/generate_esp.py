import psi4
import resp
from resp.utils import compare_grid_esp, create_difference_esp, grid_to_mol2

mol = psi4.geometry(""" 
 O  -5.53102237   2.9394029    1.97762624
 O  -2.40708886   0.9081773    0.04905909
 H  -0.50304887   1.3757920    4.63880908
 H   0.00000000   4.2510428    3.09072238
 H  -2.38948215   4.0191677    5.40195041
 H  -3.74771411   0.4714793   -1.21557512
 C  -1.41859782   2.9850050    3.84405190
 C  -3.32263781   2.2862881    1.91221413
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
charges = [-0.53839, -0.50853, 0.09919, 0.09918, 0.08675, 0.36501, -0.33301, 0.72979]

# Electrostatic Potential Charges
# charges = [-0.64961348, -0.64445712, 0.13718876, 0.13746974, 0.11109979, 0.44304754, -0.40171205, 0.86697681]

# RESP Charges
# charges = [-0.84523957, -0.07718493, -0.07718493, -0.07718493, 0.06968509, 0.44599342, -0.22813929, 0.78925513]

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
