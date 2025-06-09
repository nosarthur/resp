import psi4
import resp
from resp.utils import grid_to_mol2

mol = psi4.geometry(
""" 
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

# Call for first stage fit
charges1 = resp.resp([mol], options)
print('Electrostatic Potential Charges')
print(charges1[0])
print('Restrained Electrostatic Potential Charges')
print(charges1[1])

# Change the value of the RESP parameter A
options['RESP_A'] = 0.001

# Add constraint for atoms fixed in second stage fit
constraint_charge = []
for i in range(4, 8):
    constraint_charge.append([charges1[1][i], [i+1]])
options['constraint_charge'] = constraint_charge
options['constraint_group'] = [[2, 3, 4]]
options['grid'] = ['1_%s_grid.dat' %mol.name()]
options['esp'] = ['1_%s_grid_esp.dat' %mol.name()]

# Call for second stage fit
charges2 = resp.resp([mol], options)

# Get RESP charges
print("\nStage Two:\n")
print('RESP Charges')
print(charges2[1])

# Convert to MOL2 format for visualization
print("\nConverting RESP ESP to MOL2 format...")
grid_to_mol2(grid_file='1_default_grid.dat',
             esp_file='1_default_grid_esp.dat',
             output_mol2='resp_esp.mol2',
             normalize=False)

# Also create a PDB file of the molecule
from resp.utils import molecule_to_pdb
molecule_to_pdb(mol, filename='pyrazine.pdb', res_name='PYZ')
