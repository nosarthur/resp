import psi4
import resp

mol = psi4.geometry(
""" 
 N   6.27981210  -4.25303372   0.19313865
 N   6.06870777  -2.28125697   1.80544317
 H   3.66871156  -5.64628941  -1.96684001
 H   0.69034521  -2.37545694  -0.42864236
 H   3.35917442   0.17650601   3.09265017
 H   7.35266782  -1.99461311   2.90680702
 C   3.77952586  -1.21508453   1.76858172
 C   2.38242998  -2.55505538   0.04607681
 C   4.01771468  -4.38410490  -0.88761863
""")

# MPFIT-generated charges 
q = [-0.39242, -0.02616, 0.10538, 0.13998, 0.12150, 0.25454, 0.00422, -0.30743, 0.10040]


mol.update_geometry()

options = {'VDW_SCALE_FACTORS'  : [1.4, 1.6, 1.8, 2.0],
           'VDW_POINT_DENSITY'  : 1.0,
           'RESP_A'             : 0.0005,
           'RESP_B'             : 0.1,
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
