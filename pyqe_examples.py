import os
os.system('python pyqe.py')
import matplotlib.pyplot as plt
import numpy as np

from pyqe import *


file = File()
plot = Plot()



names = {}
names['files_folder'] = 'Examples/BN/bands_k/USPP'
names['scf_in_name'] = 'bn.scf.in'
names['scf_out_name'] = 'bn.scf.out'
names['bands_in_name'] = 'bn.bands.in'
names['bands_out_name'] = 'bn.bands.out'
names['projwfc_out_name'] = 'bn.projwfc.out'
names['pdos_prefix'] = 'bn'
file.Set_files_atributes(names)
file.Load()


file.Bands_pdos()
file.orbitals_numbers()
plot.bands_project(file.bands_proj,states = (9,13),E_min= -20, E_max= 20, dE = 5, cmap = 'bone_r')
plot.bands_project(file.bands_proj, corretalion = True, states2 = range(9,17), states1 = range(1,9),E_min= -20, E_max= 20, dE = 5, cbar_label=['B','N'],cbar_tick= [-1,1])



file.Pdos_kresolved()


plot.kresolved(file.pdos_kresolved_plot,states = (10,11,12,14,15,16),E_min= -18, E_max= 14, dE = 2,cmap = 'magma')

file.k_axis

file.orbitals_numbers()


names = {}
names['files_folder'] = 'Examples/graphene/pp'
names['ppin_name'] = 'pp2D.in'
names['charge2D_dat_name'] = 'charge-density2D.out'
file.Set_files_atributes(names)
file.Load()
file.Charge2D()

plot.charge_density2D(file.charge2D_data, xlim = (0,15,5),ylim = (0,15,5))



names = {}
names['files_folder'] = 'Examples/GeS/bulk/bands/'
file.Set_files_atributes(names)
file.Load()
file.Bands_files()
plot.bands(file.bands, file.k_points_letter, file.k_path, E_min = -3, E_max = 3, dpi = 300, title = ' GeS Band Structure ')


names = {}
names['files_folder'] = 'Examples/GeS/bulk/bands_k/'
file.Set_files_atributes(names)
file.Load()
file.Bands_pdos()
file.orbitals_numbers()
state = [2,3,4,12,12,13,15,16,17,19,20,21,28,29,30,37,38,39,41,42,43,45,46,47]

plot.bands_project(file.bands_proj,states = state,E_min= -15, E_max= 15, dE = 5, cmap = 'bone_r')
Ge = [1,2,3,4,5,6,7,8,9,18,19,20,21,22,23,24,25,44,45,46,47,48,49,50,51,52 ]
S = [10,11,12,13,14,15,16,17,36,37,38,39,40,41,42,43]

plot.bands_project(file.bands_proj, corretalion = True, states2 = S, states1 = Ge,E_min= -5, E_max= 5, dE = 1, 
                   cbar_label=['S','Ge'],cbar_tick= [-1,1], dpi = 300, font_axis= 40, font_label = 40, linewidth = 12, subplotsize = (20,10))



names = {}
names['files_folder'] = 'Examples/GeS/distorcido'
names['scf_in_name'] = 'scf.in'
names['scf_out_name'] = 'scf.out'
names['nscf_out_name'] = 'nscf.out'
names['bands_in_name'] = 'bands.in'
names['bands_out_name'] = 'bands.out'
names['projwfc_out_name'] = 'projwfc.out'
file.Set_files_atributes(names)
file.Load()

file.Pdos_files()
plot.pdos_atoms(file.pdos_per_atoms,subplotsize=(8,6),dos_max= 100, E_max = 8)
