#! /usr/bin/env python

# ===========================================================================
# File: calc_density.py
# Author: Kelly Tran
# Purpose: This python script reads in the volume of the simulation box
#     computed at each timestep from molecular dynamics simulations
#     and calculates the density of water.
# ===========================================================================

from contextlib import contextmanager
from glob import iglob
from math import sqrt
from sys import exit
import numpy as np

@contextmanager
def multi_file_manager(files, mode='rt'):
    files = [open(file, mode) for file in files]
    yield files
    for file in files:
        file.close()

# Specify format style and units
fmt_print = '#{:>11} {:>12} {:>12} {:>12}'
fmt_header = '\n#{:>11}' + ' {:>12}'*4
fmt_data = '\n' + '{:>12.1f} '*3 + '{:>12.5f} {:>12.5f}'
units = '(g/cm^3)'

# Print density to screen
print()
print(fmt_print.format('Temp.' , 'Pressure' , '<Volume>' , '<Density>'))
print(fmt_print .format('(K)' , '(atm)' , '(A^3)' , units ))

# Write density to data file for plotting
with open('dens_ssmp_m4_1atm.dat', 'wt') as header:
    header.write('# ssmp_m4\n')
    header.write(fmt_header \
    .format('Temp.' , 'Pressure' , '<Volume>' , '<Density>' , 'SD'))
    header.write(fmt_header \
        .format('(K)' , '(atm)' , '(A^3)' , units , ''))

# Primary loop: loop over all temperatures
# ===========================================================================
for temp in [238, 258, 268, 278, 298, 318, 338]:
    
    DATA_FILE_PATTERN = 'inp/volume/m4-t' + str(temp) +'-p00001*'
    MIN_DATA_FILES = 2
    
    TEMPERATURE = temp                      # (K)
    PRESSURE = 1                            # (atm)
    CONV_FACT = 1.66                        # (amu/A^3) to (g/cm^3)
    MASS = (512 * (15.9994 + 1.008 +1.008)) # (512 molecules * (MW of water))

    with multi_file_manager(iglob(DATA_FILE_PATTERN)) as datfiles:
        num_files = len(datfiles)
        if num_files < MIN_DATA_FILES:
            print('Less than {} .dat files were found to process, '
                  'terminating.'.format(MIN_DATA_FILES))
            exit(1)
    
        vol_avgs = []
        densities = []
        divisor = float(num_files-1)

        # Secondary loop: loop over all data files
        # ===================================================================
        for file in datfiles: 

            # Calculate average volume from 4ns trajectory
            volumes = [float(n) for n in file.read().split()]
            vol_avg = float(sum(volumes)) / len(volumes)
            vol_avgs.append(vol_avg)

            # Calculate density from 4ns trajectory
            density = CONV_FACT * MASS / vol_avg
            densities.append(density)
            print('{:>12.1f} {:>12.1f} {:>12.1f} {:>12.5f}' \
                .format(TEMPERATURE , PRESSURE , vol_avg , density))
            
            # Calculate mean volume from all runs
            vol_mean= float(sum(vol_avgs)) / len(vol_avgs)
            
            # Calculate average density and standard dev from all runs
            mean= float(sum(densities)) / len(densities)
            means_diff_sq = ((density-mean)**2 for density in densities)
            std_dev = sqrt(sum(means_diff_sq) / divisor)
        # ===================================================================

    print('')
    with open('dens_ssmp_m4_1atm.dat', 'at') as dens_avg:
        dens_avg.write(fmt_data \
            .format(TEMPERATURE , PRESSURE , vol_mean, mean , std_dev))
  
# ===========================================================================