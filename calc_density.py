#! /usr/bin/env python

# ========================================================================================
# File: calc_density.py
# Author: Kelly Tran
# Purpose:
# This python script reads in the volume of the simulation box
# computed at each timestep from molecular dynamics simulations
# and calculates the density of water.
# ========================================================================================

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

# Write density to data file for plotting
with open('dens_ssmp_m4_1atm.dat', 'wt') as header:
    header.write('# ssmp-m4\n')
    header.write('#     Temp       Pres        Volu_avg         Box Len        Dens_avg\n')
    header.write('#      (K)      (atm)          (A**3)             (A)       (g/cm**3)\n')

# Print density to screen
print('#     Temp       Pres        Volu_avg        Dens_avg')
print('#      (K)      (atm)          (A**3)       (g/cm**3)')

for temp in [238, 258, 268, 278, 298, 318, 338]:
    
    DATA_FILE_PATTERN = 'inp_volu/m4-t' + str(temp) +'-p00001*'
    MIN_DATA_FILES = 2
    
    TEMPERATURE = temp                       # (K)
    PRESSURE = 1                             # (atm)
    CONV_FACT = 1.66                         # (amu/A**3) to (g/cm**3)
    MASS = (512 * (15.9994 + 1.008 +1.008))  # (512 molecules * (15.9994 + 1.008 +1.008))
    
    with multi_file_manager(iglob(DATA_FILE_PATTERN)) as datfiles:
        num_files = len(datfiles)
        if num_files < MIN_DATA_FILES:
            print('Less than {} .dat files were found to process, '
                  'terminating.'.format(MIN_DATA_FILES))
            exit(1)
    
        vol_avgs = []
        densities = []
        divisor = float(num_files-1)  # Bessel's correction for sample standard dev
        for file in datfiles:  # main processing loop
            # Calculate average volume from 4ns trajectory
            volumes = [float(n) for n in file.read().split()]
            vol_avg = float(sum(volumes)) / len(volumes)
            vol_avgs.append(vol_avg)

            # Calculate density from 4ns trajectory
            density = CONV_FACT * MASS / vol_avg
            densities.append(density)
            print('{:>10.1f} {:>10.1f} {:>15.1f} {:>15.5f}' \
                .format(TEMPERATURE , PRESSURE , vol_avg , density))
            
            # Calculate mean volume from all runs
            vol_mean= float(sum(vol_avgs)) / len(vol_avgs)
            box_len = vol_mean**(1./3.)
            
            # Calculate average density and standard dev from all runs
            mean= float(sum(densities)) / len(densities)
            means_diff_sq = ((density-mean)**2 for density in densities)
            std_dev = sqrt(sum(means_diff_sq) / divisor)
    
    print('')
    with open('dens_ssmp_m4_1atm.dat', 'at') as dens_avg:
        dens_avg.write('{:>10.1f} {:>10.1f} {:>15.1f} {:>15.5f} {:>15.5f} {:>15.5f}\n' \
            .format(TEMPERATURE , PRESSURE , vol_mean, box_len , mean , std_dev))
  
