#! /usr/bin/env python

# ===========================================================================
# File: calc_diffusion.py
# Author: Kelly Tran
# Purpose: This python script reads in the mean-square-displacement (MSD)
#     of oxygen atoms computed at each timestep from molecular dynamics 
#     simulations and calculates the diffusion coefficient of water 
#     using the Einstein relation.
# ===========================================================================

from contextlib import contextmanager
from glob import iglob
from sys import exit
from math import sqrt
import numpy as np
import math

@contextmanager
def multi_file_manager(files, mode='rt'):
    files = [open(file, mode) for file in files]
    yield files
    for file in files:
        file.close()

# Specify format style and units
fmt_print = '#{:>11} {:>12} {:>12}'
fmt_header = '\n#{:>11}' + ' {:>12}'*5
fmt_data = '\n{:>12.1f} {:>12.1f}' + ' {:>12.5f}'*4
units = '(1E-09 m^2/s)'

# Read temperature dependent viscosity and write to dictionary
# ===========================================================================
viscs = {}
CONV_FACT = 1E-3    # (cP to kg/(m*s))
with open('inp/corr_diff/visc.dat') as file:
    lines = file.readlines()[2:]
    for line in lines:
        (T, visc) = line.split()
        viscs[int(T)] = float(visc)*CONV_FACT
# ===========================================================================

# Print diffusion coefficient to screen
print()
print(fmt_print.format('Temp.' , 'Pressure' , '<D>'))
print(fmt_print.format('(K)' , '(atm)' , units ))

# Write diffusion coefficient to data file for plotting
with open('diff_ssmp_m4_1atm.dat', 'wt') as header:
    header.write('# ssmp_m4\n')
    header.write(fmt_header \
    .format('Temp.' , 'Pressure' , '<D>' , 'C' , '<D>+C' , 'SD'))
    header.write(fmt_header \
        .format('(K)' , '(atm)' , units , units , units , ''))

# Primary loop: loop over all temperatures
# ===========================================================================
for temp in [238, 258, 268, 278, 298, 318, 338]:
    
    DATA_FILE_PATTERN = 'inp/msd/m4-t' + str(temp) +'-p00001*'
    MIN_DATA_FILES = 2
    
    TEMPERATURE = temp    # (K)
    PRESSURE = 1          # (atm)
    k_B = 1.380648E-23    # Boltzmann constant (J/K)
    L = 24.83689E-10      # box length (m)
    CONV_FACT3 = 100      # (1E-11 m^2/s to 1E-9 m^2/s)
                                                  
    # Correction term for the simulation size dependence (10^-9 m^2/s)
    constant = (2.837297*k_B)/(6*math.pi)
    C = constant * (TEMPERATURE/(viscs[TEMPERATURE]*L)) / 1E-9                          
     
    with multi_file_manager(iglob(DATA_FILE_PATTERN)) as datfiles:
        num_files = len(datfiles)
        if num_files < MIN_DATA_FILES:
            print('Less than {} .dat files were found to process, '
                  'terminating.'.format(MIN_DATA_FILES))
            exit(1)
    
        diffs = []
        divisor = float(num_files-1)

        # Secondary loop: loop over all data files
        # ===================================================================
        for file in datfiles:
            
            # Compute the slope of MSD vs. time
            data = [line.split() for line in file]
            xd = [float(value[0])*10 for value in data]    # Convert to ns
            yd = [float(value[1]) for value in data]

            # Determine the best fit line
            par = np.polyfit(xd, yd, 1, full=True)
            slope=par[0][0]
            intercept=par[0][1]
            
            # Calculate diffusion coefficient
            diff = slope*CONV_FACT3/6
            print('{:>12.1f} {:>12.1f} {:>12.5f}' \
                .format(TEMPERATURE , PRESSURE , diff))
            diffs.append(diff)
        # ===================================================================

        # Calculate average diffusion coefficient with correction term and SD
        mean= float(sum(diffs)) / len(diffs)
        means_diff_sq = ((diff-mean)**2 for diff in diffs)
        std_dev = sqrt(sum(means_diff_sq) / divisor)
        mean_corr = mean + C
    
    print('')
    with open('diff_ssmp_m4_1atm.dat', 'at') as diff_avg:
        diff_avg.write(fmt_data \
            .format(TEMPERATURE , PRESSURE , mean , C , mean_corr , std_dev))

# END
