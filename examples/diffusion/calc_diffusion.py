#! /usr/bin/env python

# ========================================================================================
# File: calc_diffusion.py
# Author: Kelly Tran
# Purpose:
# This python script reads in the mean-square-displacement (MSD) of oxygen atoms
# computed at each timestep from molecular dynamics simulations
# and calculates the diffusion coefficient of water using the Einstein relation.
# ========================================================================================

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

# Read temperature dependent viscosity and write to dictionary
viscs = {}
CONV_FACT = 1E-3                                # (cP to kg m**-1 s**-1)
with open('inp_corr_diff/visc.dat') as file:
    lines = file.readlines()[2:]
    for line in lines:
        (T, visc) = line.split()
        viscs[int(T)] = float(visc)*CONV_FACT

# Write diffusion coefficient to data file for plotting
with open('diff_ssmp_m4_1atm.dat', 'wt') as header:
    header.write('# ssmp-m4\n')
    header.write('#     Temp       Pres          <diff>               C        <diff>+C\n')
    header.write('#      (K)      (atm)  (1E-09 m**2/s)  (1E-09 m**2/s)  (1E-09 m**2/s)\n')

# Print diffusion coefficient to screen
print('#     Temp       Pres            diff')
print('#      (K)      (atm)  (1E-09 m**2/s)')

for temp in [238, 258, 268, 278, 298, 318, 338]:
    
    DATA_FILE_PATTERN = 'inp_msd/m4-t' + str(temp) +'-p00001*'
    MIN_DATA_FILES = 2
    
    TEMPERATURE = temp                             # (K)
    PRESSURE = 1                                   # (atm)
    k_B = 1.380648E-23                             # Boltzmann constant (J/K)
    L = 24.83689E-10                               # box length (m)
    CONV_FACT3 = 100                               # (1E-11 m**2/s to 1E-9 m**2/s)
                                                  
    # Correction term for the simulation size dependence
    constant = (2.837297*k_B)/(6*math.pi)
    C = constant * (TEMPERATURE/(viscs[TEMPERATURE]*L)) / 1E-9  # correction term (10**-9 m**2/s)                        
     
    with multi_file_manager(iglob(DATA_FILE_PATTERN)) as datfiles:
        num_files = len(datfiles)
        if num_files < MIN_DATA_FILES:
            print('Less than {} .dat files were found to process, '
                  'terminating.'.format(MIN_DATA_FILES))
            exit(1)
    
        diffs = []
        divisor = float(num_files-1)  # Bessel's correction for sample standard dev
        for file in datfiles:  # main processing loop
            
            # Compute the slope of MSD vs. time
            data = [line.split() for line in file]
            xd = [float(value[0])*10 for value in data]            # (to ns)
            yd = [float(value[1]) for value in data]
            # Determine the best fit line
            par = np.polyfit(xd, yd, 1, full=True)
            slope=par[0][0]
            intercept=par[0][1]
            
            # Calculate diffusion coefficient
            diff = slope*CONV_FACT3/6
            print('{:>10.1f} {:>10.1f} {:>15.5f}' \
                .format(TEMPERATURE , PRESSURE , diff))
            diffs.append(diff)
            
        # Calculate average diffusion coefficient with correction term and standard dev
        mean= float(sum(diffs)) / len(diffs)
        means_diff_sq = ((diff-mean)**2 for diff in diffs)
        std_dev = sqrt(sum(means_diff_sq) / divisor)
        mean_corr = mean + C
    
    print('')
    with open('diff_ssmp_m4_1atm.dat', 'at') as diff_avg:
        diff_avg.write('{:>10.1f} {:>10.1f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f}\n' \
            .format(TEMPERATURE , PRESSURE , mean , C , mean_corr , std_dev))
  
