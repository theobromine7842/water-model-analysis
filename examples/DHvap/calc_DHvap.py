#! /usr/bin/env python

# ========================================================================================
# File: calc_DHvap.py
# Author: Kelly Tran
# Purpose:
# This python script reads in the potential energy (U) of the system
# computed at each timestep from molecular dynamics simulations
# and calculates the heat of vaporization of water.
# ========================================================================================

from contextlib import contextmanager
from glob import iglob
from math import sqrt
from sys import exit

@contextmanager
def multi_file_manager(files, mode='rt'):
    files = [open(file, mode) for file in files]
    yield files
    for file in files:
        file.close()

# Read in temperature dependent corrections and write to dictionary
# Correction term for vibrational effects
C_vibs = {}
with open('inp_corr_DHvap/C_vib.dat') as file:
    lines = file.readlines()[2:]
    for line in lines:
        (T, C_vib) = line.split()
        C_vibs[int(T)] = float(C_vib)

# Correction term for nonideal gas effects
C_nis = {}
with open('inp_corr_DHvap/C_ni.dat') as file:
    lines = file.readlines()[2:]
    for line in lines:
        (T, C_ni) = line.split()
        C_nis[int(T)] = float(C_ni)

# Write heat of vaporization to data file for plotting
with open('DHvap_ssmp_m4_1atm.dat', 'wt') as header:
    header.write('# ssmp-m4\n')
    header.write('#     Temp       Pres        <DH_vap>               C      <DH_vap+C>\n')
    header.write('#      (K)      (atm)      (kcal/mol)      (kcal/mol)      (kcal/mol)\n')

# Print heat of vaporization to screen
print('#     Temp       Pres             <U>          DH_vap')
print('#      (K)      (atm)      (kcal/mol)      (kcal/mol)')

for temp in [238, 258, 268, 278, 298, 318, 338]:
    
    DATA_FILE_PATTERN = 'inp_U/m4-t' + str(temp) +'-p00001*'
    MIN_DATA_FILES = 2
    
    TEMPERATURE = temp                                            # (K)
    PRESSURE = 1                                                  # (atm)
    N = 512                                                       # number of molecules
    R = 0.0019872036                                              # gas constant (kcal K**-1 mol**-1)
                                              
    # Calculate DHvap correction term for rigid, non-polarizable water models                             
    mu_gas = 1.854989                                             # experimental dipole gas (D)
    mu_liquid = 2.28                                              # SSMP-M# dipole (D)
    alpha_gas = 1.6633E-40                                        # mean polarizability of water molecule in the gas phase (F m**2)
    CONV_FACT = (1/(3.33564E-30)**2) * (4184/6.022E+23)           # (C m)**2 to (D**2) AND (J) to (kcal/mol)    
    C_pol = -0.5*((mu_gas - mu_liquid))**2/(alpha_gas*CONV_FACT)  # (kcal/mol)
    C_vib = C_vibs[TEMPERATURE]                                   # (kcal/mol)
    C_ni = C_nis[TEMPERATURE]                                     # (kcal/mol) @T=238K
    C_x = 0                                                       # small and therefore neglected
    C = C_pol + C_vib + C_ni + C_x                                # correction term
     
    with multi_file_manager(iglob(DATA_FILE_PATTERN)) as datfiles:
        num_files = len(datfiles)
        if num_files < MIN_DATA_FILES:
            print('Less than {} .dat files were found to process, '
                  'terminating.'.format(MIN_DATA_FILES))
            exit(1)
    
        DH_vaps = []
        divisor = float(num_files-1)  # Bessel's correction for sample standard dev
        
        for file in datfiles:  # main processing loop
            # Calculate average U
            Us = [float(n) for n in file.read().split()]
            U_avg = float(sum(Us)) / len(Us)
            
            # Calculate DHvap
            DH_vap = - U_avg / N + R * TEMPERATURE
            print('{:>10.1f} {:>10.1f} {:>15.1f} {:>15.5f}' \
                .format(TEMPERATURE , PRESSURE , U_avg , DH_vap))
            DH_vaps.append(DH_vap)
            
        # Calculate average DH_vap with correction term and standard dev
        mean= float(sum(DH_vaps)) / len(DH_vaps)
        means_diff_sq = ((DH_vap-mean)**2 for DH_vap in DH_vaps)
        std_dev = sqrt(sum(means_diff_sq) / divisor)
        mean_corr = mean + C

    print('')
    with open('DHvap_ssmp_m4_1atm.dat', 'at') as DH_vap_avg:
        DH_vap_avg.write('{:>10.1f} {:>10.1f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f}\n' \
            .format(TEMPERATURE , PRESSURE , mean , C , mean_corr , std_dev))
        
