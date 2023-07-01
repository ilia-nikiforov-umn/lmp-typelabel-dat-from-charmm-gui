"""
Created on Sun Aug 7 2021

Reposition MoS2 for clearance

Input:
    Name of LAMMPS data file
    
Output: 
    Structure data
    
@author: Moon-ki Choi, Ilia Nikiforov
"""
import numpy as np

##################################################################
# READ LAMMPS BONDED DATA FILE
##################################################################
def read_LAMMPS_bonded(filename):
    """ Read LAMMPS topology file """
    fopen = open(filename, "r")
    
    # End of box data flag 
    eof_box = 0
    while (eof_box == 0):
        curr_line = fopen.readline()
        curr_line_split = curr_line.split()
        if (len(curr_line_split) > 1):
            # Number of atom 
            if (curr_line_split[1] == 'atoms'): 
                num_atom = int(curr_line_split[0])
            # Number of bonds
            if (curr_line_split[1] == 'bonds'): 
                num_bond = int(curr_line_split[0])
            # Number of angles
            if (curr_line_split[1] == 'angles'): 
                num_angle = int(curr_line_split[0])
            # Number of dihedral
            if (curr_line_split[1] == 'dihedrals'): 
                num_dihedral = int(curr_line_split[0])
            # Number of atom types 
            if ((curr_line_split[1] == 'atom') and (curr_line_split[2] == 'types')): 
                num_atom_type = int(curr_line_split[0])
            # Number of bond types 
            if ((curr_line_split[1] == 'bond') and (curr_line_split[2] == 'types')): 
                num_bond_type = int(curr_line_split[0])
            # Number of bond types 
            if ((curr_line_split[1] == 'angle') and (curr_line_split[2] == 'types')): 
                num_angle_type = int(curr_line_split[0])
            # Number of bond types 
            if ((curr_line_split[1] == 'dihedral') and (curr_line_split[2] == 'types')): 
                num_dihedral_type = int(curr_line_split[0]) 
        if (len(curr_line_split) > 2):
            # Read box size  
            if ((curr_line_split[2] == 'xlo') and (curr_line_split[3] == 'xhi')): 
                xlo = np.double(curr_line_split[0]); xhi = np.double(curr_line_split[1])
            if ((curr_line_split[2] == 'ylo') and (curr_line_split[3] == 'yhi')): 
                ylo = np.double(curr_line_split[0]); yhi = np.double(curr_line_split[1])
            if ((curr_line_split[2] == 'zlo') and (curr_line_split[3] == 'zhi')): 
                zlo = np.double(curr_line_split[0]); zhi = np.double(curr_line_split[1])
                # End of reading box size 
                eof_box = 1
        
    # Read mass data
    eof_mass = 0 # End of mass data flag 
    while (eof_mass == 0):
        curr_line = fopen.readline()
        curr_line_split = curr_line.split()
        # Masses d
        if (len(curr_line_split) > 0): 
            if (curr_line_split[0] == 'Masses'): 
                fopen.readline()
                mass = np.zeros(num_atom_type)
                labels = []
                for i in range(num_atom_type):
                    curr_line = fopen.readline()
                    curr_line_split = curr_line.split()
                    mass[i] = np.double(curr_line_split[1])
                    labels.append(curr_line_split[3])
                # End of mass data
                eof_mass = 1
    
    # Read pair coefficients
    eof_pair_coeffs = 0 # End of pair coefficients 
    while (eof_pair_coeffs == 0): 
        curr_line = fopen.readline()
        curr_line_split = curr_line.split()
        # Masses d
        if (len(curr_line_split) > 0): 
            if (curr_line_split[0] == 'Pair'): 
                fopen.readline()
                pair_coeff = np.zeros((num_atom_type,2))
                for i in range(num_atom_type):
                    curr_line = fopen.readline()
                    curr_line_split = curr_line.split()
                    pair_coeff[i,0] = np.double(curr_line_split[1]) # Epsilon 
                    pair_coeff[i,1] = np.double(curr_line_split[2]) # Sigma 
                # End of pair coefficients 
                eof_pair_coeffs = 1    
                
    # Read atom data 
    eof_atom = 0 # End of atom 
    while (eof_atom == 0): 
        curr_line = fopen.readline()
        curr_line_split = curr_line.split()
        # Masses d
        if (len(curr_line_split) > 0): 
            if (curr_line_split[0] == 'Atoms'): 
                fopen.readline()
                atom = np.zeros((num_atom,6))
                for i in range(num_atom):
                    curr_line = fopen.readline()
                    curr_line_split = curr_line.split()
                    atom[i,0] = np.double(curr_line_split[1]) # molecule-tag 
                    atom[i,1] = np.double(curr_line_split[2]) # type 
                    atom[i,2] = np.double(curr_line_split[3]) # charge
                    atom[i,3] = np.double(curr_line_split[4]) # x_position
                    atom[i,4] = np.double(curr_line_split[5]) # y_position
                    atom[i,5] = np.double(curr_line_split[6]) # z_position
                # End of atoms
                eof_atom = 1  

    if num_bond_type == 0:                   
        bond_coeff = []
        bond_label = []
        bond = []
    else:
        # Read bond coefficients 
        eof_bond_coeff = 0 # End of bond_coeffs 
        while (eof_bond_coeff == 0): 
            curr_line = fopen.readline()
            curr_line_split = curr_line.split()
            # Masses d
            if (len(curr_line_split) > 0): 
                if (curr_line_split[0] == 'Bond'): 
                    fopen.readline()
                    bond_coeff = np.zeros((num_bond_type,2))
                    bond_label = []
                    for i in range(num_bond_type):
                        curr_line = fopen.readline()
                        curr_line_split = curr_line.split()
                        bond_coeff[i,0] = np.double(curr_line_split[1]) # Bond coefficient 
                        bond_coeff[i,1] = np.double(curr_line_split[2]) # Eq_distance 
                        bond_label.append([curr_line_split[4],curr_line_split[5]])  # Label connection for bond type
                    # End of bond_coeffs 
                    eof_bond_coeff = 1  
        
        
        # Read bonds 
        eof_bond = 0 # End of bond_coeffs 
        while (eof_bond == 0): 
            curr_line = fopen.readline()
            curr_line_split = curr_line.split()
            if (len(curr_line_split) > 0): 
                if (curr_line_split[0] == 'Bonds'): 
                    fopen.readline()
                    bond = np.zeros((num_bond,3))
                    for i in range(num_bond):
                        curr_line = fopen.readline()
                        curr_line_split = curr_line.split()
                        bond[i,0] = np.double(curr_line_split[1]) # Bond type
                        bond[i,1] = np.double(curr_line_split[2]) # ID_1
                        bond[i,2] = np.double(curr_line_split[3]) # ID_2
                    # End of bond_coeffs 
                    eof_bond = 1  
                
    if num_angle_type == 0:
        angle_coeff=[]
        angle_label=[]
        angle=[]
    else:
        # Read angle coefficients 
        eof_angle_coeff = 0 # End of angle coeffs 
        while (eof_angle_coeff == 0): 
            curr_line = fopen.readline()
            curr_line_split = curr_line.split()
            if (len(curr_line_split) > 0): 
                if (curr_line_split[0] == 'Angle'): 
                    fopen.readline()
                    angle_coeff = np.zeros((num_angle_type,4))
                    angle_label = []
                    for i in range(num_angle_type):
                        curr_line = fopen.readline()
                        curr_line_split = curr_line.split()
                        angle_coeff[i,0] = np.double(curr_line_split[1]) # Angle coefficient
                        angle_coeff[i,1] = np.double(curr_line_split[2]) # Eq_angle 
                        angle_coeff[i,2] = np.double(curr_line_split[3]) # CHARMM_coeff1
                        angle_coeff[i,3] = np.double(curr_line_split[4]) # CHARMM_coeff2
                        angle_label.append([curr_line_split[6],curr_line_split[7],curr_line_split[8]]) # Label connection for angle type
                    # End of angle coeffs 
                    eof_angle_coeff = 1  
        
        # Read angle 
        eof_angle = 0 # End of angle
        while (eof_angle == 0): 
            curr_line = fopen.readline()
            curr_line_split = curr_line.split()
            if (len(curr_line_split) > 0): 
                if (curr_line_split[0] == 'Angles'): 
                    fopen.readline()
                    angle = np.zeros((num_angle,4))
                    for i in range(num_angle):
                        curr_line = fopen.readline()
                        curr_line_split = curr_line.split()
                        angle[i,0] = np.double(curr_line_split[1]) # Angle type
                        angle[i,1] = np.double(curr_line_split[2]) # ID_1 
                        angle[i,2] = np.double(curr_line_split[3]) # ID_2
                        angle[i,3] = np.double(curr_line_split[4]) # ID_3
                    # End of angle
                    eof_angle = 1 
                
    if num_dihedral_type == 0:
        dihedral_coeff=[]
        dihedral_label=[]
        dihedral=[]
    else:
        # Read dihedral coefficients  
        eof_dihedral_coeff = 0 # End of dihedral coefficient flag
        while (eof_dihedral_coeff == 0): 
            curr_line = fopen.readline()
            curr_line_split = curr_line.split()
            if (len(curr_line_split) > 0): 
                if (curr_line_split[0] == 'Dihedral'): 
                    fopen.readline()
                    dihedral_coeff = np.zeros((num_dihedral_type,4))
                    dihedral_label = []
                    for i in range(num_dihedral_type):
                        curr_line = fopen.readline()
                        curr_line_split = curr_line.split()
                        dihedral_coeff[i,0] = np.double(curr_line_split[1]) # Coeff1
                        dihedral_coeff[i,1] = np.double(curr_line_split[2]) # Coeff2 
                        dihedral_coeff[i,2] = np.double(curr_line_split[3]) # Coeff3
                        dihedral_coeff[i,3] = np.double(curr_line_split[4]) # Coeff4
                        dihedral_label.append([curr_line_split[6],curr_line_split[7],curr_line_split[8],curr_line_split[9]]) # Label connection for angle type
                    # End of dihedral coefficients flag 
                    eof_dihedral_coeff = 1  
        
        # Read dihedral
        eof_dihedral = 0 # End of dihedral flag
        while (eof_dihedral == 0): 
            curr_line = fopen.readline()
            curr_line_split = curr_line.split()
            if (len(curr_line_split) > 0): 
                if (curr_line_split[0] == 'Dihedrals'): 
                    fopen.readline()
                    dihedral = np.zeros((num_dihedral,5))
                    for i in range(num_dihedral):
                        curr_line = fopen.readline()
                        curr_line_split = curr_line.split()
                        dihedral[i,0] = np.double(curr_line_split[1]) # Dihedral coefficient
                        dihedral[i,1] = np.double(curr_line_split[2]) # ID_1
                        dihedral[i,2] = np.double(curr_line_split[3]) # ID_2
                        dihedral[i,3] = np.double(curr_line_split[4]) # ID_3
                        dihedral[i,4] = np.double(curr_line_split[5]) # ID_3
                    # End of dihedral flag 
                    eof_dihedral = 1              
        
    fopen.close()

    return num_atom,num_bond,num_angle,num_dihedral,num_atom_type,num_bond_type,num_angle_type,num_dihedral_type, \
           xlo,xhi,ylo,yhi,zlo,zhi,mass,labels,pair_coeff,atom,bond_coeff,bond_label,bond,angle_coeff,angle_label,angle,dihedral_coeff,dihedral_label,dihedral