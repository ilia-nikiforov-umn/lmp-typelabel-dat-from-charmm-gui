# -*- coding: utf-8 -*-

"""
Created on Sun Aug 7 2021

Reconstruction of data file with labels 
    
@author: Moon-ki Choi
"""

from .read_LAMMPS_data import read_LAMMPS_bonded
from .write_LAMMPS_data import write_LAMMPS_bonded_label_v2

def dump_dat(dat_in,dat_out):
    """ Read LAMMPS bonded file """ 
    print_flag = 0
    num_atom,num_bond,num_angle,num_dihedral,num_atom_type,num_bond_type,num_angle_type,num_dihedral_type, \
            xlo,xhi,ylo,yhi,zlo,zhi,mass,labels,pair_coeff,atom,bond_coeff,bond_label,bond,angle_coeff,angle_label,angle,dihedral_coeff,dihedral_label,dihedral,element_iff_dict = read_LAMMPS_bonded(dat_in,print_flag)

    """ Create label list for atoms """
    label_list = []
    for i in range(num_atom):
        # Find label for each atom based on atom type 
        label = labels[int(atom[i,1]-1)] 
        label_list.append(label)

    """ Overlapping label checking process """
    """ This section will find overlapping label type and remove it """
    check_flag = 1
    if (check_flag == 1): 
        # Bond check 
        bond_label_tmp = []
        bond_coeff_tmp = []
        num_bond_type_tmp = num_bond_type
        for i in range(num_bond_type_tmp):
            ovr_flag = 0 
            curr1 = bond_label[i][0]
            curr2 = bond_label[i][1]
            coeff1 = bond_coeff[i][0]
            coeff2 = bond_coeff[i][1]
            if (i > 0):
                for j in range(i):
                    tmp1 = bond_label[j][0]
                    tmp2 = bond_label[j][1]
                    # Overlapped case 
                    if curr1 == tmp1 and curr2 == tmp2:
                        # If label is already replaced, continue
                        if ovr_flag == 1:
                            continue
                        ovr_flag = 1
                        # Reduce num_bond
                        num_bond_type = num_bond_type - 1
            if (ovr_flag==0):
                bond_label_tmp.append([curr1,curr2])
                bond_coeff_tmp.append([coeff1,coeff2])
        # Replace type_id in bond with new one
        for i in range(num_bond):
            tmd_id = int(bond[i][0])
            curr1 = bond_label[tmd_id-1][0]
            curr2 = bond_label[tmd_id-1][1]
            for j in range(len(bond_label_tmp)):
                tmp1 = bond_label_tmp[j][0]
                tmp2 = bond_label_tmp[j][1]
                if curr1 == tmp1 and curr2 == tmp2:
                    bond[i][0] = j+1
        bond_label = bond_label_tmp
        bond_coeff = bond_coeff_tmp
                        
        # Angle check 
        angle_label_tmp = []
        angle_coeff_tmp = []
        num_angle_type_tmp = num_angle_type
        for i in range(num_angle_type_tmp):
            ovr_flag = 0 
            curr1 = angle_label[i][0]
            curr2 = angle_label[i][1]
            curr3 = angle_label[i][2]
            coeff1 = angle_coeff[i][0]
            coeff2 = angle_coeff[i][1]
            coeff3 = angle_coeff[i][2]
            coeff4 = angle_coeff[i][3]
            if (i > 0):
                for j in range(i):
                    tmp1 = angle_label[j][0]
                    tmp2 = angle_label[j][1]
                    tmp3 = angle_label[j][2]
                    # Overlapped case 
                    if curr1 == tmp1 and curr2 == tmp2 and curr3 == tmp3:
                        # If label is already replaced, continue
                        if ovr_flag == 1:
                            continue
                        ovr_flag = 1
                        # Reduce num_angle
                        num_angle_type = num_angle_type - 1
            if (ovr_flag==0):
                angle_label_tmp.append([curr1,curr2,curr3])
                angle_coeff_tmp.append([coeff1,coeff2,coeff3,coeff4])
        # Replace type_id in angle with new one
        for i in range(num_angle):
            tmd_id = int(angle[i][0])
            curr1 = angle_label[tmd_id-1][0]
            curr2 = angle_label[tmd_id-1][1]
            curr3 = angle_label[tmd_id-1][2]
            for j in range(len(angle_label_tmp)):
                tmp1 = angle_label_tmp[j][0]
                tmp2 = angle_label_tmp[j][1]
                tmp3 = angle_label_tmp[j][2]
                if curr1 == tmp1 and curr2 == tmp2 and curr3 == tmp3:
                    angle[i][0] = j+1
        angle_label = angle_label_tmp
        angle_coeff = angle_coeff_tmp
        
        # Dihedral check 
        dihedral_label_tmp = []
        dihedral_coeff_tmp = []
        num_dihedral_type_tmp = num_dihedral_type
        for i in range(num_dihedral_type_tmp):
            ovr_flag = 0 
            curr1 = dihedral_label[i][0]
            curr2 = dihedral_label[i][1]
            curr3 = dihedral_label[i][2]
            curr4 = dihedral_label[i][3]
            coeff1 = dihedral_coeff[i][0]
            coeff2 = dihedral_coeff[i][1]
            coeff3 = dihedral_coeff[i][2]
            coeff4 = dihedral_coeff[i][3]
            if (i > 0):
                for j in range(i):
                    tmp1 = dihedral_label[j][0]
                    tmp2 = dihedral_label[j][1]
                    tmp3 = dihedral_label[j][2]
                    tmp4 = dihedral_label[j][3]
                    # Overlapped case 
                    if curr1 == tmp1 and curr2 == tmp2 and curr3 == tmp3 and curr4 == tmp4:
                        # If label is already replaced, continue
                        if ovr_flag == 1:
                            continue
                        ovr_flag = 1
                        # Reduce num_dihedral
                        num_dihedral_type = num_dihedral_type - 1
            if (ovr_flag==0):
                dihedral_label_tmp.append([curr1,curr2,curr3,curr4])
                dihedral_coeff_tmp.append([coeff1,coeff2,coeff3,coeff4])
        # Replace type_id in dihedral with new one
        for i in range(num_dihedral):
            tmd_id = int(dihedral[i][0])
            curr1 = dihedral_label[tmd_id-1][0]
            curr2 = dihedral_label[tmd_id-1][1]
            curr3 = dihedral_label[tmd_id-1][2]
            curr4 = dihedral_label[tmd_id-1][3] 
            for j in range(len(dihedral_label_tmp)):
                tmp1 = dihedral_label_tmp[j][0]
                tmp2 = dihedral_label_tmp[j][1]
                tmp3 = dihedral_label_tmp[j][2]
                tmp4 = dihedral_label_tmp[j][3]
                if curr1 == tmp1 and curr2 == tmp2 and curr3 == tmp3 and curr4 == tmp4:
                    #print("found",j+1)
                    dihedral[i][0] = j+1
        dihedral_label = dihedral_label_tmp
        dihedral_coeff = dihedral_coeff_tmp
        # Improper check 

    """ Write LAMMPS data file with labels """     
    write_LAMMPS_bonded_label_v2(dat_out,num_atom,num_bond,num_angle,num_dihedral,num_atom_type,num_bond_type,num_angle_type,num_dihedral_type, \
            xlo,xhi,ylo,yhi,zlo,zhi,mass,labels,atom,bond_label,bond,angle_label,angle,dihedral_label,dihedral)