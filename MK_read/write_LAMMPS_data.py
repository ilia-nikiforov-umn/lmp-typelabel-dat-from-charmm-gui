"""
Creatd on Sun Aug 7 2021

Write LAMMPSdata file 

Input: 
    Structuredata 

Output: 
    LAMMPSdata file 
    
@author: Moon-ki Choi, ilia Nikiforov
"""
##################################################################
# WRITE LAMMPS BONDED DATA FILE AFTER REPLACING WITH LABELS (version 2)
##################################################################
def write_LAMMPS_bonded_label_v2(filename,num_atom,num_bond,num_angle,num_dihedral,num_atom_type,num_bond_type,num_angle_type,num_dihedral_type, \
       xlo,xhi,ylo,yhi,zlo,zhi,mass,labels,atom,bond_label,bond,angle_label,angle,dihedral_label,dihedral):
    with open(filename,'w') as fout:
        fout.write('LAMMPS data with labels (Moon-ki Choi/Ilia Nikiforov) \n\n')
        # Box information 
        fout.write('   {:d} atoms\n'.format(num_atom))
        fout.write('   {:d} bonds\n'.format(num_bond))
        fout.write('   {:d} angles\n'.format(num_angle))
        fout.write('   {:d} dihedrals\n\n'.format(num_dihedral))
        fout.write('   {:d} atom types\n'.format(num_atom_type))
        fout.write('   {:d} bond types\n'.format(num_bond_type))
        fout.write('   {:d} angle types\n'.format(num_angle_type))
        fout.write('   {:d} dihedral types\n\n'.format(num_dihedral_type))
        fout.write(' {:17.15e}    {:17.15e} xlo xhi\n'.format(xlo,xhi))
        fout.write(' {:17.15e}    {:17.15e} ylo yhi\n'.format(ylo,yhi))
        fout.write(' {:17.15e}    {:17.15e} zlo zhi\n\n'.format(zlo,zhi))
        fout.write('Atom Type Labels\n\n')
        for i in range(num_atom_type):
            fout.write('  {:d} {:s}\n'.format(i+1,labels[i]))
        if num_bond_type > 0:
            fout.write('\nBond Type Labels\n\n')
            for i in range(num_bond_type):
                fout.write('  {:d} {:s}-{:s}\n'.format(i+1,bond_label[i][0],bond_label[i][1]))
        if num_angle_type > 0:
            fout.write('\nAngle Type Labels\n\n')
            for i in range(num_angle_type):
                fout.write('  {:d} {:s}-{:s}-{:s}\n'.format(i+1,angle_label[i][0],angle_label[i][1],angle_label[i][2]))
        if num_dihedral_type > 0:
            fout.write('\nDihedral Type Labels\n\n')
            for i in range(num_dihedral_type):
                fout.write('  {:d} {:s}-{:s}-{:s}-{:s}\n'.format(i+1,dihedral_label[i][0],dihedral_label[i][1],dihedral_label[i][2],dihedral_label[i][3]))
        # NOTE: Current version of the code does not write improper type labels because target system MoS2 does not have improper information
        #fout.write('\nImproper Type Labels\n\n')
        
        fout.write('\nAtoms\n\n')
        # NOTE Current version of the code does not write image flag because imag flags for target system MoS2 are just 0 0 0
        for i in range(num_atom):
            fout.write('  {:d} {:d} {:d} {:17.15e} {:17.15e} {:17.15e} {:17.15e}   0   0   0\n'.format(i+1,int(atom[i,0]),int(atom[i,1]),atom[i,2],atom[i,3],atom[i,4],atom[i,5]))    
        if num_bond_type > 0:
            fout.write('\nBonds\n\n')    
            for i in range(num_bond):
                fout.write('  {:d}  {:d} {:d} {:d}\n'.format(i+1,int(bond[i,0]),int(bond[i,1]),int(bond[i,2])))       
        if num_angle_type > 0:
            fout.write('\nAngles\n\n') 
            for i in range(num_angle):
                fout.write('  {:d} {:d} {:d} {:d} {:d}\n'.format(i+1,int(angle[i,0]),int(angle[i,1]),int(angle[i,2]),int(angle[i,3])))
        if num_dihedral_type > 0:                
            fout.write('\nDihedrals\n\n')         
            for i in range(num_dihedral):
                fout.write('  {:d}  {:d}  {:d}  {:d}  {:d}  {:d}\n'.format(i+1,int(dihedral[i,0]),int(dihedral[i,1]),int(dihedral[i,2]),int(dihedral[i,3]),int(dihedral[i,4]))) 
    fout.close    