import os
from math import *

"""Implements the new GAFF hydoxyl parameterization of Fennell, Wymer, and Mobley (2014), which involves scaling partial charges on hydroxyl and some surrounding atoms, and new LJ parameters for hydroxyl oxygens.

Written by David Mobley, modeled after hydroxynator.pl by Chris Fennell. This means the current implementation is relatively un-Pythonic, but it works.

Change Log:
- First version, 11/20/13
- 11/21/13: Fixed bug with tolerance for hydrogen bond lengths; fixed output file name; added option for output topology file name to have different name from input topology file.
"""

sigmaScale = 3.21990e-01
epsilonScale = 8.45476e-01
chargeScale = 1.20905
hydroxyl_o = 'oh'
hydroxyl_h = 'ho'
tol_bond_h = 0.12
charge_tol = 0.00001

#Configure input options
from optparse import OptionParser
parser = OptionParser(usage = "Converts sigma, epsilon, and charge OH values in a GROMACS topology to dielectric corrected values.\nUsage: [-options] [topology file name] \n\n", epilog = "Note: Assumes hydroxyl oxygens and hydrogens follow standard GAFF naming ('%s' and '%s' respectively; if you have hydroxyls with other atom names you will need to adjust the source code.)." % (hydroxyl_o, hydroxyl_h))
#Check on -h
parser.add_option('-e', help='OH epsilon conversion value, if other than standard. Default: %.5g' % epsilonScale, default = epsilonScale, type = "float", dest = 'epsilonScale')
parser.add_option('-q', help='OH environment charge scaling fraction, if other than standard. Default: %.5f' % chargeScale, default = chargeScale, type = "float", dest = 'chargeScale')
parser.add_option('-s', help='OH sigma conversion value, if other than standard. Default: %.5g' % sigmaScale, default = sigmaScale, type = "float", dest = 'sigmaScale')
parser.add_option('-o', help='Output topology file name. Default: Edit input topology file. If specified, instead creates new output topology file.', dest = 'outtop', type = "string" )

(options, args) = parser.parse_args()
epsilonScale = options.epsilonScale
sigmaScale = options.sigmaScale
chargeScale = options.chargeScale
topfile = args[0]
if not options.outtop:
    outtop = topfile
else:
    outtop = options.outtop

if not os.path.isfile( topfile):
    parser.error('"%s" is not a topology file I can find. Please enter the name of a valid topology file.' % topfile )


#Now process, basically following Fennell algorithm from hydroxynator.pl
file = open(topfile, 'r')
text = file.readlines()
file.close()
types_sec = False
atoms_sec = False
bonds_sec = False
stop_output = False
outtext = []
atom_columns = {}
bond_i = []
bond_j = []
bond_lengths = []
for line in text:
    if line[0]==';': #Skip lines beginning with comments
        outtext.append(line)
        continue
    if '[ atomtypes ]' in line or '[atomtypes]' in line:
        types_sec = True
    elif '[ atoms ]' in line or '[atoms]' in line:
        atoms_sec = True
        stop_output = True
    elif '[ bonds ]' in line or '[bonds]' in line:
        bonds_sec = True
        atoms_sec = False
    elif '[ pairs ]' in line or '[pairs]' in line:
        bonds_sec = False
    elif '[ moleculetype ]' in line or '[moleculetype]' in line:
        types_sec = False
    elif '[ constraints ]' in line or '[constraints]' in line:
        bonds_sec = True
 #Why is constraints treated like a bonds section?
        atoms_sec = False 


    tmp = line.split()
    if len(tmp)>1 and tmp[0][0] <> ';' and tmp[0][0] <> '[':
        if types_sec and sigmaScale <> 1. and epsilonScale <> 1.: #If we're scaling LJ and we're in the atom types section
            if len(tmp) >= 6: #Update values
                if tmp[0] == hydroxyl_o: #If we have a hydroxyl oxygen
                    tmp[5] = sigmaScale
                    tmp[6] = epsilonScale
                #Re-build line
                line = "%-5s%11s%12.4f%8.4f%3s%14.5e%13.5e\n" % (tmp[0], tmp[1], float(tmp[2]), float(tmp[3]), tmp[4], float(tmp[5]), float(tmp[6])) 
        elif atoms_sec:
            #Save the atoms section for processing    
            tmp = line.split()
            for (idx, elem) in enumerate(tmp):
                if not atom_columns.has_key(idx): atom_columns[idx] = []
                if idx==6: #If the charge, go ahead and make it a float
                    atom_columns[idx].append( float(elem) )
                else:
                    atom_columns[idx].append( elem )
        elif bonds_sec:
            bond_i.append( tmp[0] )
            bond_j.append( tmp[1] )
            bond_lengths.append( float(tmp[3]) )
        #elif not stop_output:
        #    outtext.append(line)
    #elif not stop_output:
    #    outtext.append(line)
    outtext.append(line)

#Find hydroxyl O and H atoms, track their atom numbers
n_atom_lines = len(atom_columns[0] ) #Number of atom lines
full_scale = [] #Atoms with charges which will be fully scaled
hydroxyl_o_numbers = [] #Hydroxyl oxgen indices
hydroxyl_h_numbers = []
for n in range(n_atom_lines):
    atnr = atom_columns[0][n]
    if atom_columns[1][n] == hydroxyl_o:
        full_scale.append( atnr )
        hydroxyl_o_numbers.append( atnr )
    elif atom_columns[1][n] == hydroxyl_h:
        hydroxyl_h_numbers.append( atnr )
        full_scale.append( atnr )


hydroxyl_count = len(hydroxyl_o_numbers)
print "\n%s hydroxyl moieties identified...\n\n" % hydroxyl_count
if not len(hydroxyl_h_numbers) == len(hydroxyl_o_numbers):
    print "WARNING: Unequal numbers of hydroxyl oxygens and hydrogens found."

#Now identify hydroxyl bonded carbon atoms
other_bonded = []
for i in range(hydroxyl_count): #Loop over all oxygens
    current_o = hydroxyl_o_numbers[i]
    for j in range( len(bond_i)): #Find bond involving this oxygen
        if bond_i[j] == current_o:
            is_hydrogen = False
            if bond_j[j] in hydroxyl_h_numbers:
                is_hydrogen = True
            if not is_hydrogen:
                other_bonded.append( bond_j[j] )
        #Similarly check bond_j to see if one of them is this oxygen
        if bond_j[j] == current_o:
            is_hydrogen = False
            if bond_i[j] in hydroxyl_h_numbers:
                is_hydrogen = True
            if not is_hydrogen:
                other_bonded.append(bond_i[j] )

#Fill in the full charge scaling array
neutralize_atoms = []
for i in range(len(other_bonded)):
    full_scale.append(other_bonded[i])

#Loop over the other bondeds to identify hydrogens
for i in range(len(other_bonded)):
    for j in range(len(bond_i)):
        if bond_i[j] == other_bonded[i]:
            if bond_lengths[j] < tol_bond_h:
                full_scale.append( bond_j[j] )
            else:
                is_oxygen = False 
                if bond_j[j] in hydroxyl_o_numbers:
                    is_oxygen = True
                if not is_oxygen:
                    neutralize_atoms.append( bond_j[j] )
        if bond_j[j] == other_bonded[i]:
            if bond_lengths[j] < tol_bond_h:
                full_scale.append( bond_i[j] )
            else:
                is_oxygen = False
                if bond_i[j] in hydroxyl_o_numbers:
                    is_oxygen = True
                if not is_oxygen:
                    neutralize_atoms.append( bond_i[j] )

#Sort arrays
full_scale.sort()
neutralize_atoms.sort()

#Scan for conflicts in the full_scale and neutralize_atoms lists
dont_neutralize = []
if len(full_scale) >= len(neutralize_atoms):
    for i in range(len(full_scale)):
        if full_scale[i] in neutralize_atoms:
            dont_neutralize.append( full_scale[i] )

for i in range(len(dont_neutralize)):
    print "WARNING: atom %s is at a junction between hydroxyl charge scaling groups.\n\t This atom will be fully scaled rather than used as a neutralization site.\n" % dont_neutralize[i]


#Make a new neutralization list which removes these atoms from the neutralization list
if len(dont_neutralize)> 0:
    neutralize_atoms = [ n for n in neutralize_atoms if not n in dont_neutralize ]

neutralize_atoms.sort()

#Scale charges
old_charge = 0.;
new_charge = 0.;
for i in range(len(full_scale)):
    n = int(full_scale[i])-1
    old_charge += atom_columns[6][n]
    atom_columns[6][n] *= chargeScale
    new_charge += atom_columns[6][n]

#Distribute remaining charge to neutralization atoms
charge_diff = new_charge - old_charge
if len(neutralize_atoms)>0:
    charge_diff /= float(len(neutralize_atoms))

    for i in range(len(neutralize_atoms)):
        n = int(neutralize_atoms[i])-1
        atom_columns[6][n] -= charge_diff 
else: #If no neutralization atoms, subtract charge equally from all
    charge_diff /= float(len(atom_columns[6]))
    for i in range(len(atom_columns[6])):
        atom_columns[6][i] -= charge_diff

#Now check if the molecule is neutral and warn if not
total_charge = 0.
for i in range(len(atom_columns[6])):
    total_charge += atom_columns[6][i] 

if abs(total_charge) >= charge_tol:
    print "WARNING: After scaling, the molecule had a net charge of %s.\n\t If you want a neutral molecule, redistribute this charge manually." % total_charge


#Now modify the atoms section in our stored topology file
#Find where atoms section starts
idx = 0
line = outtext[idx]
while '[ atoms ]' not in line and '[atoms]' not in line:
    idx+=1
    line = outtext[idx]
atomstart = idx+1
#Find where atoms section ends
line = outtext[idx]
while '[ bonds ]' not in line and '[bonds]' not in line:
    idx+=1
    line = outtext[idx]
atomend = idx

outct = 0
for idx in range(atomstart, atomend):
    tmp = outtext[idx].split()
    if len(tmp)>1:
        if not ';' in tmp[0] and not '[' in tmp[0]:
            outtext[idx] = "%6d%11s%7d%8s%6s%7d%11.7f%11.6f\n" % ( int(atom_columns[0][outct]), atom_columns[1][outct], int(atom_columns[2][outct]), atom_columns[3][outct], atom_columns[4][outct], int(atom_columns[5][outct]), round(atom_columns[6][outct],7), float(atom_columns[7][outct]) )
            outct+=1
    
file = open(outtop, 'w')
file.writelines(outtext)
file.close()
