#!/bin/env python

import os
from math import *
# Read in topology files and make changes
import parmed        
import sys
import numpy as np

"""Implements the new GAFF hydoxyl parameterization of Fennell, Wymer, and Mobley (2014), which involves scaling partial charges on hydroxyl and some surrounding atoms, and new LJ parameters for hydroxyl oxygens.

Written by David Mobley, modeled after hydroxynator.pl by Chris Fennell. This means the current implementation is relatively un-Pythonic, but it works.

Change Log:
- First version, 11/20/13
- 11/21/13: Fixed bug with tolerance for hydrogen bond lengths; fixed output file name; added option for output topology file name to have different name from input topology file.

- 9/11/2015: Complete rewrite using ParmEd tools started
"""
# Default Constants
sigmaScale = 3.21990        # Changed to match parmed units (A)
epsilonScale = 0.20207      # Changed to match parmed units (kcal/mol)
chargeScale = 1.20905
hydroxyl_o = 'oh'
hydroxyl_h = 'ho'
tol_bond_h = 0.12
charge_tol = 0.00001

#Configure input options
from optparse import OptionParser
parser = OptionParser(usage = "Converts sigma, epsilon, and charge OH values in a GROMACS topology to dielectric corrected values.\nUsage: [-options] [topology file name] \n\n", epilog = "Note: Assumes hydroxyl oxygens and hydrogens follow standard GAFF naming ('%s' and '%s' respectively; if you have hydroxyls with other atom names you will need to adjust the source code.)." % (hydroxyl_o, hydroxyl_h))

# Set options 
parser.add_option('-e', 
        help='OH epsilon conversion value, if other than standard. Default: %.5g' % epsilonScale, 
        default = epsilonScale, 
        type = "float", 
        dest = 'epsilonScale')

parser.add_option('-q', 
        help='OH environment charge scaling fraction, if other than standard. Default: %.5f' % chargeScale, 
        default = chargeScale, 
        type = "float", 
        dest = 'chargeScale')

parser.add_option('-s', 
        help='OH sigma conversion value, if other than standard. Default: %.5g' % sigmaScale, 
        default = sigmaScale, 
        type = "float", 
        dest = 'sigmaScale')

parser.add_option('-o', 
        help='Output topology file name. Default: Edit input topology file. If specified, instead creates new output topology file.', 
        dest = 'outtop', 
        type = "string" )

# Load Options
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
    parser.error('ERROR: "%s" is not a topology file I can find. Please enter the name of a valid topology file.' % topfile )


def getTotalCharge(system):
    charge = 0
    for a in system.atoms:
        charge += a.charge
    return charge

def findHydroxylsAlphaCarbons(system):
    # Find 'oh' adjust sigma, epsilon, and charges
    # Save alpha 'heavy atom' (is it possible this isn't carbon?)
    # Find 'ho' adjust charge
    alpha_carbons = []
    oxygen = 0
    hydrogen = 0
    for a in system.atoms:
        if a.type == hydroxyl_h:
            c = a.charge
            a.charge = c * chargeScale
            hydrogen += 1
        elif a.type == hydroxyl_o:
            oxygen += 1
            c = a.charge
            a.charge = c * chargeScale
            a.atom_type.sigma = sigmaScale
            a.atom_type.epsilon = epsilonScale
            if len(a.bond_partners) != 2:
                print "ERROR: hydroxyl oxygen (%s) has the wrong number of bonded neighbors check your topology file!" % str(a)
                sys.exit()
            for neighbor in a.bond_partners:
                if (not neighbor.type == hydroxyl_h) and (not neighbor in alpha_carbons):
                    alpha_carbons.append(neighbor)

    # Check for 'ho' if there is a 'oh'
    if oxygen != hydrogen:
        print "ERROR: There are a different number of hydroxyl hydrogens and hydroxyl oxygens; check your topology file!"
        sys.exit
    return system, alpha_carbons, oxygen, hydrogen

def neutralize(system, alpha_carbons):
    # Scale alpha carbon/heavy atom and attached hydrogens. 
    # Make a neutralize list
    neutralize = []
    scales = 0
    junction = []
    for a in alpha_carbons:
        c = a.charge
        a.charge = c * chargeScale
        scales += 1
        for n in a.bond_partners:
            if n.type[0] == 'h':
                n.charge *= chargeScale
                scales += 1
            elif n.type == hydroxyl_o:
                continue
                # Skip the hydroxyl oxygen, it has already been scaled
            elif n in alpha_carbons:
                if not n in junction:
                    junction.append(n)
                    print "\tWARNING: atom %i is at a juction between hydroxyl charge scaling groups. \n\t\tThis atom will be fully scaled rather than used as a neutralization site" % (n.idx+1)
            else:
                neutralize.append(n)
    # If there are no neutralizing atoms return the system as is
    if len(neutralize) == 0:
        print "WARNING: There were no neutralizing atoms found"
        return system, len(neutralize), scales
    # Now neutralize these atoms
    totalCharge = getTotalCharge(system)
    neutralFactor = totalCharge / float(len(neutralize))
    for a in neutralize:
        c = a.charge
        a.charge = c - neutralFactor
   
    return system, len(neutralize), scales

systems = parmed.load_file(topfile)
components = systems.split()
molecules = []
numbers = []
names = []
for c in components:
    molecules.append(c[0])
    numbers.append(c[1])
    names.append(c[0].residues[0].name)

print "Found %s molecule(s)" %str(len(components))

for i, molecule in enumerate(molecules):
    print "molecule", i+1
    molecule, alpha_carbons, num_hydroxyl_o, num_hydroxyl_h = findHydroxylsAlphaCarbons(molecule)
    print "\tFound %s hydroxyl groups" % (str(num_hydroxyl_o))
    if num_hydroxyl_o == 0:
        print "No hydroxyl groups found in this molecule"
        continue
    molecule, num_neutral, num_scale = neutralize(molecule, alpha_carbons)
    print "\t%s atoms were fully scaled" % str(num_hydroxyl_o + num_hydroxyl_h + num_scale)
    print "\tUsed %s atom(s) to neutralize the charge on this molecule" %str(num_neutral)
    totalCharge = getTotalCharge(molecule)
    print "\tThe total charge is %.2E" % totalCharge
    if np.abs(totalCharge) > charge_tol:
        print "\tWARNING: After scaling, the net charge on the molecule is not within the tolerance (%.2E). If you want a neutral molecule, redistribute this charge manually." % charge_tol
    molecules[i] = molecule
    print

outputSys = molecules[0] * len(numbers[0])
for idx in range(1, len(molecules)):
    outputSys += molecules[idx] * len(numbers[idx])
outputSys.write(outtop)
