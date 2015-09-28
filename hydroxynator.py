#!/bin/env python

import os
from math import *
# Read in topology files and make changes
import parmed        
import sys
import numpy as np

"""Implements the new GAFF hydoxyl parameterization of Fennell, Wymer, and Mobley (2014), which involves scaling partial charges on hydroxyl and some surrounding atoms, and new LJ parameters for hydroxyl oxygens.

Written by Caitlin Bannan, modeled after hydroxynator.py by David Mobley and hydroxynator.pl by Chris Fennell. Updated using ParmEd tools to read in topology files. 

Modules needed:
    NumPy
    ParmEd - available at http://github.com/ParmEd/ParmEd 

Change Log:
- First version, 11/20/13
- 11/21/13: Fixed bug with tolerance for hydrogen bond lengths; fixed output file name; added option for output topology file name to have different name from input topology file.

- 9/11/2015: Complete rewrite by Caitlin Bannan using ParmEd tools started 
- 9/22/2015: Added to MobleyLab GitHub repository where changes will be tracked. 
"""

def getTotalCharge(system):
    """
    Calculates total charge on a molecule or system of molecules from parmed

    input: 
        system =  could be a single molecule or system of molecules from parmed

    output:
        charge = float is the net charge on the system
    """
    charge = 0
    for a in system.atoms:
        charge += a.charge
    return charge

def findHydroxylsAlphaCarbons(system,
        hydroxyl_o = 'oh',
        hydroxyl_h = 'ho'):

    """
    Finds the hydroxyl oxygens, hydrogens, and alpha carbons (or heavy atoms) in a parmed molecule

    input:
        system = parmed molecule
        hydroxyl_o = string, atom type for hyroxyl oxygen, default = 'oh' from Amber
        hydroxyl_h = string, atom type for hydroxyl hydrogen, default = 'ho' from Amber
    
    output:
        oxygens = list of atoms in the system that have the type hydroxyl oxygen
        hydrogens = list of atoms that have the type hydroxyl hydrogens
        alpha_carbons = list of atoms that are in the alpha position to a hydroxyl group
    
    Raises an error if the number of hydroxyl oxygens and hydrogens are not equal. 
    Prints a warning if there are multiple hydroxyls on one carbon or if the alpha atom is a non-carbon. These cases are not well documented. Prints a warning if a found hydroxyl oxygen is a part of a peroxide group, peroxides are not scaled as they are believed to behave differently from hydroxyl groups.
    """
    alpha_carbons = []  # save alpha atoms
    oxygens = []        # save hydroxyl oxygens
    hydrogens = []      # save hydroxyl hydrogens

    # Parse through atoms in the system
    for a in system.atoms:

        # If hydroxyl group is found investigate
        if a.type == hydroxyl_o:
            neighbors = a.bond_partners
            partner_types = [n.type for n in neighbors]

            # hydroxyl oxygen should always have 2 neighbors
            if len(neighbors) != 2:
                oxygenBondError = Exception("ERROR: hydroxyl oxygen (%s) has the wrong number of bonded neighbors check your topology file!" % str(a))
                raise oxygenBondError
            
            # hydroxyl oxygens should always be bound to a hydroxyl hydrogen
            if not hydroxyl_h in partner_types:
                noHydroxylHydrogen = Exception("ERROR: One of the hydroxyl oxygens is not bound to a hydroxyl hydrogen. Please check your topology file")
                raise noHydroxylHydrogen

            # Skip it if the hydroxyl oxygen is a part of a peroxide group
            if 'os' in [n.type for n in a.bond_partners]:
                print "WARNING: peroxide functional group found. No scaling or LJ parameter adjustment done for peroxide groups."
                continue

            # If it passed all the checks add oxygen to hydroxyl oxygen list
            oxygens.append(a)    
            
            for n in neighbors:
                # If it's the hydroxyl hydrogen, add to hydrogen list
                if n.type == hydroxyl_h:
                    hydrogens.append(n)
                
                # Otherwise check it and add to alpha_carbons list
                else:
                    # diols on single atom have not been documented
                    if n in alpha_carbons:
                        print "WARNING: diols with two hydroxyl groups on the same carbon have not been well documented. For now, alpha carbons (or other atoms) will be scaled for each hydroxyl group attached to it."
                    
                    # Non-carbon alpha atoms have not been documented
                    elif n.type[0] != 'c':
                        print "WARNING: hydroxyl groups attached to non-carbon alpha atoms has not been well documented. For now, these will be treated the same as if the alpha atom was a carbon unless the hydroxyl is a part of a peroxide group in which case no scaling will occur."

                    # add neighbor to alpha carbon list
                    alpha_carbons.append(n)

    return oxygens, hydrogens, alpha_carbons

def scaleAndNeutralize(system, oxygens, hydrogens, alpha_carbons,
        sigmaScale = 3.21990,
        epsilonScale = 0.20207,
        chargeScale = 1.20905,
        hydroxyl_o = 'oh',
        hydroxyl_h = 'ho'):
    """
    Scales all hydroxyl oxygens, hydrogens, alpha carbons (heavy atoms), and hydrogens attached to alpha atoms. 
    Changes sigma and epsilon values on hydroxyl oxygens
    Finds heavy atoms attached to alpha atoms that will be used to neutralize any excess charge so that the total charge on the molecule is not changed

    input:
        system = parmed molecule
        oxygens = list of parmed atoms that are hydroxyl oxygens
        hydrogens = list of parmed atoms that are hydroxyl hydrogens
        alpha_carbons = list of parmed atoms that are in alpha positions to hydroxyl groups
        sigmaScale = float, LJ parameter, default = 3.21990 Angstroms
        epsilonScale = float, LJ parameter, default = 0.20207 kcal/mol
        chargeScale = float, amount the scaled atoms are scaled by, default = 1.20905
        hydroxyl_o = string, atom type for hyroxyl oxygen, default = 'oh' from Amber
        hydroxyl_h = string, atom type for hydroxyl hydrogen, default = 'ho' from Amber

    output:
        system = parmed molecule system with the changes in charge for all scaled and neutralizing atoms and the sigma and epsilon values changed for hydroxyl oxygen
        len(neutralize) = integer, number of atoms used to neutralize the system
        scales = integer, number of times an atom was scaled (this could be multiple if alpha atom has multiple hydroxyl groups)
    """

    initial_charge = getTotalCharge(system)
    scales = 0 # Used to track how many total atoms were scaled

    # Scale all of the hydroxyl oxygens and hydrogens
    for a in (oxygens + hydrogens):
        a.charge *= chargeScale
        scales += 1
        # If oxygen change epsilon and sigma
        if a.type == hydroxyl_o:
            a.atom_type.sigma = sigmaScale
            a.atom_type.epsilon = epsilonScale
     
    # Scale alpha carbon/heavy atom and attached hydrogens. 
    # Make a neutralize list
    neutralize = []
    junction = []
    for a in alpha_carbons:
        a.charge *= chargeScale
        
        # Look at neighbors and scale or add to neutralizing list
        for n in a.bond_partners:
            if n.type[0] == 'h':
                n.charge *= chargeScale
                scales += 1
            
            # Skip hydroxyl oxygen, it was scaled above
            elif n.type == hydroxyl_o:
                continue

            # If neighbor is also an alpha carbon, don't scale it
            elif n in alpha_carbons:
                if not n in junction:
                    junction.append(n)
                    print "\tWARNING: atom %i is at a juction between hydroxyl charge scaling groups. \n\t\tThis atom will be fully scaled rather than used as a neutralization site" % (n.idx+1)
            
            # If it is not already in neutralizing list add to list
            elif not n in neutralize:
                neutralize.append(n)

    # If there are no neutralizing atoms return the system as is
    if len(neutralize) == 0:
        return system, len(neutralize), scales

    # Now neutralize these atoms
    charge_off = getTotalCharge(system) - initial_charge
    neutralFactor = charge_off / float(len(neutralize))
    for a in neutralize:
        a.charge -= neutralFactor
   
    return system, len(neutralize), scales

def changeMolecule(molecule,
        sigmaScale = 3.21990,
        epsilonScale = 0.20207,
        chargeScale = 1.20905,
        hydroxyl_o = 'oh',
        hydroxyl_h = 'ho',
        charge_tol = 0.00001):
    """
    Identifies hydroxyl groups, if found, changes sigma and epsilon values on the hydroxyl oxygen and scales charges near the hydroxyl groups and neutralizes molecule so the total charge is unchanged. 

    input:
        molecule = parmed molecule to be edited
        sigmaScale = float, LJ parameter, default = 3.21990 Angstroms
        epsilonScale = float, LJ parameter, default = 0.20207 kcal/mol
        chargeScale = float, amount the scaled atoms are scaled by, default = 1.20905
        hydroxyl_o = string, atom type for hyroxyl oxygen, default = 'oh' from Amber
        hydroxyl_h = string, atom type for hydroxyl hydrogen, default = 'ho' from Amber
        charge_tol = float, warning if the final charge is not within this tolerance from the original, default = 0.00001

    output:
        molecule = parmed moelcule with the LJ parameters changed for hydroxyl oxygens and charges scaled
    """

    initial_charge = int(getTotalCharge(molecule)) 
    print "The initial charge on this molecule is %.2E" % initial_charge
    oxygens, hydrogens, alpha_carbons = findHydroxylsAlphaCarbons(molecule, hydroxyl_o, hydroxyl_h)
    
    # If no hydroxyl group found, return with no changes made
    if len(oxygens) == 0:
        print "No hydroxyl groups found in this molecule, no changes were made"
        return molecule
    else:
        print "\tFound %i hydroxyl groups" % (len(oxygens))

    molecule, num_neutral, num_scale = scaleAndNeutralize(molecule, oxygens, hydrogens, alpha_carbons, sigmaScale, epsilonScale, chargeScale, hydroxyl_o, hydroxyl_h)

    print "\t%i number of times an atom was fully scaled" % num_scale
    print "\tUsed %i atom(s) to neutralize the charge on this molecule" % num_neutral
    totalCharge = getTotalCharge(molecule)
    print "\tThe total charge is %.2E" % totalCharge
    if np.abs(totalCharge-initial_charge) > charge_tol:
        print "\tWARNING: After scaling, the net charge on the molecule is not equal to the initial charge within the tolerance (%.2E). If you want the molecule to have the intial charge (%i), redistribute this charge manually." % (charge_tol, initial_charge)
    return molecule

def hydroxynate(topfile,
        outtop = None,
        sigmaScale = 3.21990,
        epsilonScale = 0.20207,
        chargeScale = 1.20905,
        hydroxyl_o = 'oh',
        hydroxyl_h = 'ho',
        charge_tol = 0.00001):
    """
    Parses a topology file using ParmEd tools, changes any molecules with hydroxyl groups. Outputs a topology file with the changes 

    input:
        topfile = string, input file that can be read with ParmEd tools
        outtpu = string, output topology file to be created, if not provided it will write over the topfile
        sigmaScale = float, LJ parameter, default = 3.21990 Angstroms
        epsilonScale = float, LJ parameter, default = 0.20207 kcal/mol
        chargeScale = float, amount the scaled atoms are scaled by, default = 1.20905
        hydroxyl_o = string, atom type for hyroxyl oxygen, default = 'oh' from Amber
        hydroxyl_h = string, atom type for hydroxyl hydrogen, default = 'ho' from Amber
        charge_tol = float, warning if the final charge is not within this tolerance from the original, default = 0.00001

    output:
        outputSys = parmed system of molecules with changes for all hydroxyl groups and no change in net charge (within tolerance)
    """

    if outtop == None:
        outtop = topfile

    if outtop.split('.')[1] != topfile.split('.')[1]:
        wrongOutputFileType = Exception('ERROR: input and output files must both be the same file type. Please change your output file extension to match the input file.')
        raise wrongOutputFileType
    
    systems = parmed.load_file(topfile)
    components = systems.split()
    molecules = []
    numbers = []
    for c in components:
        molecules.append(c[0])
        numbers.append(len(c[1]))
        
    print "Found %s molecule(s)" %str(len(components))

    for i, molecule in enumerate(molecules):
        print "molecule", i+1
        molecules[i] = changeMolecule(molecule, sigmaScale, epsilonScale, chargeScale, hydroxyl_o, hydroxyl_h, charge_tol)
        print
    outputSys = molecules[0] * numbers[0]
    for idx in range(1, len(molecules)):
        outputSys += molecules[idx] * numbers[idx]
    outputSys.write(outtop)
    return outputSys

# If being run from the command line import Option Parser and use methods above
if __name__ == '__main__':
    #Configure input options
    from optparse import OptionParser

    # Default Constants
    sigmaScale = 3.21990        # Changed to match parmed units (A)
    epsilonScale = 0.20207      # Changed to match parmed units (kcal/mol)
    chargeScale = 1.20905
    hydroxyl_o = 'oh'
    hydroxyl_h = 'ho'
    tol_bond_h = 0.12
    charge_tol = 0.00001

    parser = OptionParser(usage = "Converts sigma, epsilon, and charge OH values in a GROMACS topology to dielectric corrected values.\nUsage: [-options] [topology file name] \n\n", epilog = "Note: Assumes hydroxyl oxygens and hydrogens follow standard GAFF naming ('%s' and '%s' respectively; if you have hydroxyls with other atom names you will need to adjust the source code.)." % (hydroxyl_o, hydroxyl_h))
    # Set options 
    parser.add_option('-e', 
            help='OH epsilon conversion value, if other than standard. Default: %.5g kcal/mol' % epsilonScale, 
            default = epsilonScale, 
            type = "float", 
            dest = 'epsilonScale')

    parser.add_option('-q', 
            help='OH environment charge scaling fraction, if other than standard. Default: %.5f Angstroms' % chargeScale, 
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

    parser.add_option('-O',
            help = "Hydroxyl oxygen atom type. The Default uses Amber atom types or 'oh'",
            default = hydroxyl_o,
            type = "string",
            dest = 'hydroxyl_o')

    parser.add_option('-H',
            help = "Hydroxyl hydrogen atom type. The Default uses Amber atom types of 'ho'",
            default = hydroxyl_h,
            type = "string",
            dest = 'hydroxyl_h')

    parser.add_option('-c',
            help = 'Charge tolerance, this is how different the charge on the scaled molecules can be from the net charge on the molecule to begin with. This should always be no bigger than rounding error. Default = 0.00001',
            default = charge_tol,
            type = "float",
            dest = 'charge_tol')

    # Load Options
    (opt, args) = parser.parse_args()
    topfile = args[0]

    if not opt.outtop:
        outtop = topfile
    else:
        outtop = opt.outtop
        
    if not os.path.isfile( topfile):
        parser.error('ERROR: "%s" is not a topology file I can find. Please enter the name of a valid topology file.' % topfile )

    hydroxynate(topfile, outtop, opt.sigmaScale, opt.epsilonScale, opt.chargeScale, opt.hydroxyl_o, opt.hydroxyl_h, opt.charge_tol)
