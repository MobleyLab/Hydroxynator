Implements the new GAFF hydoxyl parameterization of Fennell, Wymer, and Mobley (2014), which involves scaling partial charges on hydroxyl and some surrounding atoms, and new LJ parameters for hydroxyl oxygens.

Written by Caitlin Bannan, modeled after hydroxynator.py by David Mobley and hydroxynator.pl by Chris Fennell. Updated using ParmEd tools to read in topology files. With this version forward

Modules needed:
    NumPy
    ParmEd - available at http://github.com/ParmEd/ParmEd 

This can be implemented as a module in python or ran from the command line. 

Here are the options for running from the command line: 

    Usage: Converts sigma, epsilon, and charge OH values in a GROMACS topology to dielectric corrected values.
    Usage: [-options] [topology file name]

    Options:
      -h, --help       show this help message and exit
      -e EPSILONSCALE  OH epsilon conversion value, if other than standard.
                       Default: 0.20207 kcal/mol
      -q CHARGESCALE   OH environment charge scaling fraction, if other than
                       standard. Default: 1.20905 
      -s SIGMASCALE    OH sigma conversion value, if other than standard. Default:
                       3.2199 Angstroms
      -o OUTTOP        Output topology file name. Default: Edit input topology
                       file. If specified, instead creates new output topology
                       file.
      -O HYDROXY_O    Hydroxyl oxygen atom type. The Default uses Amber atom
                       types or 'oh'
      -H HYDROXYL_H    Hydroxyl hydrogen atom type. The Default uses Amber atom
                       types of 'ho'
      -c CHARGE_TOL    Charge tolerance, this is how different the charge on the
                       scaled molecules can be from the net charge on the molecule
                       to begin with. This should always be no bigger than
                       rounding error. Default = 0.00001

Here are the methods available if loaded as a module: 

hydroxynate
    
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
        outputSys = parmed system of molecules with changes for all hydroxyl groups

changeMolecule: 
    
    Identifies hydroxyl groups, if found, changes sigma and epsilon values on the hydroxyl oxygen and scales charges near the hydroxyl groups and neutralizes molecule so the total charge is unchanged. 

getTotalCharge:
    
    Calculates total charge on a molecule or system of molecules from parmed

findHydroxylsAlphaCarbons:
    
    Finds the hydroxyl oxygens, hydrogens, and alpha carbons (or heavy atoms) in a parmed molecule

scaleAndNeutralize:
    
    Scales all hydroxyl oxygens, hydrogens, alpha carbons (heavy atoms), and hydrogens attached to alpha atoms. 
    Changes sigma and epsilon values on hydroxyl oxygens
    Finds heavy atoms attached to alpha atoms that will be used to neutralize any excess charge so that the total charge on the molecule is not changed


