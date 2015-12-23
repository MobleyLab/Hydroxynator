#Hydroxynator

Implements the new GAFF hydoxyl parameterization of Fennell, Wymer, and Mobley (2014), which involves scaling partial charges on hydroxyl and some surrounding atoms, and new LJ parameters for hydroxyl oxygens.

Written by Caitlin Bannan, modeled after hydroxynator.py by David Mobley and hydroxynator.pl by Chris Fennell. Updated using ParmEd tools to read in topology files. It has only been tested with GROMACS topology files, but the ParmEd methods used should be able to handle other topology file types.

####Modules needed:
    * NumPy
    * ParmEd version 2.0.4 or later
        * available at http://github.com/ParmEd/ParmEd 

This can be implemented as a module in python or ran from the command line. 

#### From the Command Line
Here are the options for running from the command line: 

    Usage:  Converts sigma, epsilon, and charge for molecules with hydroxyl group
            in topology file to dielectric corrected values.
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
 
#### Import as a python module
Here are the methods available if loaded as a module:


    hydroxynate
        Parses a topology file using ParmEd tools
        changes any molecules with hydroxyl groups. 
        Outputs a topology file with the changes 
        input:
            topfile = string, input file that can be read with ParmEd tools
            outtop = string, output topology file to be created
                     if not provided it will write over the topfile
            sigmaScale = float, LJ parameter, default = 3.21990 Angstroms
            epsilonScale = float, LJ parameter, default = 0.20207 kcal/mol
            chargeScale = float, scaling fact for charge, default = 1.20905
            hydroxyl_o = string, atom type for hyroxyl oxygen default = 'oh'
            hydroxyl_h = string, atom type for hydroxyl hydrogen default = 'ho'
            charge_tol = float, tolerance for change in final-initial charge
                         default = 0.00001
        output:
            outputSys = parmed system of molecules with changes for all hydroxyl groups

    changeMolecule: 
        Identifies hydroxyl groups
        changes sigma and epsilon values on the hydroxyl oxygen 
        Scales charges near the hydroxyl groups 
        neutralizes molecule so the total charge is unchanged. 

    getTotalCharge:
        Calculates total charge on a molecule or system of molecules from parmed

    findHydroxylsAlphaCarbons:
        Finds the hydroxyl oxygens, hydrogens, and alpha atoms in a parmed molecule

    scaleAndNeutralize:
        Scales all hydroxyl oxygens, hydroxyl hydrogens, alpha carbons (heavy atoms)
        Looks at neighbors on alpha atom 
        scales hydrogens and uses others to neutralize change in net charge
        Changes sigma and epsilon values on hydroxyl oxygens


