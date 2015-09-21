# Hydroxynator

Implements the new GAFF hydoxyl parameterization of Fennell, Wymer, and Mobley (2014), which involves scaling partial charges on hydroxyl and some surrounding atoms, and new LJ parameters for hydroxyl oxygens.

Written by David Mobley, modeled after hydroxynator.pl by Chris Fennell. This means the current implementation is relatively un-Pythonic, but it works.

Change Log:
- First version, 11/20/13
- 11/21/13: Fixed bug with tolerance for hydrogen bond lengths; fixed output file name; added option for output topology file name to have different name from input topology file.

- 9/11/2015: Complete rewrite using ParmEd tools started

