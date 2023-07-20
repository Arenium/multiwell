Gauss2Multi program

In order to convert Gaussian or Aces2 (German version: www.aces2.de) output files into multiwell input files you need:
1) Gaussian output files with frequency calculation.
2) gauss2multi.cfg file. It must be in the same directory of gaussian output files.
 You can create and edit by yourself or use the step-by-step script to create it:
 The file contains:
1 line - Energy unit: KCAL , KJOU or CM-1  (use capital letter)
2 - number of temperatures
3 - list of temperatures separated by spaces
4 - Pressure unit: BAR , TOR, ATM , MCC (use capital letter)
5 - number of pressures
6 - list of pressures separated by spaces
7 - Egrain, imax1,  Isize, Emax2      (for more information see the multiwell documentation)
8 - increasing number,  names of gaussian output files (.log or .out), type of structure: WELL, PROD or TS (capital letter)
... 

The length of gaussian output file names hasn't to be more than 10 characters.

 
to run gauss2multi (or aces2multi)  programm, type:
(in the same directory of  gaussian/aces2 files)

'gauss2multi' or 'aces2multi'
In this case, the difference energies will be calculated with respect to the first well in gauss2multi.cfg file.

or 

'gauss2multi Gaussian_Name_File" to process a single file


For each gaussian file will be create:
name.coords        | mominert input  file
name.coords.out    | mominert output file
name.vibs          | densum   input  file
name.dens          | densum   output file
name.therm         | thermo   input  file
name.therm.out     | thermo   output file
multiwell.dat      | sample of multiwell input file

IMPORTANT:
The conversion to multiwell files cannot be not completely automated:
all these files could require manual changes! 
Be careful.


gauss2multi.cfg sample:
-------------------------------------
KCAL
 3
 298  398  498
 TOR
 2
 760  1200
 10.    400     500     50000.
 1  B1.log   WELL
 2  B2.log  WELL
 3  TS-B1-B2.log  TS
 4  C.log  PROD
 5  Product.log  PROD
---------------------------------------


ACES2MULTI NON RICONOSCE ESTERNAL SYMMETRY NUMBER !!!!!!!!!!!!!!!!!!!!!!

