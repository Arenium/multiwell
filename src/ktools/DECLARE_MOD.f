c
c john r. barker
c jrbarker@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c (734) 763 6239
c
c Copyright 2019
c 
c Jason A. Sonk and john R. Barker
c University of Michigan
c Ann Arbor, MI 48109
c
c Contact: jsonk@umichedu and jrbarker@umich.edu
c
c this program is free software; you can redistribute it and/or
c modify it under the terms of the gnu general public license (version 2)
c as published by the free software foundation.
c
c this program is distributed in the hope that it will be useful,
c but without any warranty; without even the implied warranty of
c merchantability or fitness for a particular purpose. see the
c gnu general public license for more details.
c
c see the 'readme' file for a copy of the gnu general public license,
c or contact:
c
c free software foundation, inc.
c 59 temple place - suite 330
c boston, ma 02111-1307, usa.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE DECLARE_MOD
!
!    Module for sharing data among subroutines in program KTOOLS.
!
      implicit  none
      save

      character(6), parameter :: aversion = '2023'
      character(8), parameter :: adate = 'Mar 2023'
      character(6), parameter :: prog = 'ktools'
      character(46), parameter ::
     &    cut = '**************INPUT DATA SUMMARY**************'
     
      integer, parameter :: datpts =   50	! max number of trial TSs
      integer, parameter :: ndofmax = 100	! max number of d.o.f. for each chemical species
      integer, parameter :: mxcoef =   20	! max number of coefficients for hind. rotor type 'hrd'
      integer, parameter :: hrlevs = 1000	! max number of hindered rotor levels
      integer, parameter :: jtop =   1001	! max J quantum number

      integer :: maxj3		! max J quantum number for fully repulsive Veff ("testj3" in subroutine find_jmax)
c
c densum parameters
c
      real (8) :: degrain1
      real (8) :: demax2
      integer :: dimax1
      integer :: disize
c
c found in read_input
c
      character (80) :: 				fileroot	! input file name
      character (80) :: 				inputfile	! input file name root
      character (80), dimension(10) :: 	names		! file names (in main.f)

      character (180), dimension(20) ::	commentline ! up to 20 comment lines after TITLE
      integer :: nlines								! number of commentlines

      integer :: binmax			! maximum number of energy bins above TS energy
      integer :: dj				! user-defined step size for J
      integer :: etopuser		! user-defined max energy above TS
      integer :: jtopuser		! user-defined max J, the rot quantum number
      integer :: nt				! user-defined number of temperatures
      integer :: ntts			! user-defined number of trial TSs
      integer :: numnames		! in main.f
      integer :: pcnt			! user-defined number of products
      integer :: rcnt			! user-defined number of reactants

      real (8) ::  de			! user-defined step size for E (i.e. energy grain)
      real (8) ::  maxtemp		! max temperature in user-defined list of temps

      character(10) ::						tstmp
      character (150), dimension(datpts) :: 		title1
      character (150), dimension(datpts) :: 		title2
      character (150), dimension(datpts) :: 		title3
      character (180) :: 					title
      character (3) :: 						sunits		! standard state
      character (6) :: 						backup
      character (3), dimension(datpts,ndofmax) :: idofl		! type of d.o.f.
      character (4) :: 						eunits		! energy unit
      character (4), dimension(datpts,2) :: 		keyword		! 'har' vs. 'obs' and rotation units
      character (4), dimension(datpts) :: 		molname		! species name
      character (4), dimension(datpts) :: 		reprod		! reac or prod
      character (49), dimension(datpts) :: 		formulal	! chemical formula
      character (9) :: 						whatdo		! key word ('savefiles')
      real (8), dimension(datpts) :: 			sym, sopt	! symmetry no. and optical isomers
      real (8), dimension(datpts) :: 			temps		! temperatures

      integer, dimension(datpts,ndofmax) ::	mode
      real (8), dimension(datpts) :: 		vimag
      real (8), dimension(datpts) :: 		vvr

      integer, dimension(datpts) ::         	ndof			! number of d.o.f. for each species
      integer, dimension(datpts,10) :: 	    	Natom			! total number of atoms in each species
      character (6), dimension(datpts,10) ::	ATYPE			! atom type (up to 10) in chemical formula
      integer, dimension(datpts) :: 	    	Nelement		! number of atoms of each element
      integer, dimension(datpts,10) ::      	gele			! electonic level degeneracy; up to 10 levels for each species
      integer, dimension(datpts,ndofmax) :: 	ngl				! vib degeneracy
      integer, dimension(datpts) :: 	    	nele			! number of electronic levels for each species
      real (8), dimension(datpts,ndofmax) :: wel				! vib frequency
      real (8), dimension(datpts,ndofmax) :: anhl			! anharmonicity
      real (8), dimension(datpts,10) :: 	elev			! electonic energy levels; up to 10 for each species
      real (8), dimension(datpts) :: 		delh			! species enthalpy
      real (8), dimension(datpts) :: delhf
      real (8), dimension(datpts) :: delhr
      real (8), dimension(datpts) :: distl
c
c	Input for hindered rotors of type 'hrd'
c
      integer, dimension(datpts,ndofmax) :: 		ncbl		! hind rot 'hrd': NUMBER of rot coefficients
      integer, dimension(datpts,ndofmax) :: 		ncvl		! hind rot 'hrd': NUMBER of pot enrg coefficients
      character (5), dimension(datpts,ndofmax) :: vhrl		! hind rot 'hrd": TYPE of pot enrg function
      integer, dimension(datpts,ndofmax) :: 		nsvl		! hind rot 'hrd': SYMMETRY of pot enrg function
      real (8), dimension(datpts,ndofmax) :: 	phavl		! hind rot 'hrd': PHASE of pot enrg function
      real (8), dimension(datpts,ndofmax,mxcoef) :: cvl			! hind rot 'hrd': pot enrg COEFFICIENTS
      character (5), dimension(ndofmax,datpts) :: bhrl		! hind rot 'hrd": TYPE of rotation function
      integer, dimension(datpts,ndofmax) :: 		nsbl		! hind rot 'hrd': SYMMETRY of rot function
      real (8), dimension(datpts,datpts) :: 		phabl		! hind rot 'hrd': PHASE of rot function
      real (8), dimension(datpts,datpts,mxcoef) ::cbl			! hind rot 'hrd': rot COEFFICIENTS
c
c part of hrd
c
      integer, dimension(datpts,ndofmax) :: 			valmax  	! no. of levels returned: nimax (in tshrlev.f) and nmax (in tghrlev.f)
      real (8), dimension(datpts,ndofmax,hrlevs) ::	evhl		! hindered rotor energy levels (zpe_calc.f)
      real (8), dimension(hrlevs) :: 				evh		! hindered rotor energy levels (canon_rates.f)
c
c found in bsort
c
      real (8), dimension(datpts) :: barr
      real (8), dimension(datpts) :: aarr
c
c found in partition_funcs
c
      real (8), dimension(4,datpts) :: ratekinf
      real (8), dimension(datpts)   :: interr
      integer, dimension(datpts)    :: minvec
      integer, dimension(datpts)    :: starvec
      integer, dimension(datpts)    :: divide
c
c found in find_jmax
c
      integer ::  maxj			! After tests, max value of J used in calcs; not necessarily equal to jtopuser 
c
c found in jloop/veff
c
      real (8), dimension(datpts) :: vef, epsil
      real (8), dimension(jtop+1) :: v1l
      real (8), dimension(jtop+1) :: vlist
      real (8), dimension(jtop+1) :: vefmax
      real (8) :: vmax
      integer ::  vmaxi
      integer ::  epsmax
c
c found in sterab or sterabj
c
      integer, dimension(jtop+1) :: nbinl			! max number of energy bins for reactant, given a value of J
c
c super_mol
c
      real (8), dimension(datpts) :: amu
      real (8), dimension(2)      :: tamu
      real (8), dimension(2)      :: ramu
c
c found in calc_rates
c
      real (8), dimension(datpts) :: cratef
      real (8), dimension(datpts) :: ucratef
      real (8), dimension(datpts) :: crater
      real (8), dimension(datpts) :: ucrater
      real (8), dimension(datpts) :: mrate
      real (8), dimension(datpts) :: umrate
      real (8), dimension(datpts) :: queues
      logical :: manymins

      END MODULE DECLARE_MOD
      
