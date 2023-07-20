!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!                                                                      !   
!  MultiWell: a code for master equation simulations.                  !   
!  Copyright (C) 2001 - 2021 John R. Barker                            !   
!                                                                      !   
!  John R. Barker                                                      !   
!  jrbarker@umich.edu                                                  !   
!  University of Michigan                                              !   
!  Ann Arbor, MI 48109-2143                                            !   
!  (734) 763 6239                                                      !   
!                                                                      !   
!  This program is free software; you can redistribute it and/or       !   
!  modify it under the terms of the GNU General Public License         !   
!  (version 2) as published by the Free Software Foundation.           !   
!                                                                      !   
!  This program is distributed in the hope that it will be useful,     !   
!  but WITHOUT ANY WARRANTY; without even the implied warranty of      !   
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the        !   
!  GNU General Public License for more details.                        !   
!                                                                      !   
!  See the 'ReadMe' file for a copy of the GNU General Public License, !   
!  or contact:                                                         !   
!                                                                      !   
!  Free Software Foundation, Inc.                                      !   
!  59 Temple Place - Suite 330                                         !   
!  Boston, MA 02111-1307, USA.                                         !   
!                                                                      !   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

      PROGRAM MultiWell

      USE declare_mod
      USE utili_mod
      USE subs_mod
      USE I_O_mod
      USE bookstep_mod

      IMPLICIT NONE

      CHARACTER(len=20) :: NameFile, NameOUT, NameSMM, NameARY,
     &     NameDIST, NameFLX, NameRAT
      INTEGER :: Nforward , M , Igo , IPP , iP 
      INTEGER :: Iwell , ichan , i, j, k, ki,  Np , jstat , kkmax , NCUM
      REAL(8) :: Test , NumTot , Tread , numbdens
      REAL(8) :: DD , Colli , Pfact , Delt, Trials
c
c     FOR USER-DEFINED FILENAME ("NameFile")
c
      CALL get_command_argument( 1, NameFile )
      NameFile=NameFile(1:len_trim(NameFile))
      IF(NameFile.eq.'') NameFile="multiwell.dat"         ! default filename
      NameOUT=NameFile(1:len_trim(NameFile)-4)//".out"
      NameSMM=NameFile(1:len_trim(NameFile)-4)//".sum"
      NameARY=NameFile(1:len_trim(NameFile)-4)//".array"
      NameDIST=NameFile(1:len_trim(NameFile)-4)//".dist"
      NameRAT=NameFile(1:len_trim(NameFile)-4)//".rate"
      NameFLX=NameFile(1:len_trim(NameFile)-4)//".flux"
c
c     I/O Unit numbers
c
C      KIN   = 2  ! Data file
C      KSTD  = 6  ! standard output (screen)
C      KSMM  = 7  ! Short summary of output file 
C      KARY  = 8  ! Summary of array variables
C      KOUT  = 9  ! General output
C      KRAT  = 10 ! average rate 'constants'
C      KFLX  = 11 ! relative reactive flux via a given path
C      KDIS  = 12 ! vib distribution

      OPEN (KIN,FILE=NameFile,STATUS='OLD')          ! Data file
      OPEN (KOUT,FILE=NameOUT,STATUS='UNKNOWN')      ! General time-dependent output
      OPEN (KSMM,FILE=NameSMM,STATUS='UNKNOWN')      ! Short summary output for fall-off calculations
      OPEN (KARY,FILE=NameARY,STATUS='UNKNOWN')      ! Summary of input, output, and intermediate array variables
      OPEN (KRAT,FILE=NameRAT,STATUS='UNKNOWN')      ! average rate 'constants'
      OPEN (KFLX,FILE=NameFLX,STATUS='UNKNOWN')      ! relative reactive flux via a given reaction path
  
      READ (KIN,99003) TITLE

      READ (KIN,*) Egrain1 , imax1 , Isize , Emax2 , IDUM

      IF ( Isize .GT. Imax ) THEN
         WRITE (KSTD,99002) Isize
         WRITE (KOUT,99002) Isize
         STOP
      ENDIF
      IF ( imax1 .GT. Isize-2 ) THEN
       WRITE (KSTD,*) '*** FATAL:  Imax1 must be ² (Isize-2) ***'
       WRITE (KOUT,*) '*** FATAL:  Imax1 must be ² (Isize-2) ***'
       STOP
      ENDIF      
      
      Emax1 = Egrain1*(imax1-1)
      Egrain2 = Emax2/(Isize-imax1-1)
c.......................................................................
c     Energy grain parameters for interpolation 'double arrays'
c
c     'Double arrays' have two sections:
c          segment 1 consists of equally spaced (Egrain1) values ranging
c               from E=0 to Emax1
c          segment 2 consists of equally spaced values from E=0 to Emax2;
c               the spacing depends on the number of array elements
c               remaining: Isize-Emax1/Egrain1-1
c
c     Egrain1     = energy grain of lower segments in 'double array' (units: cm-1)
c     Emax1     = maximum energy of 1st segment of double array (units: cm-1)
c     Emax2     = maximum energy of 2nd segment of double array (units: cm-1)
c     IDUM     = random number seed (positive, non-zero integer)
c.......................................................................

      DistStep = Emax2/(Ndist-1)                ! Energy step for binned distributions; used in book1 and output
      DistStep = INT( DistStep/Egrain1 + 1 )*Egrain1

      IF ( IDUM.EQ.0 ) IDUM = 1990339   ! Random number seed
      IDUM0 = IDUM
      IDUM = -ABS(IDUM)
      XRANST = RAN1(IDUM)    ! initialize RAN1

c.......................................................................
c     Key words for units
c    READ Punits , Eunits , Rotatunits
c
c     Punits       = pressure units key word [use upper case CHARACTER*3]
c                 = 'BAR', 'ATM', 'TOR' [also will accept 'TORR'], or 'MCC'
c     Eunits       = energy units keyword [use upper case CHARACTER*4]
c                 = 'CM-1','KCAL', or 'KJOU' (for cm-1, kcal/mole, or kJ/mole)
c     Rotatunits  = units for moment of inertia 
c                = 'AMUA', 'GMCM', 'CM-1', MHz, 'GHz'
c.......................................................................
      READ (KIN,*) DUM(1) , DUM(2) , DUM(3)
      CALL ucase ( DUM(1) )   ! convert to upper case
      CALL ucase ( DUM(2) )   ! convert to upper case
      CALL ucase ( DUM(3) )   ! convert to upper case
      CALL UnitTest ( DUM , Punits , Eunits , Rotatunits , KOUT )
c.......................................................................
c     Physical parameters
c
c     Temp     = translational temperature (units: degrees Kelvin)
c     Tvib     = vibrational temperature (units: degrees Kelvin)
c     Np     = number of pressures
c     PP     = vector of Np pressures
c.......................................................................
      READ (KIN,*) Temp , Tvib
      READ (KIN,*) Np              ! number of pressures (or number densities)
      READ (KIN,*) (PP(i),i=1,Np)
      IF ( Np.GT.50 ) THEN
         write(*,*) 'FATAL: Only 50 pressures allowed....'
         write(KOUT,*) 'FATAL: Only 50 pressures allowed....'
         STOP
      ENDIF
c.......................................................................
c
c     Convert pressures to bars
c
c.......................................................................
      IF ( Punits.EQ.atm ) THEN
         POUT = 'atm '
         Pfact = 1.01325
      ELSEIF ( Punits.EQ.tor ) THEN
         POUT = 'Torr'
         Pfact = 1.01325/760.
      ELSEIF ( Punits.EQ.bar ) THEN
         POUT = 'bar '
         Pfact = 1.0
      ELSEIF ( Punits.EQ.mcc ) THEN
         POUT = '1/cc'
         Pfact = Temp/7.2428D+21
      ELSE
         WRITE(*,*) 'FATAL: Pressure units not specified properly'
         WRITE(KOUT,*) 'FATAL: Pressure units not specified properly'
         STOP
      ENDIF
 
      DO 100 i = 1 , Np
         P(i) = PP(i)*Pfact            ! P(i) converted to bar
 100  CONTINUE
 
      READ (KIN,*) NWells , NProds

c.......................................................................
c     Reactant/Product parameters
c
c     NWells     = number of 'wells'
c     NProds     = number of entrance/exit species
c
c     IMol     = index number for molecule
c     MolName     = name (10 characters) of molecule
c     HMol     = enthalpy of formation at 0¡K
c     MolMom     = moment of inertia for 2-dim external rotor (amu Ang^2), 
c                gram*cm^2, cm-1, MHz, or GHz
c     Molsym     = external symmetry number for molecule
c     Molele     = electronic partition function for molecule
c     Molopt     = number of optical isomers for molecule
c     Iread     =1: read densities from file
c            *** filename: same name as MolName +'.dens'
c            ***    e.g. 'Benzene.dens'
c.......................................................................

      jstat = 0
      ALLOCATE( dens(Nwells,Isize) , STAT=istat)
      jstat = jstat + istat
      ALLOCATE( COLLUP(Nwells,Isize) , STAT=istat)
      jstat = jstat + istat      
      ALLOCATE( CNORM(Nwells,Isize) , STAT=istat)
      jstat = jstat + istat      
      ALLOCATE( MolName(Nwells+Nprods) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( IMol(Nwells+Nprods) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( HMol(Nwells+Nprods) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( Molsym(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( Molopt(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( MolMom(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( Molele(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( Viblo(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( Sig(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( Eps(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( Sigma(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( Epsil(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( Mass(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( klj(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( kqm(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( CollK(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( LJQM(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( Nchan(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( Var(Nwells+Nprods) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( FracMol(Nwells+Nprods) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( AveE(Nwells+Nprods) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( nsmax(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( CNORMref(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( Eref(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( iset(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( iset2(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( ITYPE(Nwells) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( DC(Nwells,8) , STAT=istat)
      jstat = jstat + istat 
      ALLOCATE( LI(Nwells,5) , STAT=istat)
      jstat = jstat + istat 

      IF ( jstat .NE. 0 ) THEN
         write(*,*) '*** FATAL array(Nwells) not allocated ***'
         write(*,*) '         STAT=',istat
         STOP
      ENDIF

      DO 200 i = 1 , NWells
         READ (KIN,*) IMol(i) , MolName(i) , HMol(i) , MolMom(i) , 
     &                Molsym(i) , Molele(i) , Molopt(i)

         IF ( Eunits.EQ.kcal ) THEN
            HMol(i) = HMol(i)*caltocm           ! Convert kcal/mole to cm-1
         ELSEIF ( Eunits.EQ.kjou ) THEN
            HMol(i) = HMol(i)*jtocm             ! Convert kJ/mole to cm-1
         ENDIF
         Nchan(i) = 0                           ! Initialize Nchan
         call rotunits( Rotatunits , MolMom(i) , 'AMUA' )     ! convert rotational units to 'AMUA'
200   CONTINUE
 
      DO 300 i = NWells + 1 , NWells + NProds
         READ (KIN,*) IMol(i) , MolName(i) , Hmol(i)
         IF ( Eunits.EQ.kcal ) THEN
            HMol(i) = HMol(i)*caltocm           ! Convert kcal/mole to cm-1
         ELSEIF ( Eunits.EQ.kjou ) THEN
            HMol(i) = HMol(i)*jtocm           ! Convert kJ/mole to cm-1
         ENDIF
300   CONTINUE
 

      CALL DensArray                            ! Calculate or read densities of states
 
      
c.......................................................................
c     Collision parameters - collider
c
c     SigM     L-J sigma (Ang) for collider
c     EpsM     L-J Eps (Kelvins) for collider
c     AmuM     Molecular weight (g/mole) of collider
c     ETKEY    Keyword to designate energy transfer treatment:
c              = OLDET     traditional treatment of energy transfer at low E
c              = NEWET     Barker's 'New Approach" to energy transfer
c                [J.R. Barker, Int. J. Chem. Kinetics, 41, 748-763 (2009)]
c.......................................................................
      READ (KIN,*) SigM , EpsM , AmuM , Amu
      
      READ (KIN, 99173)  ETKEY
      call ucase ( ETKEY)
      IF ( ETKEY.NE.NEWET .AND. ETKEY.NE.OLDET ) THEN      ! not specified by user
        BACKSPACE (KIN)
        ETKEY = xxxet
      ENDIF
c.......................................................................
c
c     klj = L-J collision rate constant (s-1)
c     kqm = collision rate constant based on the total quantum cross section
c          [Durant and Kaufman, Chem. Phys. Letters, 142, 246 (1987)]
c
c     Collision parameters
c
c     Sig     L-J sigma (Ang) for molecule
c     Eps     L-J Eps (Kelvins) for molecule
c     Amu     Molecular weight (g/mole) of molecule
c     ITYPE     selects model type in Subroutine PSTEP
c     DC     coefficients for collision model
c     LJQM     = 'LJ' for Lennard-Jones collision frequency, ALL collisions inelastic
c          = 'QM' for total cross section collision frequency
c.......................................................................
c
      DO 500 j = 1 , NWells
         READ (KIN,*) Mol , Sig(Mol) , Eps(Mol) , ITYPE(Mol) , 
     &                 (DC(Mol,i),i=1,8)
         READ (KIN,*) LJQM(Mol)
         call ucase( LJQM(Mol) )
         Sigma(Mol) = 0.5*(SigM+Sig(Mol) )
         Epsil(Mol) = SQRT(EpsM*Eps(mol) )
         Mass(Mol) = AmuM*Amu/(AmuM+Amu)
 
c         Colli = 1.0d+00/(0.636+0.246*LOG(Temp/Epsil))                ! simple collision integral [Troe, J. Chem. Phys. 66, 4758 (1977)]

         Colli = 1.16145/((Temp/Epsil(Mol))**0.14874) +                    ! better collision integral [Reid, Prausnitz, and Sherwood, The Properties
     &              0.52487*exp(-0.7732*Temp/Epsil(Mol)) +                 !   of Gases and Liquids: Their Estimation and Correlation, 3rd Edition,
     &              2.16178*exp(-2.437887*Temp/Epsil(Mol))                 !   McGraw-Hill, New York, 1977]

         klj(mol) = dpi*1.E-16*(Sigma(Mol)**2)*1.45525E+04*
     &                               SQRT( Temp/Mass(Mol) )*Colli
         
         kqm(mol) = 7.64E-12*((Epsil(Mol)*Sigma(Mol)**6)**0.4)*
     &                                      (Temp/Mass(Mol))**0.3

         IF ( LJQM(Mol).EQ.LJ ) THEN
            CollK(Mol) = klj(Mol)
         ELSE
            CollK(Mol) = kqm(Mol)
         ENDIF

 500  CONTINUE
  
         DO Mol = 1 , NWells

            CALL COLNORM(Mol,Temp)         ! Calculate COLLUP and CNORM arrays for all Wells

         END DO

      CALL write_1            ! write to output files

      IF ( ETKEY .EQ. xxxet ) ETKEY = newet 
  
      READ (KIN,*) Nforward                     ! total number of forward reactions

      IF (Nforward .GT. 0 ) THEN
        WRITE(KOUT,99074) Nforward
        
        jstat = 0
        ALLOCATE( Nreaindex(Nforward) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Nproindex(Nforward) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( WARN1(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Jrev(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Jto(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( TSsym(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( TSele(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( TSopt(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Jread(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Path(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Afac(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Eo(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Eor(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Hts(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( NCENT(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( TSmom(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Press(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( iivr(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( vivr(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( vave(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( pivr(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( tivr(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( civr(Nwells,MaxChan,3) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( V1(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( vimag(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( itun(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( TSname(Nwells,MaxChan) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( SupName(Nwells,MaxSup) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( nsup(Nwells,MaxSup) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( sto(Nwells,MaxSup) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( norder(Nwells,MaxSup) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Afc(Nwells,MaxSup) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Bsr(Nwells,MaxSup) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( dde(Nwells,MaxSup) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( ksup(Nwells,MaxSup) , STAT=istat)
        jstat = jstat + istat 

        IF ( jstat .NE. 0 ) THEN
           write(*,*) '*** FATAL: Nreaindex, etc. not allocated ***'
           write(*,*) '         STAT=',jstat
           STOP
        ENDIF       
     
        CALL rxn_input(Nforward)          ! read input parameters for forward reactions; generates params for reverse

        largest = 0
        DO i = 1 , NWells
          largest = max( largest,Nchan(i) )
        END DO        

        jstat = 0
        ALLOCATE( Rate( NWells,largest,Isize) , STAT=istat )
        jstat = jstat + istat
        ALLOCATE( KRate( NWells,largest,Isize+imax1) , STAT=istat )
        jstat = jstat + istat
        ALLOCATE( TSsum( NWells,largest,Isize) , STAT=istat )
        jstat = jstat + istat
        ALLOCATE( TRate( NWells,largest,Isize+imax1) , STAT=istat ) 
        jstat = jstat + istat
        ALLOCATE( DensTS( largest,Isize) , STAT=istat )        ! indices: Nchan,energy index); previously (NWells*largest,Isize)
        jstat = jstat + istat
        ALLOCATE( TSdensBK( largest,Isize+imax1) , STAT=istat )
        jstat = jstat + istat
        ALLOCATE( Rab(largest+nsux+10) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Flux(Nwells,largest,Ntime) , STAT=istat)
        jstat = jstat + istat 

        IF ( jstat .NE. 0 ) THEN
           write(*,*) '*** FATAL: Rate, TSsum, etc. not allocated ***'
           write(*,*) '         STAT=',jstat
           STOP
        ENDIF       

         CALL RateArray(Temp)      ! Calculate microcanonical rate constants at temperature Temp

      ENDIF  ! (Nforward .GT. 0 )
       
      EM = (Isize-imax1-2)*Egrain2           ! Top energy in double array

      DO j = 1 , NWells
         IF ( Nchan(j) .EQ. 0 ) THEN         ! New up-collision method: normalization at imax1 (if no rxn channel)
           CNORMref(j) = CNORM(j,imax1)
           Eref(j) = (imax1-1)*Egrain1
         ELSE
           CNORMref(j) = CNORM(j,imax1)
           Eref(j) = Emax1
           DD = EM                           ! Top energy in double array
           DO i = 1 , Nchan(j)
             IF ( Eor(j,i) .LT. DD .AND. Eor(j,i) .GT. Emax1 ) THEN ! New up-collision method: normalization at Eor
                DD = Eor(j,i)            ! The lowest reaction threshold that is higher than Emax1
                Eref(j) = DD
                CNORMref(j) = fxnofe( DD , j , CNORM )
             ENDIF
           END DO
         ENDIF
      END DO

c-------------------------------------------------------------------------------------
c-------------------------------------------------------------------------------------
c     Rtrials = total number of trials, expressed as real variable
c     Tspec   = time specification (CHARACTER*4)
c             = 'TIME': Tread = Tlim (max time)
c             = 'COLL': Tread = max collisions, from which Tlim is calculated 
c           = 'TLOG': Tread = Tlim (max time) for log10( time)
c     Tread   = maximum simulated time
c     KEYTYPE = KEYWORD  for type of initial distribution
c             = DELTA    for Monoenergetic at energy Einit          (Icalc=1)
c             = THERMAL  for Thermal with energy offset Einit       (Icalc=2)
c             = CHEMACT  for Chemical activation from 'product' #IR (Icalc=3)
c             = EXTERNAL for read from external file                (Icalc=4)
c     Molinit = index of initial molecule
c     IR      = origin of initial molecule, if produced via chemical activation
c     Einit   = initial energy (relative to ZPE of Molinit)
c-------------------------------------------------------------------------------------
c
      READ (KIN,*) Rtrials , Tspec , Tread , KEYTEMP , Molinit , IR , 
     &                Einit
      call ucase( Tspec )
      call ucase( KEYTEMP )
      IF ( KEYTEMP.EQ.'DELTA'   ) THEN
         Icalc = 1
        ELSEIF ( KEYTEMP.EQ.'THERMAL' ) THEN
         Icalc = 2
        ELSEIF ( KEYTEMP.EQ.'CHEMACT' ) THEN
         Icalc = 3
        ELSEIF ( KEYTEMP.EQ.'EXTERNAL' ) THEN
         Icalc = 4
        ELSE
         write(6,*)'FATAL: INITIAL ENERGY DISTRIBUTION NOT RECOGNIZED'
         write(KOUT,*)'FATAL: INITIAL ENRGY DISTRIBUTION NOT RECOGNIZED'
         STOP
      ENDIF

      IF ( Eunits.EQ.kcal ) THEN
         Einit = Einit*caltocm          ! Convert kcal/mole to cm-1
        ELSEIF ( Eunits.EQ.kjou ) THEN
         Einit = Einit*jtocm          ! Convert kJ/mole to cm-1
      ENDIF
 
!------------------------------------------------------------------------------------------
!   WRITE OUT COLLISIONAL NORMALIZATION, UP-TRANSITION PROBABILITIES

      CALL write_2            ! write to output files
!------------------------------------------------------------------------
     
      IF (Nforward .GT. 0 ) THEN
        jstat = 0
        DEALLOCATE( TSsum , STAT=istat )
        jstat = jstat + istat
        DEALLOCATE( DensTS , STAT=istat )
        jstat = jstat + istat
        DEALLOCATE( TSdensBK , STAT=istat )
        jstat = jstat + istat
        DEALLOCATE( KRate, STAT=istat )
        jstat = jstat + istat
        IF ( jstat .NE. 0 ) THEN
         write(*,*) '*** FATAL: arrays TSsum, etc. not DEallocated ***'
         write(*,*) '         STAT=',istat
         STOP
        ENDIF
      ENDIF  ! (Nforward .GT. 0 )    
 
        jstat = 0
        ALLOCATE( Emol(Nwells,Ntime), STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( EDIST(Nwells,Ndist,Mtime) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( kuni(Nwells,largest,Ntime) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Num(Nwells+Nprods,Ntime) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Mum(MaxWell,Mtime) , STAT=istat)
        jstat = jstat + istat 
        ALLOCATE( Nuni(Nwells,Ntime), STAT=istat)
        jstat = jstat + istat 

        IF ( jstat .NE. 0 ) THEN
           write(*,*) '*** FATAL: Emol, EDIST, etc. not allocated ***'
           write(*,*) '         STAT=',jstat
           STOP
        ENDIF       
     


!------------------------------------------------------------------------
     
c @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     Initialize Estart to calculate cumulative distribution
c
      write(*,*) 'Estart Initialization'
      
      Jsize = 50 + INT( 30.d+00*Etherm( Molinit,Tvib )/Egrain1 )     ! Set array length
      IF ( Jsize .GT. Imax ) Jsize = Imax

        jstat = 0
      ALLOCATE ( Pstart(Jsize), STAT=istat )
        jstat = jstat + istat
      ALLOCATE ( Nstart(Jsize), STAT=istat )
        jstat = jstat + istat
        IF ( jstat .NE. 0 ) THEN
         write(*,*) 
     &       '*** FATAL: arrays Pstart & Nstart not DEallocated ***'
         write(*,*) '         STAT=',istat
         STOP
        ENDIF

      E = Estart(0,Icalc,IR,Molinit,Einit,Tvib)                   ! invoke Estart with flag=0 to compute cumulative

      Iwell = Molinit    ! save index for the well that the reaction will be initiated
                       ! from so that the correct Eor(Iwell,1) can be printed before the
                       ! cumulative distribution function info in multiwell.array
  
c @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c     START PRESSURE STEPS & STOCHASTIC TRIALS
c
c     Igo   = flag which signals the following:
c           = 0   stop stepping for this trial
c           = 1   continue stepping for this trial
c
      write(*,*) 'Start stepping'

      DO 1200 Ip = 1 , Np               ! Pressure steps
         numbdens = P(Ip)*7.2428D+21/Temp
         CollFreq = CollK(Molinit)*numbdens

         IF ( Tspec.EQ.'TIME' ) THEN           ! Specify max time
            Tlim = Tread
            Tstep = Tlim/(Ntime-1)
         ELSEIF ( Tspec.EQ.'LOGT' ) THEN       ! Specify max time
            Tlim = Tread
            Tstep = log( Tlim / tminlogt )/(Ntime-1)  ! stepsize on log(t) scale from log(tminlogt) up to log(Tlim)
         ELSEIF ( Tspec.EQ.'COLL' ) THEN       ! Specify max number of collisions
           IF (CollFreq .GT. 0.0d+00) THEN
              Tlim = Tread/CollFreq
              Tstep = Tlim/(Ntime-1)
           ELSE
              Tlim = 100.
           ENDIF
         ELSE
           write(*,*) '** Keyword Tspec is not recognized ***'
           STOP (' *** FATAL: Tspec = TIME, COLL, or LOGT, only *** ')
         ENDIF

         IDUM0 = IDUM
 
         CALL INITIAL                   ! Initialize bookkeeping arrays
 
         WRITE (*,99164) PP(Ip)
 
         Trials = 0.0d+00
         DO WHILE ( Trials .LE. Rtrials )
            Trials = Trials + 1.0d+00
            E = Estart(1,Icalc,IR,Molinit,Einit,Tvib)
            Igo = 1
            Mol = Molinit
            Time = 0.0D+00               ! Initialize
            DO WHILE ( Igo .EQ. 1 )
               CALL STEPPER(Mol,E,numbdens,IPP,Delt)
               CALL BOOK1(Mol,E,Temp,IPP,Time,Delt,Igo)
               IF ( Time.GE.Tlim ) Igo = 0
            ENDDO
            IF ( MOD( Trials , Rtrials/20.d+00 ) .LT. 1.e-10 ) THEN
              WRITE (*,99064) Trials ,   NINT( 100*Trials/Rtrials )
            ENDIF
         END DO  ! Trials
c
c @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     Massage booked results and WRITE output for each pressure
c
c     KEYTYPE = KEYWORD  for type of initial distribution
c             = DELTA    for Monoenergetic at energy Einit          (Icalc=1)
c             = THERMAL  for Thermal with energy offset Einit       (Icalc=2)
c             = CHEMACT  for Chemical activation from 'product' #IR (Icalc=3)
c             = EXTERNAL for read from external file                (Icalc=4)

         CALL DateTime(KOUT)
         CALL DateTime(KARY)
 
         DO ki = 8 , 11
          WRITE (ki,99059) Rtrials
          IF ( Icalc.EQ.1 ) THEN                                        ! DELTA
            WRITE (ki,99065) Temp , PP(Ip) , POUT , numbdens , 
     &                    CollFreq , Eref(Molinit), 
     &                    ADJUSTR(MolName(Molinit)) , 
     &                    Einit , IDUM0 , IDUM
            ELSEIF ( Icalc.EQ.2 ) THEN                                  ! THERMAL
            WRITE (ki,99066) Tvib , Temp , PP(Ip) , POUT , numbdens , 
     &                    CollFreq , Eref(Molinit), 
     &                    ADJUSTR(MolName(Molinit)) , 
     &                         Einit , IDUM0 , IDUM
            ELSEIF ( Icalc.EQ.3 ) THEN                                  ! CHEMACT
            WRITE (ki,99067) Tvib , Temp , PP(Ip) , POUT , numbdens , 
     &                    CollFreq , Eref(Molinit), 
     &                    ADJUSTR(MolName(Molinit)) , 
     &                    ADJUSTR(MolName(IR)) , Einit, IDUM0 , IDUM
            ELSEIF ( Icalc.EQ.4 ) THEN                                  ! EXTERNAL
            WRITE (ki,99868) Temp , PP(Ip) , POUT , numbdens , 
     &                    CollFreq , Eref(Molinit), 
     &                    ADJUSTR(MolName(Molinit)) ,
     &                    IDUM0 , IDUM
         ENDIF
        END DO

         WRITE (KOUT,99063) (ADJUSTR(MolName(i)),i=1,NWells+NProds)
         WRITE (KOUT,99060) (HMol(i),i=1,NWells)
         WRITE (KOUT,99061) (i,i,i=1,NWells+NProds)
         WRITE (KFLX,99017) ((j,Jto(j,i),i=1,Nchan(j)),j=1,NWells)
         WRITE (KRAT,99017) ((j,Jto(j,i),i=1,Nchan(j)),j=1,NWells)
c
c     Calculate populations and averages of each molecule

         DO 1100 j = 1 , Ntime                         ! ************ Start of time-loop with Ntime time-steps **********
            IF ( Tspec.EQ.'TIME' .OR. Tspec.EQ.'COLL' ) THEN
                Time = (j-1)*Tstep
              ELSEIF ( Tspec.EQ.'LOGT' ) THEN
                Time = tminlogt*exp( (j-1)*Tstep )
            ENDIF
            
            NumTot = 0.0d+00
            DO 1060 i = 1 , NWells + NProds
               FracMol(i) = 0.0D+00
               AveE(i) = 0.0D+00
               NumTot = NumTot + REAL(Num(i,j),KIND=8)
 1060       CONTINUE
          IF ( NumTot .LT. 0.01 ) NumTot = 1.d+00
 
            DO 1080 Mol = 1 , NWells + NProds
               FracMol(Mol) = REAL(Num(Mol,j),KIND=8)/NumTot             ! Average population in each well
               Var(Mol) = SQRT(FracMol(Mol)*(1.0-FracMol(Mol))/NumTot)
               IF ( Mol.LE.NWells .AND. Num(Mol,j).GT.0 ) THEN
                 AveE(Mol) = EMol(Mol,j)/REAL( Num(Mol,j), KIND=8)               ! Average E in each well
                 DO 1065 k = 1 , Nchan(Mol)
                   kuni(Mol,k,j) = kuni(Mol,k,j)/REAL(Num(Mol,j),KIND=8)       ! average flux coefficients for a given path
                   Flux(Mol,k,j) = kuni(Mol,k,j)*FracMol(Mol)     ! reaction flux via a given path
 1065            CONTINUE
               ENDIF
 1080       CONTINUE
 
            WRITE (KOUT,99062) Time , Time*CollFreq , 
     &                         (FracMol(Mol),Var(Mol),AveE(Mol),Mol=1,
     &                         NWells+NProds)
            WRITE (KRAT,99019) Time , Time*CollFreq ,
     &                         ((kuni(Mol,k,j),k=1,Nchan(Mol)),Mol=1,
     &                         NWells)  ! average rate constants
            WRITE (KFLX,99019) Time , Time*CollFreq ,
     &                         ((Flux(Mol,k,j),k=1,Nchan(Mol)),Mol=1,
     &                         NWells)  ! average rate constants
 1100    CONTINUE ! End of time-step loop
 
         WRITE (KSMM,99062) numbdens , Tlim*CollFreq , 
     &                      (FracMol(Mol),Var(Mol),AveE(Mol),Mol=1,
     &                      NWells+NProds)
      

c ---------------------------------------------------------------------
c     Open and write out energy distributions
c
         OPEN (KDIS,FILE=NameDIST,STATUS='UNKNOWN')     ! vib distributions vs. time
         WRITE (KDIS,99001) AVERSION , ADATE , AVERSION , ADATE
         CALL DateTime(KDIS)
         WRITE (KDIS,99003) TITLE
         WRITE (KDIS,99044) Egrain1 , Emax1 , imax1 , Egrain2 , Emax2 , 
     &                   (Isize-imax1)

c    Write out energy distributions
         WRITE (KDIS,99070) Temp , Tvib , P(Ip) , (Mol, Mol=1,NWells)    ! Write out energy distributions
         Tstep = Tlim/(Mtime-1)
          DO j = 1 , Mtime
           kkmax = 0                           !------------- Find highest energy bin (kkmax) with non-zero population
           k = Ndist
           DO Mol = 1 , NWells
             IF (Mum(Mol,j) .EQ. 0 ) Mum(Mol,j) = 1               ! to avoid divide-by-zero
             DO WHILE ( k.GT.1 .AND. EDIST(Mol,k,j).LE.0.0D+00 )
               k = k - 1
             ENDDO  ! k
             kkmax = MAX( kkmax, k )     ! highest energy bin (kkmax) with non-zero population
           END DO ! Mol                        !------------- 

           DO k = 1 , kkmax
             WRITE(KDIS,99071) (j-1)*Tstep, (j-1)*Tstep*Collfreq, 
     &         (k-1)*DistStep , 
     &         (EDIST(Mol,k,j)/REAL(Mum(Mol,j),KIND=8) ,  Mol=1,NWells)       ! Energy of Edist bin is at the TOP of the bin; first bin top = 0.0
           END DO  ! k

          END DO  ! j

         CLOSE (KDIS)                          ! vib distributions vs. time
c ---------------------------------------------------------------------

! this do loop ensures that the correct Eor and itun information are used when
! printing out the initial energy distribution

        do i=1,Nchan(Molinit)
          IF ( Jto(Molinit,i).EQ.IR ) ichan = i
        end do

         IF ( Icalc.GT.1 ) THEN                              ! If not a delta function
            WRITE (KARY,99068)                               ! chemical activation or thermal distribution
            NCUM = 0
            II = 0
            DO WHILE ( (II .LE. Jsize) .AND. (NCUM .LT. NINT(Rtrials) ))
              II = II + 1
c              if (Icalc.eq.2 .OR. Icalc.eq.4) then                ! THERMAL or EXTERNAL
                 E = Edels + Einit + (II-1)*Hstart           ! start at Einit for THERMAL and EXTERNAL
c              else if (Icalc.eq.3)then                                  ! CHEMACT
c                if(itun(Molinit,ichan).eq.0)then
c                  E = Eor(Molinit,ichan) + (II-1)*Hstart     ! start at reaction threshold for CHEMACT w/o tunneling 
c                else
c                  Einit = Eo(Molinit,ichan) - Emax1
c                  IF ( Einit .LT. 0.0 ) Einit = 0.0
c                  E = Einit + (II-1)*Hstart                  ! start at Egrain1*imax1 below the reaction threshold for 
c                end if  ! itun                                     ! for CHEMACT with tunneling
c              else
c                continue
c              end if ! (Icalc.eq.3)
              NCUM = Nstart(II) + NCUM                                ! cumulative
              WRITE (KARY,99069) E , Pstart(II) , 
     &            DBLE(NCUM)/Rtrials
            END DO  !  DO WHILE loop
         ENDIF !  Icalc

 1200 CONTINUE                          ! End pressure step loop
 
! ----------------------------------------------------------------------
!     DEALLOCATION OF ARRAYS

      jstat = 0
      DEALLOCATE( Dens, STAT=istat)
      jstat = jstat + istat
      DEALLOCATE( COLLUP,  STAT=istat)
      jstat = jstat + istat
      DEALLOCATE( Cnorm,  STAT=istat)
      jstat = jstat + istat
      DEALLOCATE( Rate,  STAT=istat)
      jstat = jstat + istat
      DEALLOCATE( TRate , STAT=istat ) 
      jstat = jstat + istat
      DEALLOCATE( Pstart , STAT=istat ) 
      jstat = jstat + istat
       IF ( jstat .NE. 0 ) THEN
         write(*,*) '*** FATAL: arrays Dens, etc. not DEallocated ***'
         write(*,*) '         STAT=',jstat
         STOP
      ENDIF
! ----------------------------------------------------------------------
!     CLOSE OPEN FILES

      CLOSE (KIN)       ! Input data file
      CLOSE (KOUT)      ! General time-dependent output
      CLOSE (KSMM)      ! Short summary output for fall-off calculations
      CLOSE (KARY)      ! Summary of input, output, and intermediate array variables
      CLOSE (KRAT)      ! average rate 'constants'
      CLOSE (KFLX)      ! relative reactive flux via a given reaction path

c ------------------------------------------------------------------------------------------
c     FORMAT STATEMENTS
c
99001 FORMAT (
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'//8x,
     &'                                John R. Barker'/8x,
     &'      MultiWell-',A7,'         University of Michigan'/8x,
     &'                                Ann Arbor, MI 48109-2143'/8x,
     &'        ',A8,'                jrbarker@umich.edu'/8x,
     &'                                (734) 763 6239'//8x,
     &'      http://clasp-research.engin.umich.edu/multiwell/'//,
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'//8x,
     &'Suggested Literature Citations:'//4x,
     &'a) J.R. Barker, T.L. Nguyen, J.F. Stanton, C. Aieta, M. Ceotto,',
     &/7x,'F. Gabas, T.J.D. Kumar, C.G.L.Li, L.L.Lohr, A. Maranzana,',
     &/7x,'N.F. Ortiz, J.M. Preses, J.M. Simmie, J.A. Sonk, and ',
     &'P.J. Stimac; ',
     &/7x,'MultiWell-',A7,' Software Suite; University of Michigan, ',
     &/7x,'Ann Arbor, Michigan, USA, ',A8,'; ',
     &/7x,'http://clasp-research.engin.umich.edu/multiwell/.',
     &//4x,'b) J.R. Barker, Int. J. Chem. Kinetics, 33, 232-45 (2001)',
     &//4x,'c) J.R. Barker, Int. J. Chem. Kinetics, 41, 748-763 (2009)',
     &//,
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'/)
99002 FORMAT (/'Isize greater than maximum allowed: ',I5,
     &        ' elements.'/'****** Hit RETURN to terminate run *******')
99003 FORMAT (A100)
99017 FORMAT ('  Time   Colls(Mol=1)',150(4X,I2.2,'...',I2.2))
99019 FORMAT (2(1x,1pe10.3),150(1x,1pe10.3))
99044 FORMAT (/'PARAMETERS FOR DOUBLE-ARRAYS:'/
     &    '        Grain (cm-1)      Maximum Energy    Number of grains'
     &    /2(2F20.3,i20/))
99059 FORMAT (/'RESULTS FOR',1pe7.0,' TRIALS')
99060 FORMAT ('        deltaH(cm-1): ',75(5x,0pf10.1,18x))
99061 FORMAT ('   Time  Collisns(Mol=1)',
     &        75(: 3x,'Fract(',I2.2,') +/-Error  <Evib(',I2.2,')>'))
99062 FORMAT (2(1x,1pe10.3),75(3x, 1pe11.3 , 1pe9.1 , 0pf8.0,2x))
99063 FORMAT ('                Name:    ',75(5x,A10,18x))
99064 FORMAT (1x, 1pe10.3 ,'   Trials   = ',I3,' %')
99164 FORMAT (/,1x,'Pressure = ',1pe10.3)
99065 FORMAT ('****Monoenergetic****'/'    Ttrans  = ',f6.1,
     &        ' K',/,'    Press   = ',1pe10.3,1x,A4,/,'          N = ',
     &        1pe10.3,' /cm3',/,
     &        'Coll. Freq. = ',1pe10.3,' /s at E =',0pf8.1,' cm-1',/,
     &        'Initial Molecule: ',A10,/,
     &       'E(initial) =',
     &        0pf10.1,' cm-1',/,'Starting random number seed = ',
     &        I10,/,'Ending   random number seed = ',I10,/)
99066 FORMAT ('****Shifted Thermal (at Tvib)****'/,
     &        '     Tvib   = ',f6.1, ' K',/,
     &        '     Ttrans = ',f6.1,' K',/,
     &        '     Press  = ', 1pe10.3, 1x,A4,/,
     &        '          N = ',1pe10.3,' /cm3',/,
     &        'Coll. Freq. = ',1pe10.3,' /s at E =', 0pf8.1,' cm-1',/,
     &        'Initial Molecule: ',A10,/,
     &        'E(shift) =', 0pf10.1,' cm-1'/,
     &        'Starting random number seed = ', I10,/,
     &        'Ending   random number seed = ',I10,/)
99067 FORMAT ('****Chemical Activation (at Tvib)****'/
     &        '     Tvib   = ',f6.1,' K',/,
     &        '     Ttrans = ',f6.1,' K',/,
     &        '     Press  = ', 1pe10.3,1x,A4,/,
     &        '          N = ', 1pe10.3,' /cm3',/,
     &        'Coll. Freq. = ',1pe10.3,' /s at E =',0pf8.1,' cm-1',/,
     &        'Excited Molecule: ', A10,/,
     &        '   Produced from: ', A10,/,
     &        '   E(shift) =', 0pf10.1,' cm-1'/,
     &        'Starting random number seed = ', I10,/,
     &        'Ending   random number seed = ',I10,/)
99068 FORMAT (/'CUMULATIVE DISTRIBUTION FXNS FOR STARTING ENERGIES',
     &        ' (Pstart)',//,
     &        '    ENERGY   P(Theory)  P(selected) ')
99069 FORMAT (1x,f10.1 , 2(2x,f9.7) )
99070 FORMAT (/'Energy Distribution (non-zero elements)',/,
     & '  (E = Active Energy at TOP of energy bin, relative to ZPE)',
     & //,' Ttrans=',f6.1,1x,'K',5x,'Tvib=',f6.1,1x,'K',5x,
     & 'P= ',1PE10.3,1x,'bar',//,
     & ' Time       Collisions    E(cm-1)    ',75(: 'Well=',I2.2,4x))
99071 FORMAT (1x,1pe10.3,2(1x,0pf10.1),75(1x,1pe10.3) )
99173 FORMAT ( A5 )
99074 FORMAT (/,'INPUT DATA FILE CHECK: ',I2,' REACTION CHANNELS'/
     &'#  Mol ito TS          RR       j       Qel   l   AA         ',
     &'EE       KEYWORDS')
99868 FORMAT ('**** ENERGY DISTRIBUTION FROM EXTERNAL FILE ****'/,
     &        '     Ttrans = ',f6.1,' K',/,'     Press  = ',
     &        1pe10.3,1x,A4,/'          N = ',
     &        1pe10.3,' /cm3',/,
     &        'Coll. Freq. = ',1pe10.3,' /s at E =',0pf10.1,' cm-1',/,
     &        'Initial Molecule: ',
     &        A10,/,'Starting random number seed = ',
     &        I10,/,'Ending   random number seed = ',I10,/)

      END PROGRAM multiwell
