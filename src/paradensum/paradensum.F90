
      PROGRAM paradensum
!    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!    
!    LICENSE NOTICE
!
!    Copyright (C) 2019 Michele Ceotto, Chiara Aieta, Fabio Gabas, 
!                       Thanh Lam Nguyen, and John R. Barker
!
!    Contact:
!
!    Michele Ceotto  (email: michele.ceotto@unimi.it)
!    Dipartimento di Chimica
!    Università degli Studi di Milano
!    via Golgi 19, Milano 20133, ITALY
!
!    or:
!
!    John R. Barker   (email: jrbarker@umich.edu)
!    Department of Atmospheric, Oceanic, and Space Sciences
!    College of Engineering
!    University of Michigan
!    Ann Arbor, MI 48109
!
!    http://clasp-research.engin.umich.edu/multiwell
!    This program is free software; you can redistribute it and/or
!    modify it under the terms of the GNU General Public License (version 2)
!    as published by the Free Software Foundation.
!   
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!   
!    See the 'ReadMe' file for a copy of the GNU General Public License,
!    or contact:
!   
!    Free Software Foundation, Inc.
!    59 Temple Place - Suite 330
!    Boston, MA 02111-1307, USA.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                 
!                         PROGRAM paradensum    
!
!                               by      
!
!           by Michele Ceotto, Chiara Aieta, Fabio Gabas,
!               Thanh Lam Nguyen, and John R. Barker    
!                                                                  
!           ***PARallel Anharmonic DENsities and SUM of states***
!
!                             based on
!
!          The parallel implementation* of the Wang-Landau algorithms 
!          for densities of states**
!
!    Literature Citations:                                                    
!    *parallel implementation of anaharmonic density of states
!    C. Aieta, F. Gabas and M. Ceotto, J. Phys. Chem. A., 120(27), 
!    4853-4862 (2016).
!                                                                  
!    **Density of states algorithms
!    F. Wang and D. P. Landau, Phys. Rev. Letters 86, 2050-2053 (2001).           
!    M. Basire, P. Parneix, and F. Calvo, J. Chem.Phys. 129, 081101 (2008).    
!    T. L. Nguyen and J.R. Barker, J. Phys. Chem. A., 114, 3718–3730 (2010).
!    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
      USE decl_alloc
      IMPLICIT NONE
      include "mpif.h"

      ! MPI data
      !INTEGER :: my_rank, num_procs, err	!NEW
      integer, allocatable :: my_wind(:)
      INTEGER :: jj
      INTEGER :: H_size, g_size  
      REAL(8) :: time1, time2, time_init_start, time_init_end !NEW
      CHARACTER(len=4) :: str2
      CHARACTER(len=10) :: str
      CHARACTER(len=20) :: filename
      !INTEGER(4) ngrains_corr  !NEW
      INTEGER(4) seed_modifier
      INTEGER(4) ntop !NEW
      REAL(8) B1, B2, X1, X2 !NEW

!     Variabili aggiuntive per la nuova implementazione      
      INTEGER(4) chkp, count_init!, estart_check, ckenrg
      INTEGER(4) stepstoflat, counts_stepstoflat, flatcheck!, chiara, npermedia !NEW
      REAL(8) H_flatness
      INTEGER(4) writing 
      REAL(8) flatness
      !REAL(8) perc_wind_overlap, flatness !NEW
      INTEGER(4) windcount, walkercount, fcount
      INTEGER(4) iy, iv(32) 
      REAL(8) randn
!     join g from each window
      !REAL(8), DIMENSION(:), ALLOCATABLE :: head_tail_T_comp
      REAL(8) :: dTa, dTb!, wind_norm_s !NEW 
      !INTEGER(4) :: Ttot_bookmark, old_cut_point, cut_point_array(1),cut_point
      INTEGER(4) :: Ttot_bookmark, old_cut_point,cut_point
      !INTEGER(4) :: T_shift,ind, ind_old, ind_ref, ind_ref_old 
      INTEGER(4) :: T_shift 
      INTEGER(4), DIMENSION(:), ALLOCATABLE :: record_Ttot_bookmark
      !REAL (8) :: Emin_global, Emax_global, !NEW 
      !REAL(8) :: scal, sh
      !REAL(8) :: ngrains_per_chunk_small, ngrains_per_chunk_big !NEW
      REAL(8) :: deriv, deriv_old


      CHARACTER(len=7) AVERSION 
      CHARACTER(len=8) ADATE 
      PARAMETER ( AVERSION='2022', ADATE='Jan 2022' )
 
      CHARACTER(46) cut , dumdum
      PARAMETER (cut='**************INPUT DATA SUMMARY**************')
      
      CHARACTER(10) fname , fnameT
      
      CHARACTER(9) chekpoint
      CHARACTER chekname*15     
 
      INTEGER istat , itest , iunit

      INTEGER(4) ns, i, j, k!, MaxTrials
      !INTEGER(4) l , nstart , Nacc, Nrej, Nout
      INTEGER(4) l, Nacc, Nrej, Nout
      INTEGER(4) iters, nold, nnew
      !INTEGER(4) ndmax, nh, ig, nnt, nvmax, KIN
      INTEGER(4) ndmax, ig, nvmax, KIN
      PARAMETER(KIN=2)
      
      !REAL(8) Eu, trialmax, ttmax         !NEW
      !REAL(8) Emin, Emax, Estart, energy  !NEW
      !REAL(8) Estart, energy, energy_old   !NEW
      REAL(8) energy, energy_old   !NEW
      REAL(8) enrg, zpeMorse, zpeTotal
      REAL(8) boxzpe !NEW
      REAL(8) Etop
      REAL(8) fx, Ered, Q0, Q1, Q2, Cx, Hx, Sx

      REAL(8) p , p0 , f , fo , test
      !REAL(8) acc , hmax , hmin , have, hav2, y , var
      REAL(8) acc, y ! , var
 
      CHARACTER(100) TITLE1 , TITLE2
      CHARACTER(2) WW
      !CHARACTER(6) KEYWORD, KEYWORD2, Eunits !NEW togli KEYWORD
      CHARACTER(5) windbal_key !NEW

      CHARACTER(3) IDOF(100)
      CHARACTER(5) Vhr(100) , Bhr(100)
      !CHARACTER(4) VROT !NEW
      CHARACTER(4) VROTIN, VROTOUT !NEW
      
      REAL(8) DUM1, DUM2
      REAL(8) CV(100,100), CB(100,100), CVt(100), CBt(100) 
      INTEGER(4) Nsep
      
      REAL(8) WE , XE
      
      !REAL(8) BB , BG , FAC , RG , GAM
      !INTEGER N , NR

      INTEGER(4) NY, NZ 
      INTEGER(4) nvmaxxyz, JMAX 
      REAL(8) energyxyz
     
      REAL(8) B
      INTEGER(4) lenstr

      INTEGER(4) NDIM, NSYMM, JK
 
      !REAL(8) Egrain1 , Egrain2 , Emax1 , Emax2 !NEW 
      INTEGER(4) imax1 , imax2 , Isize , NMAX
      PARAMETER (NMAX=50001)  

      REAL(8) Egrain1T , Emax2T 
      INTEGER(4) imax1T , IsizeT

      INTEGER(4) ntest, Ntx
      
      !REAL(8) TT , su , beta0 , RR
      REAL(8) TT, RR
      PARAMETER (RR=0.695038916D0)
      !NEW
      EXTERNAL energy, energy_old, energyxyz, ndmax, nvmax, nvmaxxyz
      
      SAVE

!     INITIALIZE MPI PROCESSES

      call MPI_INIT( err )
      call MPI_COMM_RANK( MPI_COMM_WORLD, my_rank, err )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, num_procs, err )

!     Number of Windows set equal to numbero of MPI porcs
      nwind = num_procs

      str2='Rank'
      WRITE (str, '(I10)') my_rank
      str = ADJUSTL(str)
      filename = str2 // TRIM(str) //'.txt'

      OPEN(10+my_rank,FILE=filename,STATUS='UNKNOWN')

      OPEN(61,FILE='Windows_info.txt',STATUS='UNKNOWN')

      OPEN (2,FILE='paradensum.dat',STATUS='OLD')
      READ (2,*) fname

      OPEN (3,STATUS='UNKNOWN',FILE=fname(1:lenstr(fname))//'.dens')  ! Succinct output for MultiWell input
      OPEN (4,FILE='paradensum.out',STATUS='UNKNOWN')
      OPEN (5,FILE=fname(1:lenstr(fname))//'.qvib',STATUS='UNKNOWN')
      
      READ (2,9002)  TITLE1
      READ (2,9002)  TITLE2

      WRITE(3,99030) cut                                  ! Start of data summary block in ____.dens file
      WRITE(5,99030) cut                                  ! Start of data summary block in ____.qvib file

      DO iunit = 3 , 5
        WRITE(iunit,9911) AVERSION , ADATE
        WRITE(iunit,9912) AVERSION , ADATE 
        WRITE(iunit,99023) fname
        WRITE(iunit,9002) TITLE1 
        WRITE(iunit,9002) TITLE2 
        CALL DateTime(iunit)                               ! print output to unit=iunit
      END DO
      
      WRITE (4,9005)

      READ(2,*)  ns , NY, NZ, WW
 
      CALL alloc_3(ns)

      CALL READ_WXYZ(KIN , ns , NY, NZ, WW)                ! read frequencies and anharmonicities
      
      READ(2,*) NSEP , VROTIN                                ! No. of separable DOF and keyword for rotation units
      CALL ucase ( VROTIN )   ! convert to upper case !NEW        
      VROTOUT = 'AMUA'        !NEW perché?
    
      CALL alloc_4(NSEP,ns)
 
      ntop = 0 !NEW 
      DO I=ns+1, ns+NSEP
         READ(2,*) MODE(I), IDOF(I), DUM1, DUM2, NG(I)
         CALL ucase ( IDOF(I) )   ! convert to upper case NEW
         IF (IDOF(I).EQ.'HRA') THEN
            IF ( NG(I) .LT. 0) THEN
                READ(2,*) NSIG(I)
                NG(I) = -NG(I)
              ELSE
                NSIG(I) = NG(I)
            ENDIF
            HRfreq(I)=DUM1
            HRI(I)=DUM2
            HRB(I)=16.85763d0/HRI(I)
            HRVo(I)=((HRfreq(I)/NG(I))**2)/HRB(I)
         ELSEIF (IDOF(I).EQ.'HRB') THEN
            IF ( NG(I) .LT. 0) THEN
                READ(2,*) NSIG(I)
                NG(I) = -NG(I)
              ELSE
                NSIG(I) = NG(I)
            ENDIF
            HRfreq(I)=DUM1
            HRVo(I)=DUM2
            HRB(I)=((HRfreq(I)/NG(I))**2)/HRVo(I)
            HRI(I)=16.85763d0/HRB(I)
         ELSEIF (IDOF(I).EQ.'HRC') THEN
            IF ( NG(I) .LT. 0) THEN
                READ(2,*) NSIG(I)
                NG(I) = -NG(I)
              ELSE
                NSIG(I) = NG(I)
            ENDIF
            HRI(I)=DUM1
            HRVo(I)=DUM2
            HRB(I)=16.85763d0/HRI(I)
            HRfreq(I)=NG(I)*sqrt(HRVo(I)*HRB(I))
         ELSEIF (IDOF(I).EQ.'HRD') THEN
            NVV(I)=INT(DUM1)
            NBB(I)=INT(DUM2)
            READ(2,*) Vhr(I), NGV(I), Phav(I), (CV(I,J), J=1, NVV(I))
            READ(2,*) Bhr(I), NGB(I), Phab(I), (CB(I,J), J=1, NBB(I))
            CALL ucase ( Vhr(I) )   ! convert to upper case NEW
            CALL ucase ( Bhr(I) )   ! convert to upper case NEW    
         ELSEIF (IDOF(I).EQ.'ROT') THEN
            RI(I) = DUM1
            NRSN(I) = DUM2
            IDIM(I) = NG(I)
            CALL rotunits( VROTIN , RI(I), VROTOUT )
         ELSEIF (IDOF(I).EQ.'QRO') THEN
            RI(I) = DUM1
            NRSN(I) = DUM2
            IDIM(I) = NG(I)
            CALL rotunits( VROTIN , RI(I), VROTOUT )
         ELSEIF ( IDOF(I).EQ.'TOP' ) THEN !NEW 15lines                              ! Symmetric Top with 2D and 1D rot constants B2 and B1
            ntop = ntop + 1
            IF ( ntop .GT. 1) THEN
              write (*,*) 'FATAL: only one symmetric top (top) allowed'
              write (3,*) 'FATAL: only one symmetric top (top) allowed'
              write (4,*) 'FATAL: only one symmetric top (top) allowed'
            ENDIF
            B2 = DUM1
            X2 = B2
            call rotunits( VROTIN , B2 , 'CM-1' )     ! unit conversion
            call rotunits( VROTIN , X2 , 'AMUA' )     ! unit conversion
            B1 = DUM2
            X1 = B1
            call rotunits( VROTIN , B1 , 'CM-1' )     ! unit conversion
            call rotunits( VROTIN , X1 , 'AMUA' )     ! unit conversion
         ELSEIF (IDOF(I).EQ.'KRO') THEN
            RI(I) = DUM1
            NRSN(I) = DUM2
            IDIM(I) = NG(I)
            CALL rotunits( VROTIN , RI(I), VROTOUT )
         ELSEIF (IDOF(I).EQ.'BOX') THEN
            RI(I) = DUM1
            IDIM(I) = NG(I)
         ELSEIF (IDOF(I).EQ.'VIB') THEN
            IF ( NG(I) .EQ. 1 ) THEN
              RI(I) = DUM1
              RJ(I) = DUM2
            ELSE
              WRITE(*,*) '********************************************'
              WRITE(*,*) 'FATAL: vibs must not be degenerate'
              WRITE(*,*) '       d.o.f. #', I
              WRITE(*,*) '********************************************'
              WRITE(4,*) '********************************************'
              WRITE(4,*) 'FATAL: vibs must not be degenerate'
              WRITE(4,*) '       d.o.f. #', I
              WRITE(4,*) '********************************************'
              STOP
            ENDIF
         ELSE
            WRITE(*,*) '*********************************************'
            WRITE(*,*) 'FATAL: degree of freedom type not recognized'
            WRITE(*,*) '       d.o.f. #', I
            WRITE(*,*) '*********************************************'
            WRITE(4,*) '*********************************************'
            WRITE(4,*) 'FATAL: degree of freedom type not recognized'
            WRITE(4,*) '       d.o.f. #', I
            WRITE(4,*) '*********************************************'
            STOP
        ENDIF
      END DO
 
      READ (2,*) Egrain1 , imax1 , Isize , Emax2 !, KEYWORD
      IF (Emax2/Egrain1 .GT. NMAX) THEN
        WRITE(*,*) 'FATAL: **Grain Size Too Small****'
        WRITE(4,*) 'FATAL: **Grain Size Too Small****'
        STOP
      ENDIF
      IF ( imax1 .GT. Isize-2 ) THEN
       WRITE (*,*) '*** FATAL:  Imax1 must be ≤ (Isize-2) ***'
       WRITE (4,*) '*** FATAL:  Imax1 must be ≤ (Isize-2) ***'
       STOP
      ENDIF      
      
      READ(2,*) chekpoint, chekname
      CALL ucase ( chekpoint )   ! convert to upper case

!     New lines in paradensum.dat
      READ (2,*) !Blank line after checkpoint line
      READ (2,*) nwalkers
      READ (2,*) perc_wind_overlap
      READ (2,*) flatness
      READ (2,*) writing
      READ (2,*) seed_modifier
      READ (2,*) windbal_key
      
      CLOSE (UNIT=2)
!
!   ******* END INPUT *******
!

!---- Double Array setting from input data -------------------------------------------------------------
         
      Emax = Emax2 + Egrain1
      Emin = 0.0                               ! Bottom of Wang-Landou energy window
      Emax1 = Egrain1*(imax1-1)
      Egrain2 = Emax2/(Isize-imax1-1)
      imax2 = Isize - imax1
 
      JMAX = 1 + INT(Emax2/Egrain1)            ! JMAX compare solo qui
      JMAX = MIN(JMAX,NMAX)

!--------------------------------------------------------------------------------------------------------

      CALL convib( ns , WW )                                    ! convert input frequencies to other forms

      WRITE(4,9003) ns
      WRITE(3,9003) ns
      WRITE(5,9003) ns

      DO i = 1 , ns
        D(i) = 0.0
        IF (xa(i,i) .LT. 0.0) D(i) = -0.25D+00*w0(i)**2/xa(i,i) ! Bond dissociation energy
        WRITE(4,9004) i, wa(i), w0(i), wf(i), D(i)
        WRITE(3,9004) i, wa(i), w0(i), wf(i), D(i)
        WRITE(5,9004) i, wa(i), w0(i), wf(i), D(i)
      END DO

      WRITE(4, 9013) ave, av0, avf, zpe                         ! zero point energy and average frequencies of coupled DOF

      WRITE(4,9006)
      WRITE(5,9006)
      WRITE(3,9006)
      WRITE(4,9043) (i, i=1,ns)
      WRITE(3,9043) (i, i=1,ns)
      DO j = 1 , ns
        WRITE(4,9042) j, (xa(j,i), i=1,ns)
        WRITE(3,9042) j, (xa(j,i), i=1,ns)
        WRITE(5,9042) j, (xa(j,i), i=1,ns)
      END DO

!    Write out separable degrees of freedom

!      WRITE(4,*)
!      WRITE(3,*)
      WRITE(4,99050)
      WRITE(3,99050)
      WRITE(5,99050)
      DO I=ns+1, ns+Nsep
            IF (IDOF(I).EQ.'HRA' .OR. IDOF(I).EQ.'HRB' .OR. & 
                IDOF(I).EQ.'HRC') THEN
              DO iunit = 3, 5
              WRITE(iunit,9030) MODE(I), HRfreq(I), HRI(I), NG(I), & 
                                NSIG(I), HRB(I), HRVo(I)
             END DO
            ELSE IF (IDOF(I).EQ.'HRD') THEN
             DO iunit = 3, 5
                WRITE(iunit,9031) MODE(I), NG(I)
                WRITE(iunit,9032) Vhr(I), NGV(I), Phav(I), &
                                 (CV(I,J), J=1, NVV(I))
             END DO
             IF(Bhr(I).EQ.'BHRD1') THEN
               DO iunit = 3, 5
                  WRITE(iunit,9033) Bhr(I), NGB(I), Phab(I), &
                                   (CB(I,J), J=1, NBB(I))
               END DO
              ELSE
               DO iunit = 3, 5
                  WRITE(iunit,9033) Bhr(I), NGB(I), Phab(I), &
                                   (CB(I,J), J=1, NBB(I))
               END DO  
              END IF
            ELSEIF (IDOF(I).EQ.'ROT') THEN
                  B = 16.85763d0 / RI(I)
             DO iunit = 3, 5
                  WRITE(iunit,99044) MODE(I), RI(I), NRSN(I), IDIM(I), B
             END DO
            ELSEIF (IDOF(I).EQ.'QRO') THEN
                  B = 16.85763d0 / RI(I)
             DO iunit = 3, 5
                  WRITE(iunit,99041) MODE(I), RI(I), NRSN(I), IDIM(I), B
             END DO
            ELSEIF ( IDOF(I).EQ.'TOP' ) THEN                               ! Symmetric Top with 2D and 1D rot constants B2 and B1
               DO iunit = 3, 5
                 WRITE (iunit,99905) MODE(I) , X2 , B2 , X1 , B1 , NG(I)
               END DO
            ELSEIF (IDOF(I).EQ.'KRO') THEN
                  B = 16.85763d0 / RI(I)
             DO iunit = 3, 5
                  WRITE(iunit,99042) MODE(I), RI(I), NRSN(I), IDIM(I), B
             END DO
            ELSEIF (IDOF(I).EQ.'BOX') THEN
             DO iunit = 3, 5
                  WRITE(iunit,99043) MODE(I), RI(I), IDIM(I)
             END DO                  
            ELSEIF (IDOF(I).EQ.'VIB') THEN
             DO iunit = 3, 5
                  WRITE(iunit,99045) MODE(I), RI(I), RJ(I)
             END DO                  
            ENDIF

      END DO

      DO iunit = 3, 5
!        NEW: Write out parallel execution parameters
         WRITE (iunit,*) 
         WRITE (iunit,*) '------------------------------'
         WRITE (iunit,*) 'Parallel execution parameters:'
         WRITE (iunit,*) '------------------------------'
         WRITE (iunit,99052) nwind
         WRITE (iunit, 99053) nwalkers
         WRITE (iunit,99054) perc_wind_overlap
         WRITE (iunit,99055) flatness
         IF (writing == 0) THEN
            WRITE (iunit,*) 'Additional writinig disabled.'
         ELSE
            WRITE (iunit,*) 'Additional writinig enabled.'
         END IF
         WRITE (iunit,99056) seed_modifier
         WRITE(iunit,99051)
      END DO

!-----------------------------------------------------------------------------------------
!
!         WANG-LANDAU PARAMETERS
!
!-----------------------------------------------------------------------------------------
      
      p0 = 1.0/ns                                       ! Scaled according to number of vibrations
      IF ( p0 .GT. 0.25d+00 ) p0 = 0.25d+00
      p = p0                                            
      iters = 21                                        ! Always
      Etop = Emax + Egrain1                             
      fo = 1.0d+00
      ngrains = (Emax-Emin)/Egrain1 + 1                 ! should be same as JAMX (when Emin = 0.0)
      if(my_rank.eq.0) WRITE(61,*)'Ngrains: ', ngrains
      nold = 1 !NEW

!----- ENERGY WINDOWS SET UP -----------------------------------------------------------------------------
!     Each window must have an integer number of bins

      CALL alloc(nwind)

      CALL wind_setup(windbal_key) !NEW

!     Now gtot can be allocated
!     since the total number of grains has been computed
      call alloc_2(ngrains, ngrains_to_add_max, &
               ngrains_per_wind_max, nwalkers, nwind,ns) !NEW

!     These variable will be used for the MPI_Reduce
      !H_size = ngrains_per_wind(nwind) * nwalkers * nwind
      !g_size = ngrains_per_wind(nwind) * nwalkers * nwind
      H_size = ngrains_per_wind_max * nwalkers * nwind !NEW
      g_size = ngrains_per_wind_max * nwalkers * nwind !NEW

!     Till now there are: 
!     - an array, lowbound, with the lower limit of each window
!     - an array, upbound, with the upper limit of each window
!     - a table of histograms, H, one for each walker in each window 
!     - a table of densisty of statetes, g, one for each walker in each window 
 
!     inizialize checkpoint control variable
      chkp = 1
      IF( chekpoint .EQ. 'CHEKSTART') chkp = 0
     
      IF (chkp .eq. 1)  THEN   ! do NOT start from saved chekpoint file 

!     PARALLEL SECTION
!     Windows are distributed to MPI processes
      ALLOCATE(my_wind(nwind))
      my_wind = 0 

      DO jj=1,nwind 
         my_wind(jj) = MOD(jj,num_procs)
      END DO

      DO windcount = 1, nwind !<--------------------------------------------------------------------------------Start cycle on WINDWS 
         IF(my_wind(windcount) .eq. my_rank ) THEN      !<--------------------------------------------START SINGLE MPI THREAD SECTION

         time1 = MPI_Wtime()         

         if( writing .eq. 1) then
         WRITE(10+my_rank,*) 'I am rank number', my_rank, 'and I am computing &
                        window number ', windcount
         end if

!        Find out window boundaries
         Emin = lowbound(windcount)
         Emax = upbound(windcount)

         DO fcount = 1, 21 !<--------------------------------------------------------------------Start cycle on modification factor f

!           Inizialize Histogram
            H(:,:,windcount) = 0.d0

            DO walkercount = 1, nwalkers !<-----------------------------------------------------------------------Start walkers Cycle

!	       inizialize random number generator independently for each walker      
               idum = -1147483646 + seed_modifier + walkercount + (windcount-1)*nwalkers
               call ran1(idum, iv, iy, y)

!              Initialize quantum numbers if it is the first f-iteration !NEW 42 lines
!              Otherwise retrieve information
               IF ( fcount == 1 ) THEN

!                   init
                    f = fo !fo=1
                    count_init = 0
                    enrg = 0
                    nv = 0

                  time_init_start = MPI_Wtime()
                  DO
                     count_init = count_init + 1

!                    chose randm a vib configuration 
                     CALL ran1(idum, iv, iy, randn)
                     nv( CEILING(ns*randn) ) = nv( CEILING(ns*randn) ) + 1

!                    Energy calculation
                     IF ((NY.EQ.0).AND.(NZ.EQ.0)) THEN
                        enrg = energy(ns)
                     ELSE
                        enrg = energyxyz(ns, wa, xa, ya, za, nv)
                     ENDIF

!                    Dissociation check (ntest=0 pass, ntest=1, fail)
                     IF ((NY.EQ.0).AND.(NZ.EQ.0)) THEN
                        CALL ckderiv( ns, ntest)
                     ELSE
                        CALL ckderivxyz( ns, ntest )
                     END IF

                     IF (ntest == 1) THEN
                        nv( CEILING(ns*randn) ) = nv( CEILING(ns*randn)) - 1
                     ELSE
                        IF (enrg.ge.Emin.and.enrg.le.Emax) THEN
                           EXIT
                        ELSE IF (enrg .gt. Emax) THEN
                           nv( CEILING(ns*randn) ) = nv(CEILING(ns*randn) ) - 1
                        END IF
                     END IF

                    ! WRITE(10+my_rank,*) 'No. random calls=',count_init
                     IF (count_init == 100000000) STOP &
                        "Starting point not found"
                  END DO

                  time_init_end = MPI_Wtime()
                  WRITE(10+my_rank,*) 'Starting point selection time:' &
                          , time_init_end - time_init_start

                  nvold = nv

!-------------------------------------------------------------------
!                 After the starting point selection
!                 find the energy and the bin corresponding to
!                 the selected quantum number combination
                  
                  IF (enrg == Emax) THEN
                       nold = FLOOR( (enrg -Emin) / Egrain1, 4) 
                  ELSE
                       nold = FLOOR( (enrg -Emin) / Egrain1, 4) + 1
                  END IF

                  if( writing .eq. 1) then
                    WRITE(10+my_rank,*) 'window ', windcount, 'walker ',walkercount
                    WRITE(10+my_rank,*) 'Emin - Emax ', Emin, ' - ', Emax
                    WRITE(10+my_rank,*) 'energy ', enrg
                    WRITE(10+my_rank,*) 'quantum numbers ', nvold(:)
                    WRITE(10+my_rank,*) 'nold ', nold
                  end if
 
               ELSE

                  idum = idum_cache(walkercount)
                  iy = iy_cache(walkercount)
                  iv(:) = iv_cache(:, walkercount)
                  nold = nold_cache(walkercount)
                  nvold(:) = nvold_cache(:, walkercount)

               END IF
 
               if( writing .eq. 1) WRITE(10+my_rank,9020) fcount, iters, f, p 

!	       Start Random Walk
               stepstoflat = 0
               counts_stepstoflat = 0
               flatcheck = 0
               Nrej = 0
               Nacc = 0
               Nout = 0

               DO WHILE ( flatcheck .eq. 0 ) !<--- Random Walk

                  stepstoflat = stepstoflat + 1

!		  Select new quantum numbers
                  DO k = 1 , ns
                     call ran1(idum, iv, iy, randn)
                     test = randn
                     IF ( test .LE. p  ) THEN 
                        nv(k) = nvold(k) - 1                                ! down-step
                        IF ( nv(k) .LT. 0 ) nv(k) = 0
                     ELSE IF ( test .GT. p .AND. test .LE. 2.0d0*p ) THEN   ! up-step
                        nv(k) = nvold(k) + 1                                ! up-step  
                     ELSE
                        nv(k) = nvold(k)                                    ! sit-tight
                     END IF
                  END DO

                  IF ((NY.EQ.0).AND.(NZ.EQ.0)) THEN
                     enrg = energy( ns)
                  ELSE
                     enrg = energyxyz( ns, wa, xa, ya, za, nv )         ! new energy above zpe
                  ENDIF
  
                  IF ( enrg.GE.Emin .AND. enrg.LE.Emax ) THEN           ! check energy boundary
                     ntest = 0
                     IF ((NY.EQ.0).AND.(NZ.EQ.0)) THEN
                        CALL ckderiv( ns, ntest)           ! check derivatives for X-matrix 
                     ELSE
                        CALL ckderivxyz( ns, ntest )! check derivatives for XYZ-matrix 
                     END IF
    
                  ELSE
                    ntest = 1                                           ! fail test, ntest = 1
                  END IF

!                 Accept or reject? 
!                 - pass --> ntest = 0 
!                 - fail --> ntest = 1
                  IF ( ntest.EQ.0 ) THEN
                     IF (enrg == Emax) THEN
                         nnew = FLOOR( (enrg -Emin) / Egrain1, 4) 
                     ELSE
                         nnew = FLOOR( (enrg -Emin) / Egrain1, 4) + 1
                     END IF

                     acc = EXP( g(nold, walkercount, windcount)-g(nnew, walkercount, windcount) ) ! acceptance probability
                     call ran1(idum, iv, iy, randn)
                     test = randn

                     IF ( test .LE. acc ) THEN                  ! accepted for a move
                        Nacc = Nacc + 1     
                        H(nnew, walkercount, windcount) = H(nnew, walkercount, windcount) + 1
                        g(nnew, walkercount, windcount) = g(nnew, walkercount, windcount) + f
                        nold = nnew                             ! energy grain
                        DO l = 1 , ns
                           nvold(l) = nv(l)                     ! quantum numbers
                        END DO
                     ELSE                                       ! rejected: sits tight
                        Nrej = Nrej + 1
                        H(nold, walkercount, windcount) = H(nold, walkercount, windcount) + 1
                        g(nold, walkercount, windcount) = g(nold, walkercount, windcount) + f
                     END IF
                  ELSE                                          ! rejected: tried to step out of energy range
                     Nout = Nout + 1
                     H(nold, walkercount, windcount) = H(nold, walkercount, windcount) + 1
                     g(nold, walkercount, windcount) = g(nold, walkercount, windcount) + f
                  END IF             

!                 Flatness check every 10000 steps
                  IF ( MOD(stepstoflat, 10000) == 0) THEN        
                        IF ( MINVAL( H(:,walkercount, windcount), MASK=H(:,walkercount,windcount) .GT. 0 ) >= & !NEW
                        !( REAL( SUM(H(:, walkercount, windcount), MASK=H(:,walkercount, windcount) .GT. 0 ), 8 ) &
                        ( SUM( REAL(H(:, walkercount, windcount), 8) , MASK=H(:,walkercount, windcount) .GT. 0 ) & !NEW
                        !/ REAL( COUNT(H(:, walkercount, windcount).GT.0 ), 8 )*flatness) ) THEN
                        / COUNT( REAL( H(:, walkercount, windcount) , 8) .GT.0 ) * flatness) ) THEN
                        flatcheck = 1 !NEW
                     ELSE
                        flatcheck = 0
                     END IF
                  END IF

!                 Check if the flatness isn't satisfied after 200000000 steps
                  IF (stepstoflat == 2000000000) THEN
                  !flatcheck = 1 NEW
                  stepstoflat = 0
                  counts_stepstoflat = counts_stepstoflat + 1
                  if( writing .eq. 1) WRITE(10+my_rank,*)  &
                   'Integer limit exceeded! Exit from random walk.'
                  END IF

               END DO !Random walk

               IF ( writing .eq. 1) THEN
                  WRITE(10+my_rank,*)'rank',my_rank,'f',fcount, &
!                            'Steps to flatness: ', stepstoflat
                       'Steps to flatness: ', stepstoflat, '+ 2E+09 * ',counts_stepstoflat
                 WRITE(10+my_rank,*) 'Nout ',Nout, 'Nrej ',Nrej, 'Nacc ', Nacc
                 WRITE(10+my_rank,*)
                 WRITE(10+my_rank,*) 'Actual Flatness', REAL(MINVAL(H(:,walkercount, windcount), MASK=H(:,walkercount,windcount) &
                          .GT. 0 ) ,8) / ( SUM( REAL( H(:, walkercount,windcount),8), MASK=H(:,walkercount, windcount) .GT. 0 ) &
                                                /  COUNT( REAL(H(:,walkercount, windcount),8) .GT.0))
                 WRITE(10+my_rank,*)

               END IF

!              Save walker information
               idum_cache(walkercount) = idum
               iv_cache(:,walkercount) = iv(:)
               iy_cache(walkercount) = iy
               nold_cache(walkercount) = nold
               nvold_cache(:,walkercount) = nvold(:)

            END DO !End walkers cycle

!           check overall flatness
            DO i = 1, nwalkers
                 H_check(:) = H_check(:) + H(:, i, windcount)
            END DO
            H_flatness = MINVAL( H_check, MASK=H_check .GT. 0 ) / &
                        ( REAL( SUM(H_check, MASK=H_check .GT. 0 ), 8 ) &
                        / REAL( COUNT(H_check.GT.0 ), 8 ) )

            WRITE(10+my_rank, *) 'H_flatness', H_flatness

!           Average g among walkers in the same window
            g_ave(:) = 0.d0
            DO i = 1, ngrains_per_wind(windcount)
               DO j = 1, nwalkers
                  g_ave(i) = g_ave(i) + g(i, j, windcount)
               END DO
               g_ave(i) = g_ave(i) / nwalkers
            END DO

!           Ridistribute average g to all walkers in the window
            DO i = 1, nwalkers
               DO j = 1, ngrains_per_wind(windcount)
                  g(j, i, windcount) = g_ave(j)
               END DO
            END DO

            f = f/2.d0

         END DO !End cycle on the modificaiton factor f

!        At this stage all walkers in the same window
!        have the final g of the window

         idum_cache = 0
         iv_cache = 0
         iy_cache = 0
         nold_cache = 0
         nvold_cache = 0

!        Find the first non-zero element in g for every window in the overlap part
!        and compute T(I) = exp( g(I) - g(ref) ) for every window

          ref(windcount) = 1 

          DO WHILE ( g(ref(windcount), 1, windcount) == 0.d0 )

             ref(windcount) = ref(windcount) + 1

          END DO

          if( writing .eq. 1) WRITE(10+my_rank,*)'wind', windcount,'ref', ref(windcount)

          IF (ref(windcount) > ngrains_to_add(windcount) .AND. nwind > 1 ) &
             STOP "*** Try with wider overlap ***"

          DO I = 1, ngrains_per_wind(windcount)

             IF (g(I,1,windcount) /= 0.d0) THEN

                Tw(I, windcount) = exp( g(I,1,windcount) - g(ref(windcount),1,windcount) )

             ELSE

                Tw(I, windcount) = 0.d0

             END IF

          END DO

         time2 = MPI_Wtime()
         wind_time(windcount) = time2 - time1
         if(writing .eq.1) WRITE(10+my_rank, *) 'wind', windcount,'time', time2 - time1

         CLOSE(10+my_rank)

         END IF   !<-------------------------------------------------------------------------------- START SINGLE MPI THREAD SECTION
      END DO 

      CALL DateTime(4) ! print output to unit=4

      call MPI_Reduce(H,H_red,H_size,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,err) 
      call MPI_Reduce(g,g_red,g_size,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,err)

      call MPI_Reduce(ref,ref_red,nwind,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,err) 
      call MPI_Reduce(Tw,Tw_red,ngrains_per_wind(nwind)*nwind,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,err)


!     Retrieve absolute minumum energy
      Emin = lowbound(1)
     
      DEALLOCATE (lowbound, STAT = DeAllocateStatus)
      IF (DeAllocateStatus /= 0) STOP "*** lowbound not deallocated ***"
      DEALLOCATE (upbound, STAT = DeAllocateStatus)
      IF (DeAllocateStatus /= 0) STOP "*** upbound not deallocated ***"
      DEALLOCATE (g_ave, STAT = DeAllocateStatus)
      IF (DeAllocateStatus /= 0) STOP "*** g_ave not deallocated ***"

      IF (my_rank .eq. 0) THEN     !<------------------------------------------------- START POST PROCESSING E WRITING (ONLY MASTER)

!    Check windows load balancing
 
      IF ( 3.d0 * MINVAL(wind_time_red) .lt. MAXVAL(wind_time_red) ) THEN
         WRITE(61,*)  '*********************************************'
         WRITE(61,*)  "***                                       ***"
         WRITE(61,*)  "*** WARNING :  WINDOWS LOAD IS UNBALANCED ***"
         WRITE(61,*)  "*** for a better performance you should   ***"
         WRITE(61,*)  "*** modify the keyword 'windows balance'  ***"
         WRITE(61,*)  "***                                       ***"
         WRITE(61,*)  '*********************************************'
      END IF

!    Scale Tw ***sequentially***

      IF (nwind > 1) THEN

      DO  j = 2, nwind

         DO I = 1, ngrains_per_wind(j)

            Tw_red(I, j) = Tw_red(I, j) * Tw_red(ngrains_per_chunk(j-1) + ref_red(j), j-1)

         END DO

      END DO

      END IF

!     Join Tw from all windows to make T

        IF (nwind > 1) THEN
  
         ALLOCATE ( record_Ttot_bookmark(nwind-1), &
                       STAT = AllocateStatus )
         IF (AllocateStatus /= 0) STOP &
            "*** Not enough memory record_Ttot_bookmark ***"
 
         Ttot_bookmark = 0
         old_cut_point = 0
   
         DO i=1, nwind-1
           T_shift = ngrains_per_chunk(i) 
            DO j=1, ngrains_to_add(i)-1
               dTa = LOG( (Tw_red(T_shift+(j+1), i)-Tw_red(T_shift+j, i))/Egrain1 ) / Egrain1
               dTb = LOG( (Tw_red(j+1, i) - Tw_red(j, i))/Egrain1) / Egrain1

               deriv = ABS(dTa - dTb)
               IF (j == 1 ) THEN
                  deriv_old = deriv
                  cut_point = j
               ELSE IF (deriv < deriv_old) THEN
                  deriv_old = deriv
                  cut_point = j
               END IF

            END DO
   
!           Copy just from first window
            DO j=1, (ngrains_per_wind(i) - ngrains_to_add(i) + cut_point - old_cut_point)
               Ttot(j+Ttot_bookmark) = Tw_red(j+old_cut_point, i)
            END DO
 
            Ttot_bookmark = (j - 1) + Ttot_bookmark
            record_Ttot_bookmark(i) = Ttot_bookmark !For normalization
            old_cut_point = cut_point
 
         END DO

!        Complete with the last window
         DO i=1, (ngrains_per_wind(nwind) - old_cut_point)
            Ttot(i+Ttot_bookmark) = Tw_red(i+old_cut_point, nwind)
         END DO

       ELSE

         DO i=1, ngrains_per_wind(nwind)
            Ttot(i) = Tw_red(i, nwind)
         END DO
 
       END IF
      OPEN (7,STATUS='REPLACE',FILE=chekname,ACTION='WRITE')  ! Chekpoint file

      WRITE(7,9911) AVERSION , ADATE
      WRITE(7,9912) AVERSION , ADATE 
      
      WRITE(7,9002) TITLE1
      WRITE(7,9002) TITLE2 
      CALL DateTime(7)                         ! print output to unit=50
      WRITE(7,9005)

      WRITE (7,99030) cut
      WRITE (7,99023) fname
      WRITE (7,99001) Egrain1 , imax1 , Emax2 , Isize
      WRITE (7,99002)
      DO I = 1 , ngrains + ngrains_to_add(nwind)
         enrg = Emin + (I-1)*Egrain1
         WRITE (7,*) I , enrg , Ttot(I)
      END DO

      CLOSE (UNIT=7) ! Close chekpoint file

       END IF !<------------------------------------------------------------------------------------------ MASTER ONLY
      END IF  !<------------------------------------------------------------------------------------------ End if on checkpoint file

      IF (my_rank.eq.0) THEN !<----------------------------------------------------------------------------------------  MASTER ONLY
    
      IF( chkp .eq. 0) THEN                        !<----------------------------------------- START CONVOLUTION FROM chekpoint FILE

          OPEN (7,STATUS='OLD',FILE=chekname,ACTION='READ',IOSTAT=istat)  ! Chekpoint file

      IF ( istat .NE. 0 ) THEN
        WRITE(*,*) '********************************************'
        WRITE(*,*) 'FATAL: problem opening chekpoint file'
        WRITE(*,*) '********************************************'
        WRITE(4,*) '********************************************'
        WRITE(4,*) 'FATAL: problem opening chekpoint file'
        WRITE(4,*) '********************************************'
        STOP
      END IF
      
      dumdum='test'
      DO WHILE ( dumdum .NE. cut )
         READ(7,99030) dumdum
      END DO

      itest = 0
      READ (7,*) fnameT
      READ (7,*) Egrain1T , imax1T , Emax2T , IsizeT
     
      IF ( fnameT .NE. fname)     itest = 1
      IF ( Egrain1T .NE. Egrain1) itest = itest + 10
      IF ( imax1T .NE. imax1)     itest = itest + 100
      IF ( Emax2T .NE. Emax2)     itest = itest + 1000
      IF ( IsizeT .NE. Isize)     itest = itest + 10000
      IF ( itest .NE. 0 ) THEN
        WRITE(*,*) '********************************************'
        WRITE(*,*) 'FATAL: chekpoint PARAMETERS DO NOT MATCH'
        WRITE(*,*) '       itest =',itest
        WRITE(*,*) '********************************************'
        WRITE(4,*) '********************************************'
        WRITE(4,*) 'FATAL: problem opening chekpoint file'
        WRITE(4,*) '       itest =',itest
        WRITE(*,*) '********************************************'
        STOP
      END IF
      
      DO iunit = 3, 5
       WRITE(iunit,*) 
       WRITE(iunit,*)' CALCULATED USING chekpoint FILE', chekname
       WRITE(iunit,*) 
      END DO
      
      READ(7,*) dumdum
      READ(7,*) dumdum
      
      DO ig = 1 , ngrains + ngrains_to_add(nwind)
        READ(7,*)  I , enrg , Ttot(I)
      END DO
      
      CLOSE (UNIT=7)                                                 ! Close chekpoint file
      
      END IF                                                          !<------------- End if on checkpoint file
      
      zpeTotal = zpe
!
!     Convolution of separable modes with non-separable via
!       the Beyer-Swinehart/Stein-Rabinovitch Algorithm
!
      IF ( Nsep .GT. 0 ) THEN       ! If separable modes are to be included
      
      DO I=ns+1, ns + Nsep
          IF((IDOF(I).EQ.'HRA').OR.(IDOF(I).EQ.'HRB').OR. & !NEW
               (IDOF(I).EQ.'HRC')) THEN !NEW
            !CALL SHRLEV(Ttot,AT,Egrain1,ngrains,HRB(I),HRVo(I),NG(I),1, &
            CALL SHRLEV(Egrain1,ngrains+ngrains_to_add(nwind),HRB(I),HRVo(I),NG(I),1, &
              IMAX(I),HRzpe(I), 'Vhrd1','Bhrd1',0.0d0,0.0d0,NG(I))  ! Symmetrical hindered rotor
            zpeTotal = zpeTotal + HRzpe(I)
          ELSEIF (IDOF(I).EQ.'HRD') THEN
           DO J=1, NVV(I)
              CVt(J)=CV(I,J)
           END DO
           DO J=1, NBB(I)
              CBt(J)=CB(I,J)
           END DO
           !CALL UHRLEV(Ttot,AT,Egrain1,ngrains,NBB(I),NVV(I),CBt,CVt, &
           CALL UHRLEV(Egrain1,ngrains+ngrains_to_add(nwind),NBB(I),NVV(I),CBt,CVt, &
               NGV(I),NGB(I),IMAX(I),HRzpe(I),Vhr(I),Bhr(I) & 
                ,Phav(I),Phab(I),NG(I))                               ! General, unsymmetrical hindered rotor
            zpeTotal = zpeTotal + HRzpe(I)
         ELSEIF (IDOF(I).EQ.'ROT') THEN
            B=16.85763d0/RI(I)
            NSYMM=NRSN(I)
            NDIM=IDIM(I)
            !CALL CROTLEV(Ttot,Egrain1,ngrains,B,NDIM,NSYMM)
            CALL CROTLEV(Egrain1,ngrains+ngrains_to_add(nwind),B,NDIM,NSYMM)
         ELSEIF (IDOF(I).EQ.'QRO') THEN
            B=16.85763d0/RI(I)
            NSYMM=NRSN(I)
            NDIM=IDIM(I)
            !CALL ROTLEV(Ttot,AT,Egrain1,ngrains,B,NDIM,NSYMM)
            CALL ROTLEV(Egrain1,ngrains+ngrains_to_add(nwind),B,NDIM,NSYMM)
         ELSEIF ( IDOF(I).EQ.'TOP' ) THEN !NEW                         ! ******** Symmetric Top *********
            NSYMM = NG(I)              ! Symmetry number !NEW
            !CALL STOPLEV(T,AT,Egrain1,ngrains,B2,B1,NSYMM) !NEW
            CALL STOPLEV(Egrain1,ngrains,B2,B1,NSYMM) !NEW
         ELSEIF (IDOF(I).EQ.'KRO') THEN
            B=16.85763d0/RI(I)
            NSYMM=NRSN(I)
            JK=IDIM(I)
            !CALL KROTLEV(Ttot,AT,Egrain1,ngrains,B,JK,NSYMM)
            CALL KROTLEV(Egrain1,ngrains+ngrains_to_add(nwind),B,JK,NSYMM)
         ELSEIF (IDOF(I).EQ.'BOX') THEN
            B=RI(I)
            NDIM=IDIM(I)
            !CALL BOXLEV(Ttot,AT,Egrain1,ngrains,B,NDIM)
            CALL BOXLEV(Egrain1,ngrains+ngrains_to_add(nwind),B,NDIM)
            zpeTotal = zpeTotal + boxzpe !NEW
         ELSEIF (IDOF(I).EQ.'VIB') THEN
            WE=RI(I)
            XE=RJ(I)
            !CALL MORLEV(Ttot,AT,Egrain1,ngrains,WE,XE,zpeMorse)
            CALL MORLEV(Egrain1,ngrains+ngrains_to_add(nwind),WE,XE,zpeMorse)
            zpeTotal = zpeTotal + zpeMorse !NEW
         ELSE
            WRITE(*,*) '************************************'
            WRITE(*,*) 'FATAL: MISTAKES at convolution step'
            WRITE(*,*) 'Degree of freedom #',I
            WRITE(*,*) '************************************'
            WRITE(4,*) '************************************'
            WRITE(4,*) 'FATAL: MISTAKES at convolution step'
            WRITE(4,*) 'Degree of freedom #',I
            WRITE(4,*) '************************************'
            STOP
         ENDIF
      END DO
      ENDIF  ! separable modes

      DO iunit = 3, 5 !NEW
         WRITE (iunit,99013) zpeTotal !NEW 
      END DO !NEW

!
!     Compute sums and densities of states
!
      SS = 0
      DS = 0
      DS(1)=Ttot(1)/Egrain1
      SS(1)=Ttot(1)

         DO I=2, ngrains
            SS(I)=SS(I-1) + Ttot(I)
            DS(I)=Ttot(I)/Egrain1
         END DO

! --------------------------------------------------------------------------
!     Write adensum.out output file
!

        ngrains = ngrains - 1   !!! the boundary of interest

      WRITE(4,9017) 
      DO I=1, ngrains
       enrg = Emin + (I-1)*Egrain1
       WRITE(4,9010) I , enrg , DS(I), SS(I)
      END DO 
! --------------------------------------------------------------------------
!
!      Write to ___.dens file (UNIT=3) 
! 

      WRITE (3,99030) cut                    ! End of data summary block in ____.dens file
      WRITE (5,99030) cut                    ! End of data summary block in ____.qvib file

      WRITE (3,*) fname(1:lenstr(fname))
      WRITE (3,*) TITLE1
      WRITE (3,99001) Egrain1 , imax1 , Emax2 , Isize , Viblo , zpeTotal

      WRITE (3,99010)
      DO 200 I = 1 , imax1
         enrg = Emin + (I-1)*Egrain1
         WRITE (3,99014) I , enrg , DS(I), SS(I)
 200  CONTINUE
 
      DO 300 I = imax1 + 1 , Isize
         enrg = Emin + (I-imax1-1)*Egrain2
         j = INT(enrg/Egrain1) + 1
         WRITE (3,99014) I , enrg , DS(j), SS(j)
 300  CONTINUE

      WRITE (5,*) fname(1:lenstr(fname))
      WRITE (5,*) TITLE1
      WRITE (5,99901) Egrain1 , Emax2 , zpeTotal !, KEYWORD
99901 FORMAT (3(3x,f10.2),3x,A6)


      Ntx = 150        ! Number of temperatures to be passed to THERMO
      WRITE(5,*) Ntx
      WRITE(5,*)
      WRITE(5,*) "==================================================="
      WRITE(5,*) "     VIBRATIONAL PARTITION FUNCTION"
      WRITE(5,*) "==================================================="
      WRITE(5,*)
      WRITE(5,*) "   INDEX     T(K)       Q       C/R      H/R      S/R"

!        TT=25.0d0
        TT=1.0d0
        DO I = 1 , Ntx           ! Number of temperatures to be passed to THERMO
           Q0 = 0.0d0
           Q1 = 0.0d0
           Q2 = 0.0d0
           DO J=1, ngrains 
             enrg = (J-1)*Egrain1
             Ered = 1.4387752D+00*enrg/TT
             fx = DS(J)*exp(-Ered)*Egrain1
             Q0 = Q0 + Fx
             Q1 = Q1 + Ered*Fx
             Q2 = Q2 + Ered*Ered*Fx
           END DO   
!            Thermo formulae from Pitzer and Gwinn, J. Chem. Phys. 10, 428 (1942)
           Cx = Q2/Q0 - (Q1/Q0)**2
           Hx = Q1/Q0
           Sx = (Hx+log(Q0))
           WRITE(5,9011) I, TT, Q0 , Cx , Hx , Sx
!           TT = TT*1.050824
           TT = TT*1.05631305d+00  ! for Ntx = 150
        END DO   

      WRITE(5,*)

!
!       End of writing
!
! --------------------------------------------------------------------------


      END IF   !<----------------------------------------------------------------------- END POST PROCESSING E WRITING (ONLY MASTER)
    

9911  FORMAT ( &
      '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
      %%%%%%%%%%%%%'//,&
      '                      paradensum-',A7,/ &
      '                         ',A8,//, &
      ' Authors: adensum:',/ &   
      '          John R. Barker    jrbarker@umich.ed',/ &
      '          Thanh Lam Nguyen  nguyenlt@cm.utexas.edu',// &
      '          paradensum: ',/ &
      '          Chiara Aieta      chiara.aieta@unimi.it',/ &
      '          Fabio Gabas       fabio.gabas@unimi.it',/ &
      '          Michele Ceotto    michele.ceotto@unimi.it',// &
      '       http://clasp-research.engin.umich.edu/multiwell/',// &
      '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
      %%%%%%%%%%%%%',/)
9912  FORMAT('Suggested Literature Citations:'//4x,&
      'a) MultiWell-',A6,' Software, ',A8,', designed and'/7x,&
      'maintained by J.R. Barker with contributor T.L. Nguyen,'/7x,&
      'J.F. Stanton, C. Aieta, M. Ceotto, F. Gabas, T.J.D. Kumar'/7x,&
      'C.G.L. Li, L.L. Lohr, A. Maranzana, N.F. Ortiz, J.M. Preses'/7x,&
      'J.M. Simmie, J.A. Sonk, and P.J. Stimac, University of Michigan,'/7x,&
      'Ann Arbor, MI; http://clasp-research.engin.umich.edu/multiwell/.'//4x,&
      'b) John R. Barker, Int. J. Chem. Kinetics, 33, 232-45'/7x,&
      '(2001).'//4x,&
      'c) John R. Barker, Int. J. Chem. Kinetics, 41, 748-763'/7x,&
      '(2009).'//4x,&
      'c) M. Basire, P. Parneix, and F. Calvo, J. Chem. Phys.'/7x,&
      ' 129, 081101 (2008).'//4x,&
      'd) F. Wang and D. P. Landau, Phys. Rev. Letters'/7x,& 
      '86, 2050 (2001).'//4x,&
      'e) Thanh Lam Nguyen and John R. Barker, J. Phys. Chem. A.,'7x,&
      '114, 3718–3730 (2010)'//4x,&
      'f) C. Aieta, F. Gabas, and M. Ceotto, J. Phys. Chem. A,'7x,& 
      '120 4853 (2016).',&
      //,&
      '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
      %%%%%%%%%%%%%'/)
9002  FORMAT(A100)
9003  FORMAT(I3,' Vibrational Frequencies (cm-1)  ',/,&
       '  We; E = We*(v+1/2) +...; Wo: E = W0*v +...;  0<->1 fund.  ',/,&
       '  De: BDE for Separable Morse oscillator [0.0 if non-Morse]',//,&
       '  No.     We         Wo      0-1_Fund.   De(cm-1)' )
9004  FORMAT( I5,3(1x, F10.3),1x,F10.1 )
9042  FORMAT( I3, 100( 1x , f9.3 ) )
9043  FORMAT( 100(7x,I3) )
9005  FORMAT(/)
9006  FORMAT('Anharmonicity Matrix (cm-1)')
!9008  FORMAT(/'Monte Carlo trials: ', 1pe8.1,' (',A6,')',/, &
!             '  iterations = ',I2,/, &
!             '          p0 = ',0pf7.3,/, &
!             '          f  = ',0pf7.3,/)
9010  FORMAT( 2x,I6,2x, F10.2,5(2x,1pe10.3) )
9011  FORMAT( 2x,I6,2x, F13.5,5(2x,1pe12.5) )
9013  FORMAT(/,' Ave:', 3(1x, F10.3)/ &
        'Anharmonic ZPE (omitting VPT2 G(0) term) = ', f10.3,' (cm-1)'/)
!9016  FORMAT( '             steps up or down:',f10.3,' %  ',/, &
!             '               steps in place:',f10.3,' %  ',/, &
!             '           steps out of range:',f10.3,' %  ',/, &
!             '   average non-zero histogram: ',f10.0,' +/-', f8.0,/, &
!             '  average un-signed deviation:',0pf10.2,' %  ',/, &
!             '       max positive deviation:',0pf10.2,' %  ',/, &
!             '        max negative devation:',0pf10.2,' %',//)
9017  FORMAT( &
       '     No.        cm-1  states/cm-1  SumStates')
9020  FORMAT(' Iteration =',I3,'/',I2,5x,'f = ',1pe10.2, 5x, &
                      'p = ',1pe10.3)
9030     FORMAT(I3,2X,'Hind.Rot',2x,'Freq(har in cm-1)=',F8.2,2x, &
          'Mom(amu.A**2)=',F8.3,2x,'fold=',I2,2x,'symm=',I2,2x, &
          'B=',F8.4,' cm-1',2x,'Uo=',F7.1,' cm-1')
9031     FORMAT (I3,2X,'Non-Rigid HinRotor:',2x,'Rotor symm. =',I2)
9032     FORMAT (8X,A5,1x,' ; symmV =',I2,2x,'; Phase(rad.) = ', &
      F8.4,2x,'; Coeff. =',20(F10.4,1x))

9033     FORMAT (8X,A5,1x,' ; symmB =',I2,2x,'; Phase(rad.) = ', &
      F8.4,2x,'; Coeff. =',20(F10.4,1x))
99001 FORMAT (1x,f10.1,1x,I10,1x,f10.1,1x,I10,1x,f10.1,1x,f10.2)
99002 FORMAT (/,'T(I) array for coupled degrees of freedom:',/, &
             '        I      (cm-1)   T(I)')
99010 FORMAT ('       No.     (cm-1)   Density   Sum')
99013 FORMAT (/,'TOTAL ZPEanh: ', f10.3,' cm-1 for all coupled & ', &
       'separable D.O.F.; omits the reaction coordinate and ', &
       'Go from VPT2.',/)
99014 FORMAT (i10,1x,f10.1,2(1x,1pe10.3))
99023 FORMAT (/A10)
99030 FORMAT (A46)
99041 FORMAT (I3,2X,'Quantized Rotor:  I(Amu.A**2) = ', &
           F8.3,' ;  Sym. = ',I1,' ;  Dim. = ',I2,' ;  B(cm-1) = ',F8.3)
99042 FORMAT (I3,2X,'K-Rotor:  I(Amu.A**2) = ',F8.3,' ;  Sym. = ',I2, &
                ' ;  Quantum No. J = ',I2,' ;  B(cm-1) = ',F7.3)
99043 FORMAT (I3,2X,'Particle-in-Box:  Freq(cm-1)  = ',F7.2, &
                   ' ;  Dim. = ',I2)
99044 FORMAT (I3,2X,'Classical Rotor:  I(Amu.A**2) = ', &
           F8.3,' ;  Sym. = ',I1,' ;  Dim. = ',I2,' ;  B(cm-1) = ', &
           F8.3)
99045 FORMAT (I3,2X,'Vibration:',8x,'Freq(cm-1)  = ',F7.2, &
      ' ;  Anharm. = ', F8.3)
99050 FORMAT (//,' SEPARABLE DEGREES OF FREEDOM',/)
99051 FORMAT (/' --- END INPUT ---')
99052 FORMAT ('Number of windows: ', I3)
99053 FORMAT ('Number of walkers in each window: ', I3)
99054 FORMAT ('Overlap: ', F6.2,'%')
99055 FORMAT ('Flatness criterion: ', F5.3)
99056 FORMAT ('Random number generator seed modifier: ', I3)
99905 FORMAT (I5,2X,'Symm-Top',2x,'2D-Moment =',F8.2,' amua, B2 =',F8.4, &
            ' cm-1  ;  1D-Moment =',F8.2,' amua, B1 =',F8.4,' cm-1',2x,  &
            ';  Symm. No. = ',I1)

      call MPI_FINALIZE(err) 

      END
