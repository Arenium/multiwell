c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    
!    LICENSE NOTICE
!
!    bdens: sums and densities of coupled anharmonic vibrations
!    Copyright (C) 2016 Collin G. L. Li, Thanh Lam Nguyen, and John R. Barker
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License: <http://www.gnu.org/licenses/>.
!
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c                                                                 
c    PROGRAM bdens
c    by John R. Barker, Thanh Lam Nguyen, and Collin G. L. Li    
c                                                                  
c  ***Direct Count and Wang-Landau Algorithm for Densities of States***
c
c    Contact:
c    John R. Barker   (email: jrbarker@umich.edu)
c    Department of Atmospheric, Oceanic, and Space Sciences
c    College of Engineering
c    University of Michigan
c    Ann Arbor, MI 48109
c
c    http://clasp-research.engin.umich.edu/multiwell
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c    Literature Citations:                                                    
c                                                                  
c    1. F. Wang and D. P. Landau, Phys. Rev. Letters 86, 2050-2053 (2001).           
c                                                                  
c    2. M. Basire, P. Parneix, and F. Calvo, J. Chem.Phys. 129, 081101 (2008).    
c                                                                  
c    3. T. L. Nguyen and J. R. Barker, J. Phys. Chem. A., 114, 3718–3730 
c         (2010). DOI: 10.1021/jp100132s    
c                                                                  
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      PROGRAM bdens

      USE vibcalcs
      USE SEPSUBS
      USE HRMOD
      IMPLICIT NONE
      
      CHARACTER(len=7) AVERSION 
      CHARACTER(len=8) ADATE 
      PARAMETER ( AVERSION='2023', ADATE='Mar 2023' )
 
      CHARACTER(46) cut , dumdum
      PARAMETER (cut='**************INPUT DATA SUMMARY**************')
      
      CHARACTER(10) fname , fnameT
      
      CHARACTER(9) chekpoint
      
      INTEGER istat , itest , iunit , ntop

      INTEGER(4) :: i, j, k, MaxTrials,
     & idum, l , nstart , Nacc, Nrej, Nout,
     & iters, ngrains , nold , nnew , ngrains2 ,
     & ndmax , nh , ig , nnt, nvmax , ngrainswl
      INTEGER(4), PARAMETER :: KIN = 2       ! input UNIT for READ statements
      
      REAL(8)         Eu, trialmax, ttmax,
     &                 Emax, Estart, energy, 
     &                 RAN1, D(100),
     &                 enrg , Etop
      REAL(8) fxx , Ered, Q0, Q1, Q2, Cx, Hx, Sx

      Real(8) Eswitch                     !Energy for switch point between bdens and recursive direct counting
      Real(8) MVAL, MVALT                           !Energy for switch point if manually set

      REAL(8) H(20003), g(20003), p , p0 , f , fo , test , gmax ,
     &                 acc , hmax , hmin , have, hav2, y , var
 
      CHARACTER(100) TITLE1 , TITLE2
      CHARACTER(2) WW
      CHARACTER(6) :: KEYWORD                     ! designates no. of Wang-Landau stochastic trials; assumed to be upper case
      CHARACTER(6) :: LTMODE                      ! Eswitch calculation mode: AUTO, MAN
      CHARACTER(6) :: LTMODET                     ! Eswitch calculation mode: AUTO, MAN

      CHARACTER(4) VROTIN , VROTOUT
      
      REAL(8) HRfreq(MaxDof), HRVo(MaxDof)  , HRI(MaxDof) ,HRBb(MaxDof), 
     &        HRzp(MaxDof) , DUM1 , DUM2 , WEh, XEh ,
     &        zpetot , morzpe , boxzpe
      INTEGER(4) Nsep, NVV(100), NBB(100), 
     &     NGV(100), NGB(100)
      INTEGER, PARAMETER :: IMAX=501              ! hindered rotor grid points and maximum number of eigenvalues
      REAL(8) T(20003), AT(20003), DS(20003), SS(20003)
      REAL(8) CAD(20003), CDO(20003), SUMAD, SUMDO, CDOAD
!      REAL(8) T(30003), AT(30003), DS(30003), SS(30003)
!      REAL(8) CAD(30003), CDO(30003), SUMAD, SUMDO, CDOAD
      
      REAL(8) XE
      
      REAL(8) BB , BG , FAC , RG
      INTEGER N , NR , nvibstest

      INTEGER(4) NY, NZ
      INTEGER(4) nvmaxxyz, JMAX 
      REAL(8) energyxyz
     
      REAL(8) Br
      INTEGER(4) lenstr

      REAL(8) RI(100) , RJ(100)
      INTEGER(4) NRSN(100), IDIM(100)
      INTEGER(4) NDIM, NSYMM, JK
 
      REAL(8) Egrain1 , Egrain2 , Emax1 , Emax2
      INTEGER(4) imax1 , imax2 , Isize
      
      REAL(8) Egrain1T , Emax2T 
      INTEGER(4) imax1T , IsizeT
      
      REAL(8) X1 , X2
      INTEGER(4) ii

      INTEGER(4) ntest , kk , nn
      
      INTEGER nTmax, nT , Ntx

      REAL(8) TT , su , beta0 , RR
      PARAMETER (RR=0.695038916D0)
      
      COMMON/I/ idum
      COMMON/TBAR/ nTmax, nT(20003), Egrain1
      
      EXTERNAL energy , energyxyz, RAN1 , ndmax, nvmax, nvmaxxyz
      
      SAVE
      
      OPEN (2,FILE='bdens.dat',STATUS='OLD')
      READ (2,*) fname

      OPEN (3,STATUS='UNKNOWN',FILE=fname(1:lenstr(fname))//'.dens')  ! Succinct output for MultiWell input
      OPEN (4,FILE='bdens.out',STATUS='UNKNOWN')
      OPEN (5,FILE=fname(1:lenstr(fname))//'.qvib',STATUS='UNKNOWN')
      OPEN (12,STATUS='UNKNOWN',FILE='bdens.lev')   ! Energy levels for vibrations, rotations, and hindered rotors
      
      READ (2,9002)  TITLE1
      READ (2,9002)  TITLE2

      WRITE(3,99030) cut                                  ! Start of data summary block in ____.dens file
      WRITE(5,99030) cut                                  ! Start of data summary block in ____.qvb file

      DO iunit = 3 , 5
        WRITE(iunit,9911) AVERSION , ADATE
        WRITE(iunit,9912) AVERSION , ADATE 
        WRITE(iunit,99023) fname
        WRITE(iunit,9002) TITLE1 
        WRITE(iunit,9002) TITLE2 
        CALL DateTime(iunit)                                ! print output to unit=iunit
      END DO
      
      WRITE (4,9005)

      READ(2,*)  ns , NY, NZ, WW
      CALL ucase ( WW )

      CALL READ_WXYZ(KIN , NY, NZ, WW)                      ! read frequencies and anharmonicities
      
      READ(2,*) NSEP , VROTIN                               ! No. of separable DOF and keyword for rotation units
      CALL ucase ( VROTIN )
      VROTOUT = 'AMUA'

         ntop = 0            ! for counting the number of symmetric tops (should never be >1)      DO I=ns+1, ns+NSEP
      DO I=ns+1, ns+Nsep
         READ(2,*) MODE(I), IDOF(I), DUM1, DUM2, NG(I)
         CALL ucase ( IDOF(I) )
         IF (IDOF(I).EQ.'HRA') THEN
            HRfreq(I)=DUM1
            HRI(I)=DUM2
            HRBb(I)=16.85763d0/HRI(I)
            HRVo(I)=((HRfreq(I)/NG(I))**2)/HRBb(I)
         ELSEIF (IDOF(I).EQ.'HRB') THEN
            HRfreq(I)=DUM1
            HRVo(I)=DUM2
            HRBb(I)=((HRfreq(I)/NG(I))**2)/HRVo(I)
            HRI(I)=16.85763d0/HRBb(I)
         ELSEIF (IDOF(I).EQ.'HRC') THEN
            HRI(I)=DUM1
            HRVo(I)=DUM2
            HRBb(I)=16.85763d0/HRI(I)
            HRfreq(I)=NG(I)*sqrt(HRVo(I)*HRBb(I))
         ELSEIF (IDOF(I).EQ.'HRD') THEN
            NVV(I)=INT(DUM1)
            NBB(I)=INT(DUM2)
            READ(2,*) Vhr(I), NGV(I), Phav(I), (CV(I,J), J=1, NVV(I))
            READ(2,*) Bhr(I), NGB(I), Phab(I), (CB(I,J), J=1, NBB(I))
            CALL ucase ( Vhr(I) )
            CALL ucase ( Bhr(I) )
         ELSEIF (IDOF(I).EQ.'ROT') THEN
       IF ( NG(I) .GT. 2 ) THEN
         WRITE(*,*) '***FATAL ERROR***: execution terminated'
         WRITE(*,*) 'Classical rotor dimension can be 1 or 2, only.'
         WRITE(*,*) 'Classical spherical tops must be computed using'
         WRITE(*,*) 'two separable rotors: 1-D and 2-D with same B.'
         WRITE(3,*) '***FATAL ERROR***: execution terminated'
         WRITE(3,*) 'Classical rotor dimension can be 1 or 2, only.'
         WRITE(3,*) 'Classical spherical tops must be computed using'
         WRITE(3,*) 'two separable rotors: 1-D and 2-D with same B.'
         STOP
       ENDIF
            RI(I) = DUM1
            NRSN(I) = INT(DUM2)
            IDIM(I) = NG(I)
            CALL rotunits( VROTIN , RI(I) , VROTOUT)
         ELSEIF (IDOF(I).EQ.'QRO') THEN
            RI(I) = DUM1
            NRSN(I) = INT(DUM2)
            IDIM(I) = NG(I)
            CALL rotunits( VROTIN , RI(I) , VROTOUT)
         ELSEIF ( IDOF(I).EQ.'TOP' ) THEN                               ! Symmetric Top with 2D and 1D rot constants B2 and B1
            ntop = ntop + 1
            IF ( ntop .GT. 1) THEN
              write (*,*) 'FATAL: only one symmetric top (TOP) allowed'
              write (3,*) 'FATAL: only one symmetric top (TOP) allowed'
              write (4,*) 'FATAL: only one symmetric top (TOP) allowed'
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
            write(*,*) '****FATAL: KRO type of d.o.f. not implemented'
            write(3,*) '****FATAL: KRO type of d.o.f. not implemented'
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
 
      READ (2,*) Egrain1 , Imax1 , Isize , Emax2 , KEYWORD, LTMODE,MVAL
      CALL ucase ( LTMODE )
      CALL ucase ( KEYWORD )

      IF (Emax2/Egrain1 .GT. NMAX) THEN
        WRITE(*,*) 'FATAL: **Grain Size Too Small: too many grains****'
        WRITE(4,*) 'FATAL: **Grain Size Too Small: too many grains****'
        STOP
      ENDIF
        IF ( Imax1 .GT. Isize-2 ) THEN
         WRITE (*,*) '*** FATAL:  Imax1 must be ≤ (Isize-2) ***'
         WRITE (4,*) '*** FATAL:  Imax1 must be ≤ (Isize-2) ***'
         STOP
        ENDIF      
      
      READ(2,*) chekpoint
      CALL UCASE ( chekpoint )
      
      CLOSE (UNIT=2)
c
c   ******* END INPUT *******
c
      Emax = Emax2                                ! nominal top of Wang-Landou energy window
         
      CALL convib( WW )                           ! convert input frequencies to other forms
 
      Emax1 = Egrain1*(imax1-1)
      Egrain2 = Emax2/(Isize-imax1-1)
      imax2 = Isize - imax1
 
      JMAX = 1 + INT(Emax2/Egrain1)
      JMAX = MIN(JMAX,NMAX)

      WRITE(4,9003) ns
      WRITE(3,9003) ns
      WRITE(5,9003) ns
      nvibstest = 0
      DO i = 1 , ns
        D(i) = 0.0
        IF (xa(i,i) .LT. 0.0) D(i) = -0.25D+00*w0(i)**2/xa(i,i) 		! Bond dissociation energy
        IF ( wf(i) .LE. 0.0 .AND. nvibstest .EQ. 0 ) nvibstest = i   	! Test for unphysical fundamental frequency
        WRITE(4,9004) i, wa(i), w0(i), wf(i), D(i)
        WRITE(3,9004) i, wa(i), w0(i), wf(i), D(i)
        WRITE(5,9004) i, wa(i), w0(i), wf(i), D(i)
      END DO

      WRITE(4, 9013) ave, av0, avf, zpe                         ! zero point energy and average frequencies

      WRITE(3,9006)
      WRITE(4,9006)
      WRITE(5,9006)
      WRITE(3,9043) (i, i=1,ns)
      WRITE(4,9043) (i, i=1,ns)
      WRITE(5,9043) (i, i=1,ns)
      DO j = 1 , ns
        WRITE(3,9042) j, (xa(j,i), i=1,ns)
        WRITE(4,9042) j, (xa(j,i), i=1,ns)
        WRITE(5,9042) j, (xa(j,i), i=1,ns)
      END DO

      zpetot = zpe     ! at this point, zpetot includes coupled modes, but nothing else
      
      IF ( nvibstest .GT. 0 ) then		! At least one fundamental frequency is un-physical
         WRITE(4, 9991) nvibstest
         WRITE(6, 9991) nvibstest
         STOP
      ENDIF
9991  FORMAT(//,'*** FATAL ** unphysical fundamental frequency, i=', I2)

      Eswitch = feswitch( Emax, Egrain1, LTMODE, MVAL,KEYWORD)       ! Emin for bdens calculation
c
c    Write out separable degrees of freedom
c

C      write(4,*)
C      write(3,*)
      write(4,99050)
      write(3,99050)
      write(5,99050)
      DO I=ns+1, ns+Nsep
            IF (IDOF(I).EQ.'HRA') THEN
             DO iunit = 3, 5
             write(iunit,9030) MODE(I), HRfreq(I), HRI(I), NG(I), NG(I), 
     &               HRBb(I), HRVo(I)
             END DO
            ELSEIF (IDOF(I).EQ.'HRB') THEN
             DO iunit = 3, 5
             write(iunit,9030) MODE(I), HRfreq(I), HRI(I), NG(I), NG(I), 
     &               HRBb(I), HRVo(I)
             END DO
            ELSEIF (IDOF(I).EQ.'HRC') THEN
             DO iunit = 3, 5
             write(iunit,9030) MODE(I), HRfreq(I), HRI(I), NG(I), NG(I), 
     &               HRBb(I), HRVo(I)
             END DO
            ELSEIF (IDOF(I).EQ.'HRD') THEN
             DO iunit = 3, 5
                write(iunit,9031) MODE(I), NG(I)
                write(iunit,9032) Vhr(I), NGV(I), Phav(I),
     &                  (CV(I,J), J=1, NVV(I))
             END DO
             IF(Bhr(I).EQ.'BHRD1') THEN
               DO iunit = 3, 5
                  write(iunit,9033) Bhr(I), NGB(I), Phab(I),
     &                  (CB(I,J), J=1, NBB(I))
               END DO
              ELSE
               DO iunit = 3, 5
                  write(iunit,9033) Bhr(I), NGB(I), Phab(I),
     &                  (CB(I,J), J=1, NBB(I))
               END DO  
              ENDIF
            ELSEIF (IDOF(I).EQ.'ROT') THEN
                  Br = 16.85763d0 / RI(I)
             DO iunit = 3, 5
                  write(iunit,99044) MODE(I), RI(I),NRSN(I), IDIM(I), Br
             END DO
            ELSEIF (IDOF(I).EQ.'QRO') THEN
                  Br = 16.85763d0 / RI(I)
             DO iunit = 3, 5
                  write(iunit,99041) MODE(I), RI(I),NRSN(I), IDIM(I), Br
             END DO
            ELSEIF ( IDOF(I).EQ.'TOP' ) THEN                               ! Symmetric Top with 2D and 1D rot constants B2 and B1
               DO iunit = 3, 5
                 WRITE (iunit,99905) MODE(I) , X2 , B2 , X1 , B1 , NG(I)
               END DO
            ELSEIF (IDOF(I).EQ.'KRO') THEN
                  Br = 16.85763d0 / RI(I)
             DO iunit = 3, 5
                  write(iunit,99042) MODE(I), RI(I),NRSN(I), IDIM(I), Br
             END DO
            ELSEIF (IDOF(I).EQ.'BOX') THEN
             DO iunit = 3, 5
                  write(iunit,99043) MODE(I), RI(I), IDIM(I)
             END DO                  
            ELSEIF (IDOF(I).EQ.'VIB') THEN
             DO iunit = 3, 5
                  write(iunit,99045) MODE(I), RI(I), RJ(I)
             END DO                  
            ENDIF

      ENDDO

      DO iunit = 3, 5, 1
        WRITE(iunit, *) ' '
        WRITE(iunit, 99052) LTMODE, Eswitch
        write(iunit,99051)
      END DO

c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------

      SELECT CASE ( chekpoint .EQ. 'CHEKSTART')
      
      CASE (.FALSE.)    ! does NOT start from saved chekpoint file
      
c-----------------------------------------------------------------------------------------
c
c	recursive direct count algorithm	
c
c-----------------------------------------------------------------------------------------
c     Emax2 in bdens is Emax in recursive direct count
c     Egrain1 in bdens is DELE in recursive direct count
c     wa, ns, xa, ya, za, NY and NZ are the same in both programs
c     
c         Recursive direct count calculate up to 1.1*Eswitch because direct 
c         count algorithm is inaccurate near high end of set range

      nTmax = INT( 1.1*Eswitch/Egrain1 ) + 1
      DO I = 1, nTmax                                                ! initialize nT
        nT(I) = 0
      END DO

      WRITE(*,*)'Be patient: Direct counting (recursive algorithm)'
      CALL CALEVN( 1.1*Eswitch, wa, ns, xa, ya, za, NY, NZ )            ! 1.1*Eswitch  

      JMAX=INT(Eswitch/Egrain1) + 1
      DO I = 1, JMAX+1
        T(I) = DBLE( nT(I) )                                         ! T(I) contains direct count
      END DO

      WRITE(4,*) '   '
      WRITE(4,*) 'Finished Direct Count '
      CALL DateTime(4)                                      ! print output to unit=4
      
      IF ( Eswitch .LT. Emax2 ) THEN                        ! then do Wang-Landau
      WRITE(*,*) 'Starting Wang-Landau counting'
c-----------------------------------------------------------------------------------------
c         WANG-LANDAU PARAMETERS
c
      p0 = 1.0/ns                                           ! Scaled according to number of vibrations
      IF ( p0 .GT. 0.25d+00 ) p0 = 0.25d+00
      iters = 21                                            ! Always
      Etop = Emax2 + 10*Egrain1                             ! a little higher than the range of interest
      ngrainswl = INT( (Etop-0.9*Eswitch)/Egrain1 + 1 )   
      fo = 1.0d+00
c-----------------------------------------------------------------------------------------
c   Trials per energy grain
c
      ttmax = trialmax( KEYWORD , ngrainswl )
      MaxTrials = INT( ttmax )
c-----------------------------------------------------------------------------------------

      idum = -1147483647                                  ! Random number seed
      y = RAN1( idum )
      f  = fo
c
c     Initialize quantum numbers
c
      ntest = 1
      DO WHILE ( ntest .EQ. 1 )
       Estart = Eswitch + RAN1( idum )*(Etop-Eswitch) 
        DO l = 1, ns
          DO k=l+1, ns
               nvold(k)=0
          ENDDO
          IF((NY.EQ.0).AND.(NZ.EQ.0)) THEN
            Eu=Estart-energy( ns, wa, xa, nvold )
            IF(Eu .NE. 0.0d0) THEN
               nstart=nvmax( ns, wa, xa, nvold, l, Eu )
               nvold(l) = nstart*RAN1( idum )
            ELSE
               nvold(l) = 0
            ENDIF
          ELSE
            Eu=Estart-energyxyz( ns , wa , xa , ya, za, nvold )       ! initializg quantum numbers
            IF(Eu .NE. 0.0d0) THEN
               nstart = nvmaxxyz( ns, wa, xa, ya, za, nvold, l, Eu ) 
               nvold(l) = nstart*RAN1( idum )
            ELSE
               nvold(l) = 0
            ENDIF
          ENDIF     
        ENDDO

      IF((NY.EQ.0).AND.(NZ.EQ.0)) THEN
         CALL ckderiv( ns, nvold, wa, xa, ntest )               ! check derivatives for X-matrix
      ELSE 
         CALL ckderivxyz( ns, nvold, wa, xa, ya, za, ntest )    ! check derivatives for XYZ-matrix
      ENDIF 

      END DO                                                    ! END of DO WHILE
      
      IF((NY.EQ.0).AND.(NZ.EQ.0)) THEN
          enrg = energy( ns, wa, xa, nvold )
      ELSE
             enrg = energyxyz( ns , wa , xa , ya, za, nvold )
             ENDIF
        nold = IDNINT( (enrg - 0.9*Eswitch)/Egrain1 ) + 1              ! energy grain number: nearest INTEGER(4)
c-----------------------------------------------------------------------------------------
c     Start Iterations
c-----------------------------------------------------------------------------------------
      DO i = 1 , ngrainswl                                        ! initialize grains
        g(i) = 0.0d+00
      END DO
      p = p0
      DO i = 1 , iters
        WRITE(*,9020) i , iters, f, p
        DO l = 1 , ngrainswl                                      ! initialize Histogram
          H(l) = 0.0d+00
        END DO

        Nrej = 0
        Nacc = 0
        Nout = 0
        
        nnt = MaxTrials
        DO j = 1 , nnt

          DO k = 1 , ns                                          ! select new quantum numbers
            test = RAN1(idum)
            IF ( test .LE. p  ) THEN                             ! down-step
               nv(k) = nvold(k) - 1
               IF ( nv(k) .LT. 0 ) nv(k) = 0 
             ELSEIF (test .GT. p .AND. test .LE. 2.0d0*p )THEN   ! up-step

                nv(k) = nvold(k) + 1

             ELSE                                                ! sit-tight
               nv(k) = nvold(k)
            ENDIF
          END DO                                                 ! end quantum number selection

      IF((NY.EQ.0).AND.(NZ.EQ.0)) THEN
          enrg = energy( ns, wa, xa, nv)
      ELSE
          enrg = energyxyz( ns , wa , xa , ya, za, nv )             ! new energy above zpe
      ENDIF

      IF ( enrg.GE.0.9*Eswitch .AND. enrg.LE.Etop ) THEN      ! check energy boundary

      ntest = 0
      IF((NY.EQ.0).AND.(NZ.EQ.0)) THEN
         CALL ckderiv( ns, nv, wa, xa, ntest )            ! check derivatives for X-matrix 
      ELSE
         CALL ckderivxyz( ns, nv, wa, xa, ya, za, ntest )            ! check derivatives for XYZ-matrix 
      ENDIF

      ELSE
            ntest = 1                        ! fail test, ntest = 1
      ENDIF

      IF ( ntest.EQ.0 ) THEN            ! Accept or reject? pass: ntest = 0; fail: ntest = 1 
            nnew = IDNINT( (enrg - 0.9*Eswitch)/Egrain1 ) + 1    ! energy grain number: nearest INTEGER(4)
            acc = exp( g(nold)-g(nnew) )                         ! acceptance probability
            test = RAN1(idum)
            If ( test .LE. acc ) THEN                            ! accepted for a move
              Nacc = Nacc + 1        
              H(nnew) = H(nnew) + 1.d+00
              g(nnew) = g(nnew) + f
              nold = nnew                                        ! energy grain
                DO l = 1 , ns
                  nvold(l) = nv(l)                               ! quantum numbers
                END DO
             ELSE                                                ! rejected: sits tight
              Nrej = Nrej + 1
              H(nold) = H(nold) + 1.d+00
              g(nold) = g(nold) + f
            END IF
      ELSE                                                 ! rejected: tried to step out of energy range
              Nout = Nout + 1
              H(nold) = H(nold) + 1.d+00
              g(nold) = g(nold) + f
      ENDIF             
          
        END DO           ! j -> ntrials
        
        f = f/2.0d+00                      ! prepare for next iteration

      END DO             ! i -> iters

      WRITE(4,*) '   '
      WRITE(4,*) 'Finished Wang-Landau '
      WRITE(*,*) '***Finished Wang-Landau '
      CALL DateTime(4)                         ! print output to unit=4

      hmax = 0.0
      hmin = 1.e+10
      have = 0.0d+00
      nh = 0
      DO ig = 1 , ngrainswl-9    
        hmax = MAX( hmax, H(ig) )      ! Maximum histogram
        IF (H(ig) .GT. 1.e-2) THEN
          hmin = MIN( hmin, H(ig) )    ! minimum non-zero histogram
          have = have + H(ig)          ! sum of non-zero histogram values
          hav2 = hav2 + H(ig)**2       ! sum of squares
          nh = nh + 1
        END IF
      END DO
      have = have/nh                   ! average non-zero histogram
      hav2 = hav2/nh                   ! average non-zero sum of squares
      hmax = (hmax-have)*1.0d+02/have     ! % positive deviation
      hmin = (hmin-have)*1.0d+02/have     ! % negative deviation
      var = hav2 - have**2             ! variance

      hav2 = 0.0d+00
      DO ig = 1 , ngrainswl-9    
        IF (H(ig) .GT. 1.e-2) hav2 = hav2 + abs( H(ig) - have )       ! sum of un-signed errors
      END DO
      hav2 = hav2/nh                  ! average un-signed deviation from mean
      hav2 = hav2*1.0d+02/have        ! % average un-signed deviation from mean

c    PRINT OUTPUT

      DO iunit = 3, 5
        WRITE(iunit,9008) ttmax, KEYWORD, iters, p0 , fo
        WRITE(iunit,9016) 2100.*Nacc/nnt , 100.*Nrej/nnt , 
     &       100.*Nout/nnt ,
     &       have , SQRT(var) , hav2 , hmax , hmin
      END DO       
c --------------------------------------------------------------------------
c
c     Combining recursive direct count and bdens (first part)
c
c     Preliminary normalization
      gmax = -2000.d+00
      DO I = INT(0.9*Eswitch/Egrain1)+2, INT(Eswitch/Egrain1)+2
	kk = I-INT(0.9*Eswitch/Egrain1+1)
	IF ( g(kk) .GT. gmax ) gmax = g(kk)                           ! max g(kk)
      ENDDO

c     Sum up T of bdens and recursive direct count
      SUMAD = 0.0d0
      SUMDO = 0.0d0
      DO I = INT(0.9*Eswitch/Egrain1)+2, INT(Eswitch/Egrain1)+2
	SUMAD = SUMAD + T(I)                        ! T(I) from recursive direct count
	kk = I-INT(0.9*Eswitch/Egrain1+1)
	SUMDO = SUMDO + exp( g(kk) - gmax )                !  Wang-Landau with preliminary normalization
          enrg = (I-1)*Egrain1
      ENDDO

c     Calculate the normalization ratio
      IF (SUMAD .GT. 0.0d0 .AND. SUMDO .GT. 0.0d0) THEN
	CDOAD = SUMAD/SUMDO
      ELSE
	write(*,*)
	write(*,*)'*** FATAL: Wang-Landau Normalization failed ***'
	write(*,*)'    Too many empty grains: Try increasing Eswitch'
	write(*,*)
	write(3,*)
	write(3,*)'*** FATAL: Wang-Landau Normalization failed ***'
	write(3,*)'    Too many empty grains: Try increasing Eswitch'
	write(3,*)
	write(4,*)
	write(4,*)'*** FATAL: Wang-Landau Normalization failed ***'
	write(4,*)'    Too many empty grains: Try increasing Eswitch'
	write(4,*)
	STOP
      ENDIF
      
      DO I=INT(Eswitch/Egrain1 + 2), INT(Emax2/Egrain1 + 1)             ! normalize Wang-Landau from Eswitch to Emax2
	kk = I - INT(0.9*Eswitch/Egrain1 + 1)
     	T(I)  = exp( g(kk) - gmax )*CDOAD              ! final normalized Wang-Landau above Eswitch
	AT(I) = 0.0d0
      ENDDO
      
      ENDIF    !ENDIF for Wang-Landau section
      WRITE(*,*) '***End W-L Normalization '

c --------------------------------------------------------------------------
c   CREATE CHECKPOINT FILE
c --------------------------------------------------------------------------

      OPEN (50,STATUS='REPLACE',FILE=fname(1:lenstr(fname))//'.chk',
     &      ACTION='WRITE')  ! Chekpoint file
      
      WRITE(50,9911) AVERSION , ADATE
      WRITE(50,9912) AVERSION , ADATE 

      WRITE(50,9002) TITLE1
      WRITE(50,9002) TITLE2 
      CALL DateTime(50)                         ! print output to unit=50
      WRITE (50, 99052) LTMODE, Eswitch
      WRITE (50,9005)
      WRITE (50,99030) cut

      WRITE (50,99023) fname
      WRITE (50,99991) Egrain1 , imax1 , Emax2 , Isize , LTMODE ,Eswitch
      WRITE (50,99002)
      ngrains2 = INT(Emax2/Egrain1)+1           
      DO I = 1 , ngrains2
         enrg = (I-1)*Egrain1
         WRITE (50,*) I , enrg , T(I)
      END DO

      CLOSE (UNIT=50)                                                 ! Close chekpoint file

c --------------------------------------------------------------------------
c --------------------------------------------------------------------------
c --------------------------------------------------------------------------

      CASE (.TRUE.)                               ! START CONVOLUTION FROM chekpoint FILE
 
      OPEN (50,STATUS='OLD',FILE=fname(1:lenstr(fname))//'.chk',
     &     ACTION='READ',IOSTAT=istat)  ! Chekpoint file

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
         READ(50,99030) dumdum
      END DO

      READ (50,*) fnameT
      READ (50,*) Egrain1T , imax1T , Emax2T , IsizeT, LTMODET , MVALT
      
      itest = 0
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
       WRITE(iunit,*)' CALCULATED USING chekpoint FILE'
       WRITE(iunit,*) 
      END DO
      
      READ(50,*) dumdum
      READ(50,*) dumdum

      ngrains2 = INT(Emax2/Egrain1)+1           
      DO ig = 1 , ngrains2
        READ(50,*)  I , enrg , T(I)
      END DO
      
      CLOSE (UNIT=50)                                                 ! Close chekpoint file
      
      END SELECT
c --------------------------------------------------------------------------
c
c     Convolution of separable modes with non-separable via
c       the Beyer-Swinehart/Stein-Rabinovitch Algorithm
c
      IF ( Nsep .GT. 0 ) THEN       ! If separable modes are to be included
            
      DO I=ns+1, ns + Nsep
          IF( (IDOF(I).EQ.'HRA').OR.(IDOF(I).EQ.'HRB').OR.
     &         (IDOF(I).EQ.'HRC') ) THEN
              CALL SHRLEV( T,Egrain1,ngrains2,HRBb(I),HRVo(I), NG(I), 1,
     &        IMAX, HRzp(I) , 'VHRD1', 'BHRD1', 0.0d0, 0.0d0 )  ! rigid, symmetrical hindered rotor
              zpetot = zpetot + HRzp(I)
          ELSEIF (IDOF(I).EQ.'HRD') THEN
              DO J=1, NVV(I)
                 CVt(J)=CV(I,J)
              ENDDO
              DO J=1, NBB(I)
                 CBt(J)=CB(I,J)
              ENDDO
              CALL UHRLEV( T,Egrain1,ngrains2,NBB(I),NVV(I),CBt,CVt,
     &            NGV(I),NGB(I),IMAX,HRzp(I),Vhr(I),Bhr(I),
     &            Phav(I),Phab(I),NG(I) )   ! non-rigid, unsymmetrical hindered rotor
                 zpetot = zpetot + HRzp(I)
         ELSEIF (IDOF(I).EQ.'ROT') THEN
            Br=16.85763d0/RI(I)
            NSYMM=NRSN(I)
            NDIM=IDIM(I)
            CALL CROTLEV(T,Egrain1,ngrains2,Br,NDIM,NSYMM)
         ELSEIF (IDOF(I).EQ.'QRO') THEN
            Br=16.85763d0/RI(I)
            NSYMM=NRSN(I)
            NDIM=IDIM(I)
            CALL ROTLEV(T,Egrain1,ngrains2,Br,NDIM,NSYMM)
         ELSEIF ( IDOF(I).EQ.'TOP' ) THEN                         ! ********  Symmetric Top   *********
            NSYMM = NG(I)              ! Symmetry number
            CALL STOPLEV(T,Egrain1,ngrains2,B2,B1,NSYMM)
         ELSEIF (IDOF(I).EQ.'KRO') THEN
            write(*,*) '****FATAL: KRO type of d.o.f. not implemented'
            write(3,*) '****FATAL: KRO type of d.o.f. not implemented'
         ELSEIF (IDOF(I).EQ.'BOX') THEN
            Br=RI(I)
            NDIM=IDIM(I)
            CALL BOXLEV(T,Egrain1,ngrains2,Br,NDIM,boxzpe)
            zpetot = zpetot + boxzpe
         ELSEIF (IDOF(I).EQ.'VIB') THEN
            WEh=RI(I)
            XEh=RJ(I)
            CALL MORLEV(T,Egrain1,ngrains2,WEh,XEh,NG(I),morzpe)
            zpetot = zpetot + morzpe
         ELSE
            write(*,*) '************************************'
            write(*,*) 'FATAL: MISTAKES at convolution step'
            write(*,*) 'Degree of freedom #',I
            write(*,*) '************************************'
            write(4,*) '************************************'
            write(4,*) 'FATAL: MISTAKES at convolution step'
            write(4,*) 'Degree of freedom #',I
            write(4,*) '************************************'
            STOP
         ENDIF
      ENDDO
      ENDIF  ! separable modes

      DO iunit = 3, 5 
         WRITE (iunit,99013) zpetot
      END DO
      
c --------------------------------------------------------------------------
c
c     Compute sums and densities of states
c
      DS(1)=T(1)/Egrain1
      SS(1)=T(1)
      DO I= 2, ngrains2
         SS(I)=SS(I-1) + T(I)
         DS(I)=T(I)/Egrain1
      ENDDO

c --------------------------------------------------------------------------
c     Write bdens.out output file
c --------------------------------------------------------------------------

      WRITE(4,9017) 
      DO I=1, ngrains2-1
c       enrg = Eswitch + (I-1)*Egrain1
       enrg = (I-1)*Egrain1
       WRITE(4,9010) I,enrg,DS(I),SS(I)                  !,r_den(I),r_sum(I)
      ENDDO 
c --------------------------------------------------------------------------
c      Write to ___.dens file (UNIT=3) 
c 
      WRITE (3,9035) Eswitch, LTMODE
      WRITE (3,99030) cut                    ! end of data summary block in ____.dens file
      WRITE (5,99030) cut                    ! end of data summary block in ____.qvb file

      WRITE (3,*) fname(1:lenstr(fname))
      WRITE (3,*) TITLE1
      WRITE (3,99001) Egrain1 , imax1 , Emax2 , Isize , Viblo , zpetot

      WRITE (3,99010)
      DO 200 I = 1 , imax1
         enrg = (I-1)*Egrain1
         WRITE (3,99014) I , enrg , DS(I), SS(I)
 200  CONTINUE
 
      DO 300 I = imax1 + 1 , Isize
         enrg = (I-imax1-1)*Egrain2
         j = INT(enrg/Egrain1) + 1
         WRITE (3,99014) I , enrg , DS(j), SS(j)
 300  CONTINUE

      WRITE (5,*) fname(1:lenstr(fname))
      WRITE (5,*) TITLE1
      WRITE (5,99901) Egrain1 , Emax2 , zpetot, KEYWORD
99901 FORMAT (3(3x,f10.2),3x,A6)

      Ntx = 150
      WRITE(5,*) Ntx   ! *********** Number of temperatures to be passed to THERMO
      write(5,*)
      write(5,*) "==================================================="
      WRITE(5,*) "     VIBRATIONAL PARTITION FUNCTION"
      write(5,*) "==================================================="
      write(5,*)
      write(5,*) "   INDEX     T(K)       Q       C/R      H/R      S/R"

        TT=1.0d0
        DO I = 1 , Ntx   ! Number of Temperatures to be passed to THERMO
           Q0 = 0.0d0
           Q1 = 0.0d0
           Q2 = 0.0d0
           DO J=1, ngrains2-1 
             enrg = (J-1)*Egrain1
             Ered = 1.4387752D+00*enrg/TT
             fxx = DS(J)*exp(-Ered)*Egrain1
             Q0 = Q0 + fxx
             Q1 = Q1 + Ered*fxx
             Q2 = Q2 + Ered*Ered*fxx
           ENDDO   
c            Thermo formulae from Pitzer and Gwinn, J. Chem. Phys. 10, 428 (1942)
           Cx = Q2/Q0 - (Q1/Q0)**2
           Hx = Q1/Q0
           Sx = (Hx+log(Q0))
           write(5,9011) I, TT, Q0 , Cx , Hx , Sx
           TT = TT*1.05631305d+00  ! *********** for Ntx = 150
        ENDDO   

      write(5,*)

c
c       End of writing
c
      CLOSE (UNIT=3)
      CLOSE (UNIT=4)
c --------------------------------------------------------------------------

9911  FORMAT (
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'//,
     &'                     bdens-',A7,' (',A8,')',//,
     &'Copyright (C) 2016, 2017 Collin G.L. Li, Thanh Lam Nguyen, and'
     &' John R. Barker',//    
     &'        CONTACT:         John R. Barker',/
     &'                      (jrbarker@umich.edu)',/
     &'                     University of Michigan',/
     &'                  Ann Arbor, Michigan 48109-2143',//
     &'         http://clasp-research.engin.umich.edu/multiwell/',//
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%',/)
9912  FORMAT('Suggested Literature Citations:'//4x,
     &'a) J.R. Barker, T.L. Nguyen, J.F. Stanton, C. Aieta, M. Ceotto,',
     &/7x,'F. Gabas, T.J.D. Kumar, C.G.L.Li, L.L.Lohr, A. Maranzana,',
     &/7x,'N.F. Ortiz, J.M. Preses, J.M. Simmie, J.A. Sonk, and ',
     &'P.J. Stimac; ',
     &/7x,'MultiWell-',A7,' Software Suite; University of Michigan, ',
     &/7x,'Ann Arbor, Michigan, USA, ',A8,'; ',
     &/7x,'http://clasp-research.engin.umich.edu/multiwell/.',
     &//4x,'b) J.R. Barker, Int. J. Chem. Kinetics, 33, 232-45 (2001)',
     &//4x,'c) J.R. Barker, Int. J. Chem. Kinetics, 41, 748-763 (2009)',
     & //4x,
     &'d) F. Wang and D.P. Landau, Phys. Rev. Letters 86, ',/7x,
     &'2050-2053 (2001).',//4x,
     &'e) M. Basire, P. Parneix, and F. Calvo, J. Chem. Phys. 129,',/7x,
     &'081101 (2008).',//4x,    
     &'f) T.L. Nguyen and J.R. Barker, J. Phys. Chem. A., 114,',/7x,
     &'3718-3730 (2010).',
     &//,
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'/)
9002  FORMAT(A100)
9003  FORMAT(I3,' Vibrational Frequencies (cm-1)  ',/,
     & '  We; E = We*(v+1/2) +...; Wo: E = W0*v +...;  0<->1 fund.  ',/,
     & '  De: BDE for Separable Morse oscillator [0.0 if non-Morse]',//,
     & '  No.     We         Wo      0-1_Fund.   De(cm-1)' )
9004  FORMAT( I5,3(1x, F10.3),1x,F10.1 )
9042  FORMAT( I3, 100( 1x , f9.3 ) )
9043  FORMAT( 100(7x,I3) )
9005  FORMAT(/)
9006  FORMAT('Anharmonicity Matrix (cm-1)')
9008  FORMAT(/'Monte Carlo trials: ', 1pe8.1,' (',A6,')',/,
     &        '  iterations = ',I2,/,
     &        '          p0 = ',0pf7.3,/,
     &        '          f  = ',0pf7.3,/)
9010  FORMAT( 2x,I6,2x, F10.3,5(2x,1pe10.4) )
9011  FORMAT( 2x,I6,2x, F13.5,5(2x,1pe12.5) )
9013  FORMAT(/,' Ave:', 3(1x, F10.3)/
     &         '     Partial ZPEanh (cm-1) = ', f10.3,
     &         ' (omits Go from VPT2 and separable D.O.F.)',/)
9016  FORMAT( '             steps up or down:',f10.3,' %  ',/,
     &        '               steps in place:',f10.3,' %  ',/,
     &        '           steps out of range:',f10.3,' %  ',/,
     &        '   average non-zero histogram: ',f10.0,' +/-', f8.0,/,
     &        '  average un-signed deviation:',0pf10.2,' %  ',/,
     &        '       max positive deviation:',0pf10.2,' %  ',/,
     &        '        max negative devation:',0pf10.2,' %',//)
9017  FORMAT(
     & '     No.        cm-1  states/cm-1  SumStates')
9020  FORMAT(' Iteration =',I3,'/',I2,5x,'f = ',1pe10.2, 5x,
     &                'p = ',1pe10.4)
9030  FORMAT(I3,2X,'Hind.Rot',2x,'Freq(har in cm-1)=',F8.2,2x,
     &     'Mom(amu.A**2)=',F8.3,2x,'fold=',I2,2x,'symm=',I2,2x,
     &     'B=',F8.4,' cm-1',2x,'Uo=',F7.1,' cm-1')
9031  FORMAT (I3,2X,'General Hindered Rotor:',2x,'Rotor symm. =',I2)
9032  FORMAT (8X,A5,1x,' ; symmV =',I2,2x,'; Phase(rad.) = ',
     & F8.4,2x,'; Coeff. =',20(F10.4,1x))

9033  FORMAT (8X,A5,1x,' ; symmB =',I2,2x,'; Phase(rad.) = ',
     & F8.4,2x,'; Coeff. =',20(F10.4,1x))
9035  FORMAT('Eswitch =', f10.1, ' cm-1 ('A4')'/)
99001 FORMAT (1x,f10.1,1x,I10,1x,f10.1,1x,I10,1x,f10.1,1x,f10.2)
99002 FORMAT (/,'T(I) array for coupled degrees of freedom:',/,
     &        '        I      (cm-1)   T(I)')
99010 FORMAT ('       No.     (cm-1)   Density    Sum')
99013 FORMAT (/,'TOTAL ZPEanh: ', f10.3,' cm-1 for all coupled & ',
     &  'separable D.O.F.; omits Go from VPT2.',/)
99014 FORMAT (i10,1x,f10.1,2(1x,1pe10.4))
99023 FORMAT (A10)
99030 FORMAT (A46)
99041 FORMAT (I3,2X,'Quantized Rotor:  I(Amu.A**2) = ',
     &     F8.3,' ;  Sym. = ',I2,' ;  Dim. = ',I2,' ;  B(cm-1) = ',F8.3)
99042 FORMAT (I3,2X,'K-Rotor:  I(Amu.A**2) = ',F8.3,' ;  Sym. = ',I2,
     &          ' ;  Quantum No. J = ',I2,' ;  B(cm-1) = ',F7.3)
99043 FORMAT (I3,2X,'Particle-in-Box:  Freq(cm-1)  = ',F7.2,
     &             ' ;  Dim. = ',I2)
99044 FORMAT (I3,2X,'Classical Rotor:  I(Amu.A**2) = ',
     &     F8.3,' ;  Sym. = ',I2,' ;  Dim. = ',I2,' ;  B(cm-1) = ',
     &     F8.3)
99045 FORMAT (I3,2X,'Vibration:',8x,'Freq(cm-1)  = ',F7.2,
     &' ;  Anharm. = ', F8.3)
99050 FORMAT (//,' SEPARABLE DEGREES OF FREEDOM',/)
99051 FORMAT (/' --- END INPUT ---')
99052 FORMAT ( A6, ': Eswitch = ', F7.1 )
99905 FORMAT (I5,2X,'Symm-Top',2x,'2D-Moment =',F8.2,' amua, B2 =',F8.4,
     &       ' cm-1  ;  1D-Moment =',F8.2,' amua, B1 =',F8.4,' cm-1',2x,
     &       ';  Symm. No. = ',I1)
99991 FORMAT (1x,f10.1,1x,I10,1x,f10.1,1x,I10,1x,A10,1x,f10.2)

      END

