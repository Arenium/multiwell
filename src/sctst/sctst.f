c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Copyright (C) 2016 John R. Barker and Thanh Lam Nguyen
c
c Authors: John R. Barker and Thanh Lam Nguyen
c      jrbarker@umich.edu     LNGUYEN@cm.utexas.edu
c
c Contact:
c John R. Barker
c jrbarker@umich.edu
c University of Michigan
c Ann Arbor, MI 48109-2143
c (734) 763 6239
c
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License (version 2)
c as published by the Free Software Foundation.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
c GNU General Public License for more details.
c
c See the 'ReadMe' file for a copy of the GNU General Public License,
c or contact:
c
c Free Software Foundation, Inc.
c 59 Temple Place - Suite 330
c Boston, MA 02111-1307, USA.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c                                                                 
c                         PROGRAM sctst    
c
c                               by      
c
c          John R. Barker       and      Thanh Lam Nguyen                       
c          jrbarker@umich.edu            nguyenlt@umich.edu
c                                                                  
c            ***Semi-Classical Transition State Theory***
c
c                             based on
c
c          The theory of W. H. Miller and coworkers* and the 
c           Wang-Landau algorithms for densities of states**
c
c
c   Literature Citations (basis for this program):
c
c   T. L. Nguyen, J. F. Stanton, J. R. Barker, Chem. Phys. Letters, 
c              499, 9-15 (2010).
c
c    *Semi-Classical Transition State Theory
c    W. H. Miller, J. Chem. Phys. 62, 1899-1906 (1975).
c    W. H. Miller, Faraday Discuss. Chem. Soc. 62, 40-46 (1977).
c    W. H. Miller, R. Hernandez, N. C. Handy, D. Jayatilaka, and A. Willets,
c      Chem. Phys. Letters 172, 62-68 (1990).
c    R. Hernandez and W. H. Miller, Chem. Phys. Lett. 214, 129-136 (1993).
c    J. F. Stanton, J. Phys. Chem. Lett. 7, 2708-2713 (2016)
c                                                                  
c    **Density of states algorithms
c    F. Wang and D. P. Landau, Phys. Rev. Letters 86, 2050-2053 (2001).           
c    M. Basire, P. Parneix, and F. Calvo, J. Chem.Phys. 129, 081101 (2008).    
c    T. L. Nguyen and J.R. Barker, J. Phys. Chem. A., 114, 3718–3730 (2010).
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c Copyright (C) 2016 John R. Barker and Thanh Lam Nguyen
c
c Authors: John R. Barker and Thanh Lam Nguyen
c      jrbarker@umich.edu     nguyenlt@umich.edu
c          September 2009
c
c John R. Barker
c jrbarker@umich.edu
c University of Michigan
c Ann Arbor, MI 48109-2143
c (734) 763 6239
c
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License (version 2)
c as published by the Free Software Foundation.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
c GNU General Public License for more details.
c
c See the 'ReadMe' file for a copy of the GNU General Public License,
c or contact:
c
c Free Software Foundation, Inc.
c 59 Temple Place - Suite 330
c Boston, MA 02111-1307, USA.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      PROGRAM sctst

      IMPLICIT NONE

      CHARACTER(len=7) AVERSION 
      CHARACTER(len=8) ADATE 
      PARAMETER ( AVERSION='2023', ADATE='Mar 2023' )
      CHARACTER(46) cut , dumdum
      PARAMETER (cut='**************INPUT DATA SUMMARY**************')
      
      CHARACTER(10) fname , fnameT , chekpoint
      
      INTEGER istat , itest , iunit

      INTEGER(4) ns, nv(100), i, j, k, MaxTrials,
     & idum, l , nstart , Nacc, Nrej, Nout,
     & nvold(100), iters, ngrains , nold , nnew , 
     & ndmax , nh , ig , nnt, nvmax , KIN , Nt
      PARAMETER(KIN=2)
      
      REAL(8) wa , w0 , wf , xa , ya , za ,
     &                 Eu, trialmax, ttmax,
     &                 Emin, Emax, Estart, energy, 
     &                 RAN1, D(100),
     &                 enrg, zpe , boxzpe , morzpe , zpetot ,
     &                 ave , av0, avf , Etop
      REAL(8) fx , Ered, Q0, Q1, Q2, Cx, Hx, Sx
     
      REAL(8) H(20002), g(20002), p , p0 , f , fo , test , 
     &                 acc , hmax , hmin , have, hav2, y , var
 
      CHARACTER(100) TITLE1 , TITLE2
      CHARACTER(2) WW
      CHARACTER(6) KEYWORD1 , KEYWORD2 , Eunits

      CHARACTER(3) IDOF(100)
      CHARACTER(5) Vhr(100) , Bhr(100)
      CHARACTER(4) VROTIN , VROTOUT
      
      REAL(8) HRfreq(100), HRVo(100)  , HRI(100) , 
     &        HRB(100)   , HRzpe(100) , DUM1     , DUM2
      REAL(8) CV(100,100), CB(100,100), CVt(100) , CBt(100) 
      REAL(8) Phav(100)  , Phab(100)
      INTEGER(4) Nsep, NG(100), IMAX(100), NVV(100), NBB(100), 
     &     NGV(100), NGB(100)
      REAL(8) T(20002), AT(20002), DS(20002), SS(20002)
      
      REAL(8) WE , XE
      
      REAL(8) BB , BG , FAC , RG , GAM
      INTEGER N , NR

      INTEGER(4) NY, NZ , MODE(100)
      INTEGER(4) nvmaxxyz, JMAX 
      REAL(8) energyxyz
     
      REAL(8) B
      INTEGER(4) lenstr

      REAL(8) RI(100) , RJ(100)
      INTEGER(4) NRSN(100), IDIM(100)
      INTEGER(4) NDIM, NSYMM, JK
 
      REAL(8) Egrain1 , Egrain2 , Emax1 , Emax2 , Viblo
      INTEGER(4) imax1 , imax2 , Isize , NMAX
      PARAMETER (NMAX=50001)  

      REAL(8) Egrain1T , Emax2T 
      INTEGER(4) imax1T , IsizeT

      INTEGER(4) ntest
      
      REAL(8) Vo , FI , xFF, XI(100) , PP
      CHARACTER VPTx*5
      REAL(8) PE, eng, DE, omega, TTF(NMAX), Ev(NMAX) , gamma
      REAL(8) temp , PN
      INTEGER(4) Ntrial , Ntx

      REAL(8) B1 , B2 , X1 , X2
      INTEGER(4) ii , ntop

      REAL(8) Vfi, Vri, Vf, Vr, DelH, TT, RR, beta0, su
      PARAMETER (RR=0.695038916D0)
      
      CHARACTER chekname*15

      COMMON/I/ idum
      COMMON/wxyz/ wa(100) , w0(100) , wf(100) , xa(100,100) ,
     &   ya(100,100,100) , za(100,100,100,100) , zpe , ave ,
     &   av0 , avf , Viblo

      EXTERNAL energy , energyxyz, RAN1 , ndmax, nvmax, nvmaxxyz
      
      SAVE
      
      OPEN (2,FILE='sctst.dat',STATUS='OLD')
      READ (2,*) fname

      OPEN (3,STATUS='UNKNOWN',FILE=fname(1:lenstr(fname))//'.crp')  ! Succinct output for MultiWell input
      OPEN (4,FILE='sctst.out',STATUS='UNKNOWN')
      OPEN (5,STATUS='UNKNOWN',FILE=fname(1:lenstr(fname))//'.qcrp')  ! Succinct output for Thermo input

      READ (2,9002)  TITLE1
      READ (2,9002)  TITLE2

      WRITE(3,99030) cut                                  ! Start of data summary block in ____.crp file
      WRITE(5,99030) cut                                  ! Start of data summary block in ____.crp file

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

      CALL READ_WXYZ(KIN , ns , NY, NZ, WW)                ! read frequencies and anharmonicities
      
      READ(2,*) NSEP , VROTIN                              ! No. of separable DOF and keyword for rotation units
      CALL ucase ( VROTIN )   ! convert to upper case
      VROTOUT = 'AMUA'

      DO I=ns+1, ns+NSEP
         READ(2,*) MODE(I), IDOF(I), DUM1, DUM2, NG(I)
         CALL ucase ( IDOF(I) )   ! convert to upper case
         IF (IDOF(I).EQ.'HRA') THEN
            HRfreq(I)=DUM1
            HRI(I)=DUM2
            HRB(I)=16.85763d0/HRI(I)
            HRVo(I)=((HRfreq(I)/NG(I))**2)/HRB(I)
         ELSEIF (IDOF(I).EQ.'HRB') THEN
            HRfreq(I)=DUM1
            HRVo(I)=DUM2
            HRB(I)=((HRfreq(I)/NG(I))**2)/HRVo(I)
            HRI(I)=16.85763d0/HRB(I)
         ELSEIF (IDOF(I).EQ.'HRC') THEN
            HRI(I)=DUM1
            HRVo(I)=DUM2
            HRB(I)=16.85763d0/HRI(I)
            HRfreq(I)=NG(I)*sqrt(HRVo(I)*HRB(I))
         ELSEIF (IDOF(I).EQ.'HRD') THEN
            NVV(I)=INT(DUM1)
            NBB(I)=INT(DUM2)
            READ(2,*) Vhr(I), NGV(I), Phav(I), (CV(I,J), J=1, NVV(I))
            READ(2,*) Bhr(I), NGB(I), Phab(I), (CB(I,J), J=1, NBB(I))
            CALL ucase ( Vhr(I) )   ! convert to upper case
            CALL ucase ( Bhr(I) )   ! convert to upper case
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
            NRSN(I) = DUM2
            IDIM(I) = NG(I)
            CALL rotunits( VROTIN , RI(I) , VROTOUT)
         ELSEIF (IDOF(I).EQ.'QRO') THEN
            RI(I) = DUM1
            NRSN(I) = DUM2
            IDIM(I) = NG(I)
            CALL rotunits( VROTIN , RI(I) , VROTOUT)
         ELSEIF ( IDOF(I).EQ.'TOP' ) THEN                               ! Symmetric Top with 2D and 1D rot constants B2 and B1
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
            CALL rotunits( VROTIN , RI(I) , VROTOUT)
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
 
      READ (2,*) Egrain1 , imax1 , Isize , Emax2 , KEYWORD1
         CALL ucase ( KEYWORD1 )   ! convert to upper case
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

      READ(2,*) chekpoint , chekname          ! Start calculation from chekpoint file if chekpoint='chekstart' or 'checkstart'; chekname is checkpoint filename
          CALL ucase ( chekpoint )   ! convert to upper case
          IF ( chekpoint .EQ. 'CHECKSTART' ) chekpoint = 'CHEKSTART' 
          
      READ(2,*) KEYWORD2           ! Monte Carlo trials/bin: 10    10^2  10^3  10^4    10^5  10^6
                                   !              KEYWORD2 = poor, fair, good, better, best, extra
          CALL ucase ( KEYWORD2 )   ! convert to upper case
      READ(2,*) VPTx                 ! key word for VPT4 and/or Wagner correction
          CALL ucase ( VPTx )   ! convert to upper case
c         VPTx  = keyword ( CHARACTER*5 )
c                 VPT2: VPT2, only
c                 VPT4A: VPT2 + VPT4 correction using Eq. 42 from Stanton, J. Phys. Chem. Lett. 7, 2708-2713 (2016)
c                 VPT4B: VPT2 + VPT4 correction using Eq. 37 from Stanton, J. Phys. Chem. Lett. 7, 2708-2713 (2016)
c                 VPT4W: as in VPT4A + Wagner semiempirical correction

      READ(2,*) Vfi, Vri, Eunits     ! Barrier heights (including zpe) in forward and reverse directions, Eunits='kcal', 'kj', or 'cm-1'
          CALL ucase ( Eunits )   ! convert to upper case      
      IF ( Eunits .EQ. 'KCAL' ) THEN
        Vf = Vfi*349.754D0
        Vr = Vri*349.754D0
      ELSEIF ( Eunits .EQ. 'KJOU' ) THEN
        Vf = Vfi*83.5932D0
        Vr = Vri*83.5932D0
      ELSEIF ( Eunits .EQ. 'CM' ) THEN
        Vf = Vfi
        Vr = Vri
      ELSE
        WRITE(*,*) '***FATAL: Energy units (Eunits) not recognized***'
        WRITE(4,*) '***FATAL: Energy units (Eunits) not recognized***'
        STOP
      ENDIF

      IF(Vf.GT.Vr) THEN                     ! to be used later to set energy zero (for Q) to reactant energy (units of cm-1)
        Vo = Vr                             ! Vo = barrier height in exothermic direction
        DelH = Vf-Vr
      ELSE
        Vo = Vf
        DelH = 0.0d0
      ENDIF

      READ(2,*) FI, xFF              ! imaginary frequency (cm-1) and diagonal anharmonicity (cm-1) for rxn coord.
      READ(2,*) ( XI(I), I=1,ns)     ! off-diagonal anharmonicities (cm-1) involving rxn coord.
      
      CLOSE (UNIT=2)
c
c   ******* END INPUT *******
c
      Emax = Emax2 + Egrain1                     ! Top of Wang-Landou energy window
      Emin = 0.0                                 ! Bottom of Wang-Landou energy window
      
      Emax1 = Egrain1*(imax1-1)
      Egrain2 = Emax2/(Isize-imax1-1)
      imax2 = Isize - imax1
 
      JMAX = 1 + INT(Emax2/Egrain1)
      JMAX = MIN(JMAX,NMAX)

      CALL convib( ns , WW )                      ! convert input frequencies to other forms

      DO iunit = 3, 5

      WRITE(iunit,9003) ns
      DO i = 1 , ns
        D(i) = 0.0
        IF (xa(i,i) .LT. 0.0) D(i) = -0.25D+00*w0(i)**2/xa(i,i) ! Bond dissociation energy
        WRITE(iunit,9004) i, wa(i), w0(i), wf(i), D(i)
      END DO

      WRITE(iunit, 9013) ave, av0, avf    ! average frequencies
      WRITE(iunit,9006)
      WRITE(iunit,9043) (i, i=1,ns)

      DO j = 1 , ns
        WRITE(iunit,9042) j, (xa(j,i), i=1,ns)
      END DO

      WRITE(iunit, 9015) zpe    ! ZPEanh
      write(iunit,99050) Nsep
      IF ( Nsep .EQ. 0 ) THEN
         write(iunit,*)  '       (none)'
      ENDIF

      END DO    ! END WRITE iunit BLOCK            
c
c    Write out separable degrees of freedom
c
      DO I=ns+1, ns+Nsep
            IF (IDOF(I).EQ.'HRA') THEN
             DO iunit = 3, 5
             write(iunit,9030) MODE(I), HRfreq(I), HRI(I), NG(I), NG(I), 
     &               HRB(I), HRVo(I)
             END DO
            ELSEIF (IDOF(I).EQ.'HRB') THEN
             DO iunit = 3, 5
             write(iunit,9030) MODE(I), HRfreq(I), HRI(I), NG(I), NG(I), 
     &               HRB(I), HRVo(I)
             END DO
            ELSEIF (IDOF(I).EQ.'HRC') THEN
             DO iunit = 3, 5
             write(iunit,9030) MODE(I), HRfreq(I), HRI(I), NG(I), NG(I), 
     &               HRB(I), HRVo(I)
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
                  B = 16.85763d0 / RI(I)
             DO iunit = 3, 5
                  write(iunit,99044) MODE(I), RI(I), NRSN(I), IDIM(I), B
             END DO
            ELSEIF (IDOF(I).EQ.'QRO') THEN
                  B = 16.85763d0 / RI(I)
             DO iunit = 3, 5
                  write(iunit,99041) MODE(I), RI(I), NRSN(I), IDIM(I), B
             END DO
            ELSEIF ( IDOF(I).EQ.'TOP' ) THEN                               ! Symmetric Top with 2D and 1D rot constants B2 and B1
               DO iunit = 3, 5
                 WRITE (iunit,99905) MODE(I) , X2 , B2 , X1 , B1 , NG(I)
               END DO
            ELSEIF (IDOF(I).EQ.'KRO') THEN
                  B = 16.85763d0 / RI(I)
             DO iunit = 3, 5
                  write(iunit,99042) MODE(I), RI(I), NRSN(I), IDIM(I), B
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

      DO iunit = 3, 5
        write(iunit,99051)
      END DO
c-----------------------------------------------------------------------------------------
c
c         WANG-LANDAU PARAMETERS
C
c-----------------------------------------------------------------------------------------
      p0 = 1.0/ns                                       ! Scaled according to number of vibrations
      IF ( p0 .GT. 0.25d+00 ) p0 = 0.25d+00
      iters = 21                                        ! Always
      ngrains = (Emax-Emin)/Egrain1 + 1                 ! should be same as JAMX (when Emin = 0.0)
      Etop = Emax + Egrain1
      fo = 1.0d+00

c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
        
      SELECT CASE ( chekpoint .EQ. 'CHEKSTART')   !! NOTE SPELLING
      
      CASE (.FALSE.)    ! do NOT start from saved chekpoint file, but create a new one
      
c-----------------------------------------------------------------------------------------
c   Trials per energy grain
c
      ttmax = trialmax( KEYWORD1 , ngrains )
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
       Estart = Emin + RAN1( idum )*(Emax-Emin) 
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
        nold = IDNINT( (enrg - Emin)/Egrain1 ) + 1              ! energy grain number: nearest INTEGER(4)
            
      DO i = 1 , ngrains                                        ! initialize grains
        g(i) = 0.0d+00
      END DO
      
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
c     Start Iterations
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
      p = p0
      DO i = 1 , iters
        WRITE(*,9020) i , iters, f, p
        DO l = 1 , ngrains                                      ! initialize Histogram
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

      IF ( enrg.GE.Emin .AND. enrg.LE.Emax ) THEN      ! check energy boundary

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
            nnew = IDNINT( (enrg - Emin)/Egrain1 ) + 1           ! energy grain number: nearest INTEGER(4)
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
      ELSE                                                 ! rejected: failed tests
              Nout = Nout + 1
              H(nold) = H(nold) + 1.d+00
              g(nold) = g(nold) + f
      ENDIF             
          
        END DO           ! j -> ntrials
        
        f = f/2.0d+00                      ! prepare for next iteration

      END DO             ! i -> iters

      CALL DateTime(4)                         ! print output to unit=4

      hmax = 0.0
      hmin = 1.e+10
      have = 0.0d+00
      nh = 0
      DO ig = 1 , ngrains    
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
      DO ig = 1 , ngrains    
        IF (H(ig) .GT. 1.e-2) hav2 = hav2 + abs( H(ig) - have )       ! sum of un-signed errors
      END DO
      hav2 = hav2/nh                  ! average un-signed deviation from mean
      hav2 = hav2*1.0d+02/have        ! % average un-signed deviation from mean

c    PRINT OUTPUT

      DO iunit = 3, 5
        WRITE(iunit,9008) ttmax, KEYWORD1, iters, p0 , fo
        WRITE(iunit,9016) 100.*Nacc/nnt , 100.*Nrej/nnt , 
     &       100.*Nout/nnt ,
     &       have , SQRT(var) , hav2 , hmax , hmin
      END DO       
       
c --------------------------------------------------------------------------

      DO I=1, ngrains
          T(I)=exp(g(I)-g(1))      ! Normalize and obtain sums of states
          AT(I)=0.0d0
      ENDDO

c --------------------------------------------------------------------------

      OPEN (50,STATUS='REPLACE',FILE=chekname, ACTION='WRITE')  ! Chekpoint file
      
      WRITE(50,9911) AVERSION , ADATE
      WRITE(50,9912) AVERSION , ADATE 
      
      WRITE(50,9002) TITLE1
      WRITE(50,9002) TITLE2 
      CALL DateTime(50)                         ! print output to unit=50
      WRITE(50,9005)

      WRITE (50,99030) cut
      WRITE (50,99023) fname
      WRITE (50,99001) Egrain1 , imax1 , Emax2 , Isize
      WRITE (50,99002)
      DO I = 1 , ngrains
         enrg = Emin + (I-1)*Egrain1
         WRITE (50,*) I , enrg , T(I)
      END DO

      CLOSE (UNIT=50)                                                 ! Close chekpoint file (finished writing)

c --------------------------------------------------------------------------
c --------------------------------------------------------------------------
c --------------------------------------------------------------------------

      CASE (.TRUE.)                               ! START FROM existing chekpoint FILE

      OPEN (50,STATUS='OLD',FILE=chekname, ACTION='READ',IOSTAT=istat)  ! Chekpoint file

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

      itest = 0
      READ (50,*) fnameT
      READ (50,*) Egrain1T , imax1T , Emax2T , IsizeT
      
c      IF ( fnameT .NE. fname)     itest = 1
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
       WRITE(iunit,*)' CALCULATED USING chekpoint FILE: ', chekname
       WRITE(iunit,*) 
      END DO
      
      READ(50,*) dumdum
      READ(50,*) dumdum
      
      DO ig = 1 , ngrains
        READ(50,*)  I , enrg , T(I)
      END DO
      
      CLOSE (UNIT=50)                                                 ! Close chekpoint file
      
      END SELECT
c --------------------------------------------------------------------------
c --------------------------------------------------------------------------

      WRITE(*,*)	
      WRITE(*,*) 'sctst has finished the orthogonal coupled DOF'
      WRITE(*,*) '  and is now running the semi-classical part.'
      WRITE(*,*)	
      
      ttmax = trialmax( KEYWORD2 , ngrains )

      DO iunit = 3 , 5
        WRITE(iunit,99060) ttmax, KEYWORD2, Vfi, Vri, Eunits, Vf, Vr
        IF ( VPTx .EQ. 'VPT2' )  WRITE(iunit,99063) VPTx
        IF ( VPTx .EQ. 'VPT4A' ) WRITE(iunit,99064) VPTx
        IF ( VPTx .EQ. 'VPT4B' ) WRITE(iunit,99065) VPTx
        WRITE(iunit,99061) FI, xFF
        WRITE(iunit,99062) (I , -XI(I) , XI(I) , I=1 , ns)
      END DO

        IF ( VPTx .EQ. 'VPT2' ) THEN                        ! VPT4 corrections applied to imaginary frequency
           gamma = 1.0d+00
        ELSEIF ( VPTx .EQ. 'VPT4A' ) THEN
           gamma = SQRT( 1.e+00 - (xFF/FI)**2 )
        ELSEIF ( VPTx .EQ. 'VPT4B' ) THEN
           gamma = 1.e+00 - 0.5d+00*(xFF/FI)**2
        ELSE
           WRITE(*,*) ' '
           WRITE(*,*) '******FATAL*******'
           WRITE(*,*) 'VPTx keyword not recognized or not implemented' 
           WRITE(*,*) ' '
           STOP
        ENDIF
        FI = FI*gamma
c
c The 22nd iter
c
        DO I=1, ngrains
                TTF(I) = 0.0d0
                Ev(I) = 0.0d0
        END DO

      Emin = 0.0d0
c
c     Initialize quantum numbers
c
      ntest = 1
      DO WHILE ( ntest .EQ. 1 )
       Estart = Emin + RAN1( idum )*(Emax-Emin) 
        DO l = 1, ns
          DO k=l+1, ns
               nvold(k)=0
          ENDDO
            Eu=Estart-energy( ns, wa, xa, nvold )
            IF(Eu .NE. 0.0d0) THEN
               nstart=nvmax( ns, wa, xa, nvold, l, Eu )
               nvold(l) = nstart*RAN1( idum )
            ELSE
               nvold(l) = 0
            ENDIF
        ENDDO
         CALL ckderiv( ns, nvold, wa, xa, ntest )		         ! check derivatives returns ntest=0,1 for PASS,FAIL
      END DO                                                    ! END of DO WHILE
      
          enrg = energy( ns, wa, xa, nvold )
        nold = IDNINT( (enrg - Emin)/Egrain1 ) + 1              ! energy grain number: nearest integer
      
!!!	write(*,*) 'OUT-initial selection'
!!!	write(*,*) enrg, nold
	
      
      DO i = 1 , ngrains                                        ! initialize grains
        g(i) = 0.0d+00
      END DO
      
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
c     Start Iterations
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
      p = p0
	f = 1.0d0
        DO l = 1 , ngrains                                      ! initialize Histogram
          H(l) = 0.0d+00
        END DO

        Nrej = 0
        Nacc = 0
        Nout = 0
  
	Emin = 0.0d0

	nnt = INT(1.0E+4*ngrains)

	IF(nnt.LT.MaxTrials) nnt = MaxTrials
	
          ttmax = trialmax( KEYWORD2 , ngrains )
          nnt = INT( ttmax )
	
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

          enrg = energy( ns, wa, xa, nv)

	IF ( enrg.GE.Emin .AND. enrg.LE.Emax ) THEN	! check energy boundary
	   ntest = 0
             CALL ckderiv( ns, nv, wa, xa, ntest )		          ! check derivatives returns ntest=0,1 for PASS,FAIL
	ELSE
	   ntest = 1				! fail test, ntest = 1
	ENDIF

	IF ( ntest .EQ. 0 ) THEN            ! Accepted       pass: ntest = 0; fail: ntest = 1 
            nnew = IDNINT( (enrg - Emin)/Egrain1 ) + 1           ! energy grain number: nearest integer
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
	ELSE                                                 ! rejected:  failed tests
              Nout = Nout + 1
              H(nold) = H(nold) + 1.d+00
              g(nold) = g(nold) + f
	ENDIF             

        enrg = energy( ns, wa, xa, nvold)

        temp = omega(ns, nvold, FI, XI)                          ! imaginary frequency, corrected for anharmonic coupling to orthogonal modes
        TTF(nold) = TTF(nold) + temp                             ! sum imaginary frequency, corrected for anharmonic coupling to orthogonal modes
        Ev(nold) = Ev(nold) + enrg
        
        END DO     ! j -> ntrials

        DO I=1, ngrains
          IF( H(I).GT.1.d-02 ) THEN
             TTF(I) = TTF(I) / H(I)                 ! average imaginary frequency, corrected for anharmonic coupling to orthogonal modes
             Ev(I) = Ev(I) / H(I)                   ! average energy of samples within the bin, measured from vib zpe
          ENDIF
        ENDDO

        DO I=1, ngrains                             ! Eq. 10 in Nguyen et al., Chem. Phys. Lett. 499, 9-15, 2010
           SS(I)=0.0d0
           eng = (I-1)*Egrain1
           DO J=1, I                                       ! contributions from all lower occupied grains and the current grain I
              IF( eng.GE.Ev(J) .AND. H(J).GT.1.d-02 ) THEN 
                 DE = eng - ( Vo + Ev(J) )                 ! new definition of DE (opposite sign)
                 CALL CalPN(PP, DE, TTF(J), xFF)           ! Semiclassical tunneling transmission probability
                 SS(I) = SS(I) + PP*T(J)                   ! T(J) is the number of orthogonal states in Jth grain; SS(I) is the "cumulative rxn probability"
              ENDIF
           ENDDO
        ENDDO

        T(1) = SS(1)
        AT(1)=0.0d0
        DO I=2, ngrains
          T(I)=SS(I)-SS(I-1)
          IF(T(I).LT.0.0d0) T(I)=0.0d0
          AT(I)=0.0d0
        ENDDO
        
C
C AT() and T() array is required
C
c --------------------------------------------------------------------------
c
c     Convolution of separable modes with non-separable via
c       the Beyer-Swinehart/Stein-Rabinovitch Algorithm
c
      zpetot = zpe
      IF ( Nsep .GT. 0 ) THEN       ! If separable modes are to be included
      
      DO I=ns+1, ns + Nsep
          IF((IDOF(I).EQ.'HRA').OR.(IDOF(I).EQ.'HRB').OR.
     &         (IDOF(I).EQ.'HRC')) THEN
              CALL SHRLEV(T,AT,Egrain1,ngrains,HRB(I),HRVo(I),NG(I),1,
     &        IMAX(I),HRzpe(I),
     &        'Vhrd1','Bhrd1',0.0d0,0.0d0,NG(I))                      ! Symmetrical hindered rotor
            zpetot = zpetot + HRzpe(I)
          ELSEIF (IDOF(I).EQ.'HRD') THEN
           DO J=1, NVV(I)
              CVt(J)=CV(I,J)
           ENDDO
           DO J=1, NBB(I)
              CBt(J)=CB(I,J)
           ENDDO
           CALL UHRLEV(T,AT,Egrain1,ngrains,NBB(I),NVV(I),CBt,CVt,
     &         NGV(I),NGB(I),IMAX(I),HRzpe(I),Vhr(I),Bhr(I),
     &         Phav(I),Phab(I),NG(I))                                 ! General, unsymmetrical hindered rotor
            zpetot = zpetot + HRzpe(I)
         ELSEIF (IDOF(I).EQ.'ROT') THEN
            B=16.85763d0/RI(I)
            NSYMM=NRSN(I)
            NDIM=IDIM(I)
            CALL CROTLEV(T,Egrain1,ngrains,B,NDIM,NSYMM)
         ELSEIF (IDOF(I).EQ.'QRO') THEN
            B=16.85763d0/RI(I)
            NSYMM=NRSN(I)
            NDIM=IDIM(I)
            CALL ROTLEV(T,AT,Egrain1,ngrains,B,NDIM,NSYMM)
         ELSEIF ( IDOF(I).EQ.'TOP' ) THEN                         ! ********  Symmetric Top   *********
            NSYMM = NG(I)              ! Symmetry number
            CALL STOPLEV(T,AT,Egrain1,ngrains,B2,B1,NSYMM)
         ELSEIF (IDOF(I).EQ.'KRO') THEN
            B=16.85763d0/RI(I)
            NSYMM=NRSN(I)
            JK=IDIM(I)
            CALL KROTLEV(T,AT,Egrain1,ngrains,B,JK,NSYMM)
         ELSEIF (IDOF(I).EQ.'BOX') THEN
            B=RI(I)
            NDIM=IDIM(I)
            CALL BOXLEV(T,AT,Egrain1,ngrains,B,NDIM,boxzpe)
            zpetot = zpetot + boxzpe
         ELSEIF (IDOF(I).EQ.'VIB') THEN
            WE=RI(I)
            XE=RJ(I)
            CALL MORLEV(T,AT,Egrain1,ngrains,WE,XE,morzpe)
            zpetot = zpetot + morzpe
         ELSE
            write(*,*) '************************************'
            write(*,*) 'FATAL: MISTAKES at convolution step'
            write(*,*) 'Separable Degree of freedom #',I
            write(*,*) '************************************'
            write(4,*) '************************************'
            write(4,*) 'FATAL: MISTAKES at convolution step'
            write(4,*) 'Separable Degree of freedom #',I
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
      DO I=2, ngrains
         SS(I)=SS(I-1) + T(I)
         DS(I)=T(I)/Egrain1
      ENDDO
c --------------------------------------------------------------------------
c     Write to output files
c
      WRITE (5,99030) cut                    ! End of data summary block in ____.qcrp file

      WRITE (5,*) fname(1:lenstr(fname))
      WRITE (5,*) TITLE1
      Nt = 100                               ! Number of temperatures for qvib
      WRITE (5,99901) Egrain1 , Emax2 , Vf, Vr, zpetot, KEYWORD2 , VPTx    ! zpe = coupled & separable (total anharmonic)
99901 FORMAT(5(1x,f10.2),1x,A6,1x,A5)

! About the "Partition function"....
! DS(J) = derivative of cumulative rxn prob. By using DS(J) instead of SS(J), this 
! "partition fxn" incorporates a factor of 1/kT, which cancels the factor of kT in
! the ordinary CTST expression (the kT factor does not appear in SCTST). This 
! cancellation is convenient for use in program THERMO.
! Also...
! Energy zero (for Q) = transition state ZPEanh.
! Starting energy of Q-integral is at DelH-Vf = -Vr or at 0, because delH = Vf-Vr or 0.

      Ntx = 150        ! Number of temperatures to be passed to THERMO
      WRITE(5,*) Ntx
      write(5,*)
      write(5,*) '==================================================='
      WRITE(5,*) ' THERMAL PARTITION FUNCTION',
     &           ' (zero of energy at transition state ZPE)'
      write(5,*) '==================================================='
      write(5,*)
      write(5,*) '   INDEX     T(K)       Q(T)    C/R      H/R      S/R'

        TT=1.0d0
        I = 0
        DO I = 1 , Ntx   ! Number of temperatures to be passed to THERMO
           Q0 = 0.0d0
           Q1 = 0.0d0
           Q2 = 0.0d0
           DO J=1, ngrains 
             enrg = (J-1)*Egrain1 + DelH - Vf
             Ered = 1.4387752D+00*enrg/TT
             fx = DS(J)*exp(-Ered)*Egrain1
             Q0 = Q0 + Fx
             Q1 = Q1 + Ered*Fx
             Q2 = Q2 + Ered*Ered*Fx
           ENDDO   
c            Thermo formulae from Pitzer and Gwinn, J. Chem. Phys. 10, 428 (1942)
           Cx = Q2/Q0 - (Q1/Q0)**2
           Hx = Q1/Q0
           Sx = (Hx+log(Q0))
           write(5,9011) I, TT, Q0 , Cx , Hx , Sx
           TT = TT*1.05631305d+00   ! for Ntx = 150
        ENDDO   

      write(4,*)
      write(4,*) '==================================================='
      WRITE(4,*) '     MICROCANONICAL SUMS AND DENSITIES '
      write(4,*) '==================================================='
      write(4,*)

        ngrains = ngrains - 1   !!! the boundary of interest

      WRITE(4,9017) 
      DO I=1, ngrains
       enrg = Emin + (I-1)*Egrain1     ! zero of energy is at ZPEanh
       WRITE(4,9010) I , enrg , DS(I), SS(I)
      ENDDO 

c --------------------------------------------------------------------------
c --------------------------------------------------------------------------
c
c      Write to ___.crp file (UNIT=3) 
c 

      WRITE (3,99030) cut                    ! End of data summary block in ____.crp file

      WRITE (3,*) fname(1:lenstr(fname))
      WRITE (3,*) TITLE1
      WRITE (3,99001) Egrain1 , imax1 , Emax2 , Isize , Viblo
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
c
c       End of writing to ___.crp file
C
C     CLOSE OUTPUT FILES

      CLOSE (3)
      CLOSE (4)
      CLOSE (5)

c
c --------------------------------------------------------------------------

9911  FORMAT (
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'//,
     &'                      sctst-',A7,/
     &'            (Semi-Classical Transition State Theory) ',//,
     &'                         ',A8,//,
     &'              John R. Barker and Thanh Lam Nguyen',/
     &'          jrbarker@umich.edu     LNGUYEN@cm.utexas.edu *',//
     &'                   University of Michigan',/
     &'               Ann Arbor, Michigan 48109-2143',/
     &'            (present address: U. of Texas, Austin)'//
     &'        http://clasp-research.engin.umich.edu/multiwell',//
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%',//)
9912  FORMAT('Suggested Literature Citations:'//4x,
     &'a) J.R. Barker, T.L. Nguyen, J.F. Stanton, C. Aieta, M. Ceotto,',
     &/7x,'F. Gabas, T.J.D. Kumar, C.G.L.Li, L.L.Lohr, A. Maranzana,',
     &/7x,'N.F. Ortiz, J.M. Preses, J.M. Simmie, J.A. Sonk, and ',
     &'P.J. Stimac; ',
     &/7x,'MultiWell-',A7,' Software Suite; University of Michigan, ',
     &/7x,'Ann Arbor, Michigan, USA, ',A8,'; ',
     &/7x,'http://clasp-research.engin.umich.edu/multiwell/.',
     & //4x,
     &'b) T. L. Nguyen, J. F. Stanton, J. R. Barker, Chem. Phys.',/7x,
     &'Letters, 499, 9-15 (2010).',//4x,
     &'c) W. H. Miller, J. Chem. Phys. 62, 1899-1906 (1975).',//4x,
     &'d) W. H. Miller, Faraday Discuss. Chem. Soc. 62,',/7x,
     &'40-46 (1977).',//4x,
     &'e) W. H. Miller, R. Hernandez, N. C. Handy, D. Jayatilaka, ',/7x,
     &'and A. Willets, Chem. Phys. Letters 172, 62-68 (1990).',//4x,
     &'f) R. Hernandez and W. H. Miller, Chem. Phys. Lett. 214,',/7x,
     &'129-136 (1993).',//4x,
     &'g) F. Wang and D.P. Landau, Phys. Rev. Letters 86, ',/7x,
     &'2050-2053 (2001).',//4x,
     &'h) M. Basire, P. Parneix, and F. Calvo, J. Chem.Phys. 129,',/7x,
     &'081101 (2008).',//4x,    
     &'i) T.L. Nguyen and J.R. Barker, J. Phys. Chem. A., 114,',/7x,
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
9010  FORMAT( 2x,I6,2x, F10.2,5(2x,1pe10.3) )
9011  FORMAT( 2x,I6,2x, F13.5,5(2x,1pe12.5) )
9013  FORMAT(/,' Ave:', 3(1x, F10.3)/)
9015  FORMAT(/'             Partial ZPEanh = ', f10.3,' cm-1',
     &        '  (omits Go from VPT2 and separable D.O.F.)',/)
9016  FORMAT( '             steps up or down:',f10.3,' %  ',/,
     &        '               steps in place:',f10.3,' %  ',/,
     &        '           steps out of range:',f10.3,' %  ',/,
     &        '   average non-zero histogram: ',f10.0,' +/-', f8.0,/,
     &        '  average un-signed deviation:',0pf10.2,' %  ',/,
     &        '       max positive deviation:',0pf10.2,' %  ',/,
     &        '        max negative devation:',0pf10.2,' %',//)
9017  FORMAT(
     & '     No.        cm-1  states/cm-1  SumStates   (E=0 at ',
     &    'reactant ZPE in exothermic direction)')
9020  FORMAT(' Iteration =',I3,'/',I2,5x,'f = ',1pe10.2, 5x,
     &                'p = ',1pe10.4)
9030     FORMAT(I3,2X,'Hind.Rot',2x,'Freq(har in cm-1)=',F8.2,2x,
     &     'Mom(amu.A**2)=',F8.3,2x,'fold=',I2,2x,'symm=',I2,2x,
     &     'B=',F8.4,' cm-1',2x,'Uo=',F7.1,' cm-1')
9031     FORMAT (I3,2X,'General HindRotor:',2x,'Rotor symm. =',I2)
9032     FORMAT (8X,A5,1x,' ; symmV =',I2,2x,'; Phase(rad.) = ',
     & F8.4,2x,'; Coeff. =',20(F10.4,1x))

9033     FORMAT (8X,A5,1x,' ; symmB =',I2,2x,'; Phase(rad.) = ',
     & F8.4,2x,'; Coeff. =',20(F10.4,1x))
99001 FORMAT (1x,f10.1,1x,I10,1x,f10.1,1x,I10,1x,f10.1,1x,f10.2)
99002 FORMAT (/,'T(I) array for coupled degrees of freedom:',/,
     &        '        I      (cm-1)   T(I)')
99010 FORMAT ('       No.     (cm-1)   Density   Sum   (E=0 at ',
     &    'reactant ZPE in exothermic direction)')
99013 FORMAT (/,'TOTAL ZPEanh: ', f10.3,' cm-1 for all coupled & ',
     &  'separable D.O.F.; omits the reaction coordinate and ',
     &  'Go from VPT2.',/)
99014 FORMAT (i10,1x,f10.1,2(1x,1pe10.3))
99023 FORMAT (/A10)
99030 FORMAT (A46)
99041 FORMAT (I3,2X,'Quantized Rotor:  I(Amu.A**2) = ',
     &     F8.3,' ;  Sym. = ',I1,' ;  Dim. = ',I2,' ;  B(cm-1) = ',F8.3)
99042 FORMAT (I3,2X,'K-Rotor:  I(Amu.A**2) = ',F8.3,' ;  Sym. = ',I2,
     &          ' ;  Quantum No. J = ',I2,' ;  B(cm-1) = ',F7.3)
99043 FORMAT (I3,2X,'Particle-in-Box:  Freq(cm-1)  = ',F7.2,
     &             ' ;  Dim. = ',I2)
99044 FORMAT (I3,2X,'Classical Rotor:  I(Amu.A**2) = ',
     &     F8.3,' ;  Sym. = ',I1,' ;  Dim. = ',I2,' ;  B(cm-1) = ',
     &     F8.3)
99045 FORMAT (I3,2X,'Vibration:',8x,'Freq(cm-1)  = ',F7.2,
     &' ;  Anharm. = ', F8.3)
99050 FORMAT (//,1x,I2,' SEPARABLE DEGREES OF FREEDOM',/)
99051 FORMAT (/' --- END INPUT ---')
99060 FORMAT (5x,'REACTION PARAMETERS',/
     & 5x,'Monte Carlo trials: ', 1pe8.1,' (',A6,')',//
     & 5x,'Forward and reverse barriers (+ZPE):   Vf           Vr',/
     &  37x, 0pf10.3,2x,0pf10.3,2x, A6,/
     &  37x, 0pf10.3,2x,0pf10.3,2x,'cm-1',/)
99061 FORMAT (5x,'Imaginary frequency (cm-1)    = ',F10.2,/
     &        5x,'Diagonal anharmonicity (cm-1) = ',F10.2,/)
99062 FORMAT (5x,'Off-Diagonal anharmonicities (cm-1):',//
     & '         CFOUR style    or   GAUSSIAN style',/
     & '      k     X(k,f)      or      i*X(k,f) ',//
     &   50( 5x,I2,1x,F10.2,'i ',9x,F10.2,/) )
99063 FORMAT (5x,A5,': VPT2, only',/)
99064 FORMAT (5x,A5,': VPT2 + VPT4 (Eq. 42 in Stanton, 2016)',/)
99065 FORMAT (5x,A5,': VPT2 + VPT4 (Eq. 37 in Stanton, 2016)',/)
99905 FORMAT (I5,2X,'Symm-Top',2x,'2D-Moment =',F8.2,' amua, B2 =',F8.4,
     &       ' cm-1  ;  1D-Moment =',F8.2,' amua, B1 =',F8.4,' cm-1',2x,
     &       ';  Symm. No. = ',I1)

      END program

