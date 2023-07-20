c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c DenSum: a code for calculating sums and densities of states.
c Copyright (C) 2017 John R. Barker
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

      PROGRAM DenSum  
!
!	Main Program: For computing sums and densities of states via the
!    Beyer-Swinehart algorithm.
!
!    Revised to modernized Fortran by John R. Barker (July 26, 2017)
!
      USE SEPSUBS
      IMPLICIT NONE
      
      REAL(8) :: WOBS
      REAL(8) :: E
      REAL(8) :: Egrain1
      REAL(8) :: Egrain2
      REAL(8) :: Emax1
      REAL(8) :: Emax2
      REAL(8) :: Emin
      REAL(8), DIMENSION(50001) :: SUM
      REAL(8), DIMENSION(50001) :: DENS
      INTEGER :: i
      INTEGER :: ib
      INTEGER :: Imax1
      INTEGER :: Imax2
      INTEGER :: Isize
      INTEGER :: iv
      INTEGER :: IWR
      INTEGER :: j
      INTEGER :: JMAX
      INTEGER :: ki
      INTEGER :: L
      INTEGER :: N
      INTEGER :: NGG
      INTEGER :: ntop
      INTEGER :: NZ
      REAL(8) :: viblo
      REAL(8) :: x
      REAL(8) :: X1
      REAL(8) :: X2
      INTEGER :: lenstr       !Added JP 10/03
      CHARACTER(180) :: TITLE
      CHARACTER(10) :: FNAME, Arg
      CHARACTER(4) :: VHAR
      CHARACTER(6) :: bestImaxTXT
      
      CHARACTER(100) :: readline
      CHARACTER(20) :: filename
      CHARACTER(3), PARAMETER ::  OBS='OBS'
      CHARACTER(3), PARAMETER ::  HAR='HAR'
      INTEGER, DIMENSION(6) ::    bestImax
      CHARACTER(7), PARAMETER ::  AVERSION='2023'
      CHARACTER(12), PARAMETER :: ADATE='Mar 2023'
      CHARACTER(46), PARAMETER :: 
     &            cut='**************INPUT DATA SUMMARY**************'
      CHARACTER(4), DIMENSION(2) :: keyword
      CHARACTER(4) :: VROTIN
      CHARACTER(4) :: VROTOUT

!      EXTERNAL WRAB , WRDEN
 
      call GETARG(1,Arg)
      Arg=Arg(1:lnblnk(Arg))
      IF(Arg.eq."-batch".or.Arg.eq."-BATCH".or.Arg.eq."-Batch") THEN
       Arg="BATCH"
       OPEN(unit=1,file='densum.batch',status='old')
       READ(1,*) Egrain1, Imax1, Isize, Emax2
c       ELSE
c40      write(*,*) "densum.batch not found"
c        STOP
        IF ( Imax1 .GT. Isize-2 ) THEN
         WRITE (*,*) '*** FATAL:  Imax1 must be ≤ (Isize-2) ***'
         WRITE (4,*) '*** FATAL:  Imax1 must be ≤ (Isize-2) ***'
         STOP
        ENDIF      
      ENDIF
 
      DO 

       IF(Arg.eq."BATCH") THEN

       READ(1,'(A100)',END=10) readline
       filename=readline(1:lnblnk(readline))
       IF(readline(1:lnblnk(readline)).eq.'') goto 450
       OPEN(unit=5,file=filename,err=50,status='old')  
       GOTO 51
50     WRITE(*,*) 
       WRITE(*,*) "File ",filename(1:lnblnk(filename))," not found"
       WRITE(*,*) 
       GOTO 450
51     OPEN(unit=2,file='densum.dat',status='OLD')
       DO i=1,3
        READ(5,'(A100)') readline
        WRITE(2,*) readline(1:lnblnk(readline))
       ENDDO
       READ(5,*)
       WRITE(2,*) Egrain1, Imax1, Isize, Emax2
       DO 
        READ(5,'(A100)',END=20) readline
        WRITE(2,*) readline(1:lnblnk(readline))
       ENDDO
20    CLOSE(5)
      REWIND(2)

       ELSE

       OPEN (UNIT=2,STATUS='OLD',FILE='densum.dat')              ! Data file 

      ENDIF

      OPEN (3,STATUS='UNKNOWN',FILE='densum.out')  ! Full output
      OPEN (12,STATUS='UNKNOWN',FILE='densum.lev')   ! Energy levels for vibrations, rotations, and hindered rotors
 
      WRITE (3,99016) AVERSION , ADATE ,  AVERSION , ADATE
      CALL DateTime(3)
      READ (2,99017) TITLE
      WRITE (3,99017) TITLE
      WRITE (12,99017) TITLE
      READ (2,*) FNAME
      WRITE (3,99023) FNAME
      WRITE (12,99023) FNAME
      
      OPEN(4,STATUS='UNKNOWN',FILE=FNAME(1:lenstr(fname))//'.dens')      ! Succinct output for MultiWell input

      WRITE (4,99030) cut                                                ! Start of data summary block in ____.dens file

      WRITE (4,99016) AVERSION , ADATE ,  AVERSION , ADATE
      WRITE (12,99016) AVERSION , ADATE ,  AVERSION , ADATE
      CALL DateTime(4)

      READ (2,*) N , IWR, (keyword(i),i=1,2) ! N=no. of DoF, IWR=0 for direct count, Harmonic or observed (fundamental) frequencies
      IWR = 0  ! enforce dirct count; the Whitten-Rabinovitch option has been deleted.
         CALL ucase ( keyword(1) )   ! convert to upper case
         CALL ucase ( keyword(2) )   ! convert to upper case
c
c    KEYWORDS:
c    VHAR = 'HAR' for harmonic vibs, 
c         = 'OBS' or 'FUN' FOR observed (0-1 fundamental transitions)
c    Units for qro and rot types [NOT hindered rotor types]:
c    VROTIN = 'AMUA' for moment of inertia in units of amu*Ang^2
c         = 'GMCM' for moment of inertia units of gram*cm^2
c         = 'CM-1' for rotational constant in units of cm^-1
c         = 'MHZ' for rotational constant in units of MHz
c         
         DO ki = 1 , 2
           IF ( keyword(ki) .EQ. 'HAR' ) THEN
             VHAR = 'HAR'
           ELSEIF ( keyword(ki) .EQ. 'OBS' .OR. 
     &              keyword(ki) .EQ. 'FUN' ) THEN
             VHAR = 'OBS'
           ELSEIF ( keyword(ki) .EQ. 'AMUA' )THEN
             VROTIN = 'AMUA'
           ELSEIF ( keyword(ki) .EQ. 'GMCM' )THEN
             VROTIN = 'GMCM'
           ELSEIF ( keyword(ki) .EQ. 'CM-1' )THEN
             VROTIN = 'CM-1'
           ELSEIF ( keyword(ki) .EQ. 'MHZ'  )THEN
             VROTIN = 'MHZ'
           ELSEIF ( keyword(ki) .EQ. 'GHZ' )THEN
             VROTIN = 'GHZ'
           ELSE
             write(*,*) 'FATAL: unknown keyword: ', keyword(ki)
             STOP
           ENDIF
         END DO

      IF(Arg.eq."BATCH") THEN
       READ (2,*) 
      ELSE
       READ (2,*) Egrain1 , imax1 , Isize , Emax2
      ENDIF
 
      IF (Emax2/Egrain1 .GT. NMAX) THEN
        WRITE(*,*) '****Grain Size Too Small****'
        STOP
      ENDIF

      Emax1 = Egrain1*(imax1-1)
      Egrain2 = Emax2/(Isize-imax1-1)
      imax2 = Isize - imax1
 
      JMAX = 1 + INT(Emax2/Egrain1)
      JMAX = MIN(JMAX,NMAX)
 
      WRITE (3,99008)
      WRITE (4,99008)
 
      VROTOUT = 'AMUA'    ! unless specified otherwise
      Viblo = 10000.
      NZ = 0 
      ntop = 0            ! for counting the number of symmetric tops (should never be >1)
      DO 100 I = 1 , N
         READ (2,*) MODE(I) , IDOF(I) , WE(I) , ANH(I) , NG(I)

         CALL ucase ( IDOF(I) )   ! convert to upper case
         IF ( WE(I) .LT. 0.0 ) THEN
           OPEN( 7, status='UNKNOWN', file='densum.ERROR')
           WRITE(*,*) '***FATAL: WE(i) negative value: ', WE(I)
           WRITE(7,*) '***FATAL: WE(I) negative value: ', WE(I)
           STOP
         ENDIF
     
         IF ( IDOF(I).EQ.VIB ) THEN                                   ! Vibration: Evib = WE(v+1/2) + ANH*(v+1/2)^2
            Viblo = MIN(Viblo,WE(I))         ! lowest vib frequency
            IF ( VHAR .EQ. HAR ) THEN
              WOBS = WE(I)+2.*ANH(I)         ! WOBS = observed freq
            ELSEIF ( VHAR .EQ. OBS ) THEN
              WOBS = WE(I)
              WE(I)= WE(I)-2.*ANH(I)         ! Convert WE to harmonic frequency
            ENDIF
            WRITE (3,99002) MODE(I) , WOBS , WE(I) , ANH(I) , NG(I)
            WRITE (4,99002) MODE(I) , WOBS , WE(I) , ANH(I) , NG(I)
            WRITE (12,99002) MODE(I) , WOBS , WE(I) , ANH(I) , NG(I)

         ELSEIF ( IDOF(I).EQ.BOX ) THEN                               ! Particle-in-a-Box Vibration: Evib = WE*v^2
            Viblo = MIN(Viblo,WE(I))         ! lowest vib frequency
            IF ( VHAR .EQ. HAR ) THEN
              WOBS = 3.d+00*WE(I)            ! WOBS = observed freq
            ELSEIF ( VHAR .EQ. OBS ) THEN
              WOBS = WE(I)
              WE(I)= WE(I)/3.0d+00           ! Convert WE to frequency
            ENDIF
            WRITE (3,99022) MODE(I) , WOBS , WE(I) , NG(I)
            WRITE (4,99022) MODE(I) , WOBS , WE(I) , NG(I)
            WRITE (12,99022) MODE(I) , WOBS , WE(I) , NG(I)

         ELSEIF ( IDOF(I).EQ.KRO ) THEN                                 ! 1-dim. K-rotor
            CALL ucase ( IDOF(I) )   ! convert to upper case
            AMOM(I) = WE(I)                     ! moment of inertia or rotational constant for 2D rotor (J-rotor)
            call rotunits( VROTIN , AMOM(I) , 'AMUA')     ! unit conversion to 
            B(I) = 16.85763D+00/AMOM(I)         ! rotational constant (cm-1) for 2D rotor
            B2 = WE(I)
            X2 = WE(I)
            call rotunits( VROTIN , B2 , 'CM-1' )      ! rotational constant (cm-1) for 2D rotor
            call rotunits( VROTIN , X2 , 'AMUA' )      ! moment of inertia (AMUA) for 2D rotor
            call rotunits( VROTIN , ANH(I) , 'CM-1' )      ! rotational constant (cm-1) for K-rotor
            B1 = ANH(I)                                ! for K-rotor
            X1 = 16.85763D+00/ANH(I)                       ! moment of inertia for K-rotor
            WRITE (3,99051) MODE(I) , AMOM(I) , B2 ,  X1 , B1 , NG(I)
            WRITE (4,99051) MODE(I) , AMOM(I) , B2 ,  X1 , B1 , NG(I)
            WRITE (12,99051) MODE(I) , AMOM(I) , B2 ,  X1 , B1 , NG(I)

         ELSEIF ( IDOF(I).EQ.QRO ) THEN                                 ! Quantized rotor
            AMOM(I) = WE(I)                     ! moment of inertia or rotational constant
            call rotunits( VROTIN , AMOM(I) , VROTOUT )     ! unit conversion
            B(I) = 16.85763D+00/AMOM(I)         ! rotational constant (cm-1)
            WRITE (3,99005) MODE(I) , AMOM(I) , ANH(I) , NG(I) , B(I)
            WRITE (4,99005) MODE(I) , AMOM(I) , ANH(I) , NG(I) , B(I)
            WRITE (12,99005) MODE(I) , AMOM(I) , ANH(I) , NG(I) , B(I)
            IF ( ANH(I).LE.0 ) THEN
               WRITE (3,99006)
               WRITE (4,99006)
               STOP
            ENDIF

         ELSEIF ( IDOF(I).EQ.TOP ) THEN                                 ! Symmetric Top with 2D and 1D rot constants B2 and B1
            ntop = ntop + 1
            IF ( ntop .GT. 1) THEN
              write (*,*) 'FATAL: only one symmetric top (top) allowed'
              write (3,*) 'FATAL: only one symmetric top (top) allowed'
              write (4,*) 'FATAL: only one symmetric top (top) allowed'
            ENDIF
            B2 = WE(I)
            X2 = B2
            call rotunits( VROTIN , B2 , 'CM-1' )     ! unit conversion
            call rotunits( VROTIN , X2 , 'AMUA' )     ! unit conversion
            B1 = ANH(I)
            X1 = B1
            call rotunits( VROTIN , B1 , 'CM-1' )     ! unit conversion
            call rotunits( VROTIN , X1 , 'AMUA' )     ! unit conversion
            WRITE (3,99905) MODE(I) , X2 , B2 , X1 , B1 , NG(I)
            WRITE (4,99905) MODE(I) , X2 , B2 , X1 , B1 , NG(I)
            WRITE (12,99905) MODE(I) , X2 , B2 , X1 , B1 , NG(I)

         ELSEIF ( IDOF(I).EQ.ROT ) THEN                                 ! Classical rotor
            AMOM(I) = WE(I)                     ! moment of inertia or rotational constant
            call rotunits( VROTIN , AMOM(I) , VROTOUT )     ! unit conversion
            B(I) = 16.85763D+00/AMOM(I)         ! rotational constant (cm-1)
            WRITE (3,99003) MODE(I) , AMOM(I) , ANH(I) , NG(I) , B(I)
            WRITE (4,99003) MODE(I) , AMOM(I) , ANH(I) , NG(I) , B(I)
            WRITE (12,99003) MODE(I) , AMOM(I) , ANH(I) , NG(I) , B(I)
            IF ( ANH(I).LE.0 ) THEN
               WRITE (3,99006)
               WRITE (4,99006)
               STOP
            ENDIF

         ELSEIF ( IDOF(I).EQ.HRA ) THEN         ! INPUT MOMENT & VIB FREQ   Hindered rotor
            IF ( NG(I) .LT. 0) THEN
                READ(2,*) NSIG(I)
                NG(I) = -NG(I)
              ELSE
                NSIG(I) = NG(I)
            ENDIF
            W(I) = WE(I)                        ! small vibration frequency
            Viblo = MIN(Viblo,W(I))             ! lowest vib frequency
            AMOM(I) = ANH(I)                    ! moment of inertia
            call rotunits( VROTIN , AMOM(I) , VROTOUT)     ! unit conversion
            B(I) = 16.85763D+00/AMOM(I)         ! rotational constant (cm-1)
            VV(I) = ((W(I)/NG(I))**2)/B(I)      ! hindrance barrier (cm-1)
            WRITE (3,99004) MODE(I) , W(I) , AMOM(I) , NG(I) , NSIG(I) , 
     &             B(I) , VV(I)                         ! Hindered rotor (extended Troe method)
            WRITE (4,99004) MODE(I) , W(I) , AMOM(I) , NG(I) , NSIG(I) , 
     &             B(I) , VV(I)                         ! Hindered rotor (extended Troe method)
            WRITE(12,99004) MODE(I) , W(I) , AMOM(I) , NG(I) , NSIG(I) , 
     &             B(I) , VV(I)                         ! Hindered rotor (extended Troe method)
            IF ( ANH(I).LE.0 .OR. NG(I).EQ.0 ) THEN
               WRITE (3,99007)
               WRITE (4,99007)
               STOP
            ENDIF

         ELSEIF ( IDOF(I).EQ.HRB ) THEN         ! INPUT VIB FREQ & BARRIER   Hindered rotor
            IF ( NG(I) .LT. 0) THEN
                READ(2,*) NSIG(I)
                NG(I) = -NG(I)
              ELSE
                NSIG(I) = NG(I)
            ENDIF
            W(I) = WE(I)                        ! small vibration frequency
            Viblo = MIN(Viblo,W(I))             ! lowest vib frequency
            VV(I) = ANH(I)                      ! hindrance barrier (cm-1)
            IF (VV(I) .LT. 1.0d-05) THEN
              WRITE(*,*)
     &        'Rotation barrier too small. Hit RETURN to terminate' 
              STOP
            ENDIF

            B(I) = ((W(I)/NG(I))**2)/VV(I)      ! rotational constant (cm-1)
            AMOM(I) = 16.85763D+00/B(I)         ! moment of inertia (amu*ang^2)

            WRITE (3,99004) MODE(I) , W(I) , AMOM(I) , NG(I) , NSIG(I) ,
     &             B(I) , VV(I)                         ! Hindered rotor (extended Troe method)
            WRITE (4,99004) MODE(I) , W(I) , AMOM(I) , NG(I) , NSIG(I) ,
     &             B(I) , VV(I)                         ! Hindered rotor (extended Troe method)
            WRITE(12,99004) MODE(I) , W(I) , AMOM(I) , NG(I) , NSIG(I) ,
     &             B(I) , VV(I)                         ! Hindered rotor (extended Troe method)
            IF ( ANH(I).LE.0 .OR. NG(I).LE.0 ) THEN
               WRITE (3,99007)
               WRITE (4,99007)
               STOP
            ENDIF

         ELSEIF ( IDOF(I).EQ.HRC ) THEN         ! INPUT MOMENT & BARRIER   Hindered rotor
            IF ( NG(I) .LT. 0) THEN
                READ(2,*) NSIG(I)
                NG(I) = -NG(I)
              ELSE
                NSIG(I) = NG(I)
            ENDIF
            AMOM(I) = WE(I)                     ! moment of inertia (amu*ang^2)
            call rotunits( VROTIN , AMOM(I) , VROTOUT)     ! unit conversion
            B(I) = 16.85763D+00/AMOM(I)         ! rotational constant (cm-1)
            VV(I) = ANH(I)                      ! hindrance barrier (cm-1)
            W(I) = NG(I)*SQRT(B(I)*VV(I))       ! small vibration frequency
            Viblo = MIN(Viblo,W(I))             ! lowest vib frequency
            WRITE (3,99004) MODE(I) , W(I) , AMOM(I) , NG(I) , NSIG(I) ,
     &             B(I) , VV(I)                 ! Hindered rotor (extended Troe method)
            WRITE (4,99004) MODE(I) , W(I) , AMOM(I) , NG(I) , NSIG(I) ,
     &             B(I) , VV(I)                         ! Hindered rotor (extended Troe method)
            WRITE(12,99004) MODE(I) , W(I) , AMOM(I) , NG(I) , NSIG(I) ,
     &             B(I) , VV(I)                         ! Hindered rotor (extended Troe method)
            IF ( ANH(I).LE.0 .OR. NG(I).EQ.0 ) THEN
               WRITE (3,99007)
               WRITE (4,99007)
               STOP
            ENDIF

      ELSEIF ( IDOF(I).EQ.HRD ) THEN                                    ! General, unsymmetrical hindered rotor
        NSIG(I)=NG(I)         ! rotational symmetry number
        NCV(I)=INT(WE(I))     ! number of coefficients in rot. potential energy function  
        NCB(I)=INT(ANH(I))    ! number of coefficients in rot. constant
 
        READ(2,*) Vhr(I), NSV(I), Phav(I), (CV(I,Iv), Iv=1, NCV(I))     ! always cm-1 units
        READ(2,*) Bhr(I), NSB(I), Phab(I), (CB(I,Ib), Ib=1, NCB(I))
        call ucase ( Vhr(I) )
        call ucase ( Bhr(I) )

        IF( Bhr(I) .EQ. 'BHRD1' ) THEN                                  ! always cm-1 units
          WRITE (3,99040) MODE(I), Vhr(I), Bhr(I) 
          WRITE(3,99041) NSV(I), Phav(I), (CV(I,Iv), Iv=1, NCV(I))
          WRITE(3,99042) NSB(I), Phab(I), (CB(I,Ib), Ib=1, NCB(I))

          WRITE (4,99040) MODE(I), Vhr(I), Bhr(I)
          WRITE(4,99041) NSV(I), Phav(I), (CV(I,Iv), Iv=1, NCV(I))
          WRITE(4,99042) NSB(I), Phab(I), (CB(I,Ib), Ib=1, NCB(I))

          WRITE (12,99040) MODE(I), Vhr(I), Bhr(I)
          WRITE(12,99041) NSV(I), Phav(I), (CV(I,Iv), Iv=1, NCV(I))
          WRITE(12,99042) NSB(I), Phab(I), (CB(I,Ib), Ib=1, NCB(I))
        ELSEIF( Bhr(I) .EQ. 'IHRD1' ) THEN                              ! always amu*Ang^2 units
          WRITE (3,99040) MODE(I), Vhr(I), Bhr(I)
          WRITE(3,99041) NSV(I), Phav(I), (CV(I,Iv), Iv=1, NCV(I))
          WRITE(3,99043) NSB(I), Phab(I), (CB(I,Ib), Ib=1, NCB(I))

          WRITE (4,99040) MODE(I), Vhr(I), Bhr(I)
          WRITE(4,99041) NSV(I), Phav(I), (CV(I,Iv), Iv=1, NCV(I))
          WRITE(4,99043) NSB(I), Phab(I), (CB(I,Ib), Ib=1, NCB(I))

          WRITE (12,99040) MODE(I), Vhr(I), Bhr(I)
          WRITE(12,99041) NSV(I), Phav(I), (CV(I,Iv), Iv=1, NCV(I))
          WRITE(12,99043) NSB(I), Phab(I), (CB(I,Ib), Ib=1, NCB(I))
        ELSE
          write(*,*) "ERROR at INPUT for Bhr/Ihr"
        ENDIF
	     DO L=1, NCV(I)
	       CVt(L)=CV(I,L)
	     ENDDO
          DO L=1, NCB(I)
            CBt(L)=CB(I,L)
          ENDDO
	     NGG = NG(I)         ! Hinder rotation symmetry No. (foldedness)

        ELSEIF ( IDOF(I) .EQ. TRN ) THEN                                ! RELATIVE TRANSLATION
          AMOM(I) = WE(I)*ANH(I)/( WE(I) + ANH(I) )     ! reduced mass (amu)
          WRITE (3,99021) MODE(I) , WE(I) , ANH(I), AMOM(I)
          WRITE (4,99021) MODE(I) , WE(I) , ANH(I), AMOM(I)
          IF ( NG(I).LE.0 .OR. NG(I).GT.3 ) THEN
               WRITE (3,99022)
               STOP
          ENDIF

        ELSE
          write( 3,*) IDOF(I)
          write( 4,*) IDOF(I)
          WRITE(*,*) 'FATAL ERROR: unrecognized dof-type: ',IDOF(I)

        ENDIF   ! end of types of degrees of freedom
 
 100  CONTINUE
 
      WRITE (4,99030) cut       ! End of data summary block in ____.dens file
      
      ZPPE = 0     
      HRZPE = 0.0d+00
C
C     CALCULATE SUMS AND DENSITIES
C
      WRITE (*,99019)
      Emin = 0.0              ! zero of energy VZPE-Emin (Emin = 0 for all dof except for type KRO)
      CALL STERAB( N, JMAX, Egrain1, SUM, DENS )
c
c
c      Check for density fluctuations
c
      call chkdens(Egrain1,Imax1,Isize,DENS,Emax2,bestImax)
c
c
c     Write out results
c
c
c      Succinct output for MultiWell input
c 
      WRITE (4,*) FNAME(1:lenstr(fname))
      IF(Imax1.lt.bestImax(5)) THEN
       WRITE (4,'(A15,A49,I6)') TITLE(1:lnblnk(TITLE)),
     & "  --  WARNING. Miminum suggested value for Imax1:",bestImax(5)
      ELSE
       WRITE (4,'(A50)') TITLE
      ENDIF
      WRITE (4,99001) Egrain1 , imax1 , Emax2 , Isize , Viblo
 
      WRITE (3,99012) Egrain1   ! Exact Counts

      WRITE (3,99020) ZPPE , HRZPE 

      x=0
      WRITE(3,*) " Density fluctuations %       Imax1       E(cm-1)"
      DO i=1,5
         if (bestImax(i).ne.INT(Emax2/Egrain1)) then
          E=(bestImax(i)-1)*Egrain1
         else
          E=Emax2
         endif
!      IF(bestImax(i).ne.0) THEN
        IF(Imax1.ge.bestImax(i).AND.x.lt.1) THEN
         x=x+1
         WRITE(3,'(A10,I2,A18,I6,A6,F8.1,A12)') "          ",i,
     &"%                 ", bestImax(i), "      ", E, "  <- Current"
          ELSE
           IF(bestImax(i).eq.INT(Emax2/Egrain1)) then
           WRITE(3,'(A10,I2,A18,I6,A6,F8.1,A12)') "          ",i,
     &"%                >", bestImax(i), "     >", E
           ELSE
           WRITE(3,'(A10,I2,A18,I6,A6,F8.1,A12)') "          ",i,
     &"%                 ", bestImax(i), "      ", E
           ENDIF 
        ENDIF
!      ENDIF
      ENDDO
       WRITE(3,*)
      IF(Imax1.lt.bestImax(5)) THEN
       WRITE (3,*) "WARNING: the current density fluctuation at Imax1 is
     & larger the 5%."
       WRITE (3,'(A35,I6)') " Miminum suggested value for Imax1:",
     & bestImax(5)
       WRITE(3,*)
       DO i=1,5
        IF(bestImax(i).eq.INT(Emax2/Egrain1).AND.
     &   bestImax(i+1).le.INT(Emax2/Egrain1)) then 
       WRITE(3,'(A57,I2,A2)')" WARNING: the density fluctuation at Emax2
     & is larger than",i," %."
         WRITE (3,*) "Increase Emax2"
        ENDIF
       ENDDO
       WRITE(3,*)
       WRITE(3,*)
      ENDIF
      
      CALL DateTime(3)
 
      WRITE (3,99011) Emin
      WRITE (4,99010) Emin
 
      DO 200 i = 1 , imax1                          ! Lower part of double array; energy at TOP of grain
         E = (i-1)*Egrain1
         WRITE (4,99015) i , E , DENS(i) , SUM(i)
 200  CONTINUE
 
      DO 300 i = imax1 + 1 , Isize                  ! Upper part of double array; energy at TOP of grain
         E = (i-imax1-1)*Egrain2
         j = INT(E/Egrain1) + 1
         WRITE (4,99015) i , E , DENS(j) , SUM(j)
 300  CONTINUE
 
      DO 400 i = 1 , JMAX                           ! write to densum.out
         E = (i-1)*Egrain1
         WRITE (3,99015) i , E , DENS(i) , SUM(i)
 400  CONTINUE
 
       write(12,*) '---------------------------------------------------'

 
      IF(Arg.ne."BATCH") GOTO 500
     
      CLOSE(2) !  data input file
      CLOSE(4) ! '.dens' output file
      ENDDO

450   write(*,*)

10    CLOSE(1)

500   CONTINUE
      CLOSE(3)
      CLOSE(12)

!   FORMAT STATEMENTS
!
99001 FORMAT (1x,f10.1,1x,I10,1x,f10.1,1x,I10,1x,f10.1)
99002 FORMAT (I3,2X,'Vibrator: Fund= ',F8.2,2x,
     &        'Harm= ',F8.2,2x,'Anh(cm-1)= ',F7.2,2X,'g= ',I2)
99022 FORMAT (I3,2X,'P-in-Box',2x,'Freq(0-1)=',F8.2,2x,
     &        'Freq.Param.=',F8.2,2x,'g =',I2)
99003 FORMAT (I3,2X,'Clas.Rot',2x,'AmuAng^2 =',F8.2,2x,'Symm =',F4.0,2x,
     &        ' dim =',I2,2x,'B =',F8.4,' cm-1')
99004 FORMAT (I3,2X,'Hind.Rot',2x,'Freq(har)=',F8.2,2x,'Mom=',F8.3,2x,
     &       'fold=',I2,2x,'symm=',I2,2x,'B=',F8.4,' cm-1',2x,
     &        'Uo=',F7.1,' cm-1')
99005 FORMAT (I3,2X,'Quan.Rot',2x,'AmuAng^2 =',F8.2,2x,'Symm =',F4.0,2x,
     &        ' dim =',I2,2x,'B =',F8.4,' cm-1')
99905 FORMAT (I3,2X,'Symm-Top',2x,'2D-Moment =',F8.2,' amua, B2 =',F8.4,
     &       ' cm-1  ;  1D-Moment =',F8.2,' amua, B1 =',F8.4,' cm-1',2x,
     &       ';  Symm. No. = ',I1)
99006 FORMAT (//'Rotor symmetry number must be greater than zero')
99007 FORMAT (//'Hindered rotor parameters incomplete')
99008 FORMAT (//'Degrees of Freedom'/)
99009 FORMAT (/'  No.  Freq(cm-1)  ANH(cm-1)      Degeneracy',/)
99010 FORMAT ('       No.    E-Emin     Density     Sum   ',
     &     '[Emin = ',F6.1,' cm-1; E at TOP of energy grains]')
99011 FORMAT ('  Note: E at TOP of energy grains',/,
     &        '         Emin = ',F6.1,' cm-1', / ,
     &        '       No.    E-Emin     Density     Sum')
99012 FORMAT (//'Stein-Rabinovitch Method (grain =',F6.1,' cm-1)'//)
99013 FORMAT (//'Whitten-Rabinovitch Method'//)
99014 FORMAT (i10,1x,f10.1,2(1x,es10.3))
99015 FORMAT (i10,1x,f10.1,2(1x,es12.5))
99016 FORMAT (
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'//8x,
     &'                                John R. Barker'/8x,
     &'        DenSum-',A6,'           University of Michigan'/8x,
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
     & //4x,'b) J.R. Barker, Int. J. Chem. Kinetics, 33, 232-45 (2001)',
     &//,
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'/)
99017 FORMAT (A180)
99018 FORMAT (A10)
99019 FORMAT (/'Calculating Densities')
99020 FORMAT (/' Total Zero Point energy (ZPE): ',f12.3,' cm-1',/
     &         '    ZPE from Hindered Rotor(s): ',f12.3,' cm-1'//)
99021 FORMAT (I3,2X,'Translat',2x,'MolWt(A) =',F8.2,2x,'MolWt(B) =',
     &           F8.2,2x,'RedMass =',F8.2,1x,'g/mol')
99023 FORMAT (/A10)
99030 FORMAT (A46)
99040 FORMAT (I3,2X,'General Hind.Rot:',2x,
     &     'Vhr_Type = ',A5,';',1x,'Bhr/Ihr_Type = ',A5)
99041 FORMAT (5X,'Vhr: ',1x,'Symm. = ',I2,2x,'Phase (rad.) = ',F7.4,
     &     2x,'Coeff. (cm-1) = ',20(F10.4,1x))
99042 FORMAT (5X,'Bhr: ',1x,'Symm. = ',I2,2x,'Phase (rad.) = ',F7.4,
     &     2x,'Coeff. (cm-1) = ',20(F10.4,1x))
99043 FORMAT (5X,'Ihr: ',1x,'Symm. = ',I2,2x,'Phase (rad.) = ',F7.4,
     &  2x,'Coeff. (amu.A**2) = ',20(F10.4,1x))
99045 FORMAT(8X, 'Hindered rotor Zero Point Energy =',F7.1,' cm-1')
99051 FORMAT (I3,2X,'Kro type',2x,'2D-Moment =',F8.2,' AMUA, B2 =',F8.4,
     &       ' cm-1  ;  1D-Moment =',F8.2,' AMUA, B1 =',F8.4,' cm-1',2x,
     &       ';  JK = ',I4)

      END PROGRAM DenSum
