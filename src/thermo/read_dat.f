c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Copyright (C) 2009 Andrea Maranzana and John R. Barker
c
c Andrea Maranzana
c andrea.maranzana@unito.it
c Department of General and Organic Chemistry
c University of Torino
c Corso Massimo D'Azeglio, 48
c Torino  10125
c ITALY
c ++39-011-670-7637
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

!     READ THERMO.DAT  AND CALC ZPE

      subroutine read_dat

      include 'declare.inc'
      REAL(8) Mass(2), dum , amass , Hfxn, Tmax , VVr
      CHARACTER(len=20) NameFile , NameOUT
      CHARACTER(len=15) readpot
      CHARACTER(len=49) formula
      CHARACTER(len=1) ADUM
      CHARACTER(len=4) keyword , VROT
      DIMENSION keyword(2)
      
      CHARACTER(len=46) cut
      PARAMETER (cut='**************INPUT DATA SUMMARY**************')
      CHARACTER(len=33) cut2 , dumcut                                   ! shorter
      PARAMETER (cut2='**************INPUT DATA SUMMARY*')              ! shorter
      INTEGER(4) ki, istat , ntop
      REAL(8) Vo, Vm
      SAVE

      hindrance=0
      Gorin=0
      c=0
      a1=0
c
c     FOR USER-DEFINED FILENAME ("NameFile")
c
      CALL get_command_argument( 1, NameFile )
      NameFile=NameFile(1:len_trim(NameFile))
      IF(NameFile.eq.'') NameFile="thermo.dat"              ! default filename
      NameOUT=NameFile(1:len_trim(NameFile)-4)//".out"
      
      OPEN (UNIT=4,STATUS='OLD',FILE=NameFile)              ! Input Data file
      OPEN (UNIT=3,STATUS='UNKNOWN',FILE=NameOUT)           ! Full output

      READ (4,*) Eunits , Sunits
      CALL ucase ( Eunits )   ! convert to upper case
      CALL ucase ( Sunits )   ! convert to upper case

       IF ( Eunits.EQ.'KCAL' ) THEN
         Rgas = 8.314472D+00/4.184d+00
      ELSEIF ( Eunits.EQ.'KJOU' ) THEN
         Rgas = 8.314472D+00
      ELSEIF ( Eunits.EQ.'CM-1' ) THEN
         Rgas = 0.6950356D+00
      ELSE
         WRITE(3,*) '*** FATAL: Unknown Energy Units ***'
         write(*,*)  '*** FATAL: Unknown Energy Units ***'
         STOP  '*** FATAL: Unknown Energy Units ***'
      ENDIF
 
      IF ( Sunits.EQ.'ATM' ) THEN
         UNITS = 'atmosphere'
      ELSEIF ( Sunits.EQ.'BAR' ) THEN
         UNITS = 'bar'
      ELSEIF ( Sunits.EQ.'MCC' ) THEN
         UNITS = 'molecule/cc'
      ELSEIF ( Sunits.EQ.'NOT' ) THEN
         UNITS = 'no_translat.'
      ELSE
         WRITE(3,*) '*** FATAL: Unknown Standard State ***'
         write(*,*)  '*** FATAL: Unknown Standard State ***'
         STOP  '*** FATAL: Unknown Standard State ***'
      ENDIF

      IF ( Eunits .EQ. 'CM-1' ) THEN
        DIV = 1.d+00
      ELSE
        DIV = 1000.d+00		! for kilo-units: kJ and kcal
      ENDIF

c
c      Number of temperatures (K)
c
      READ (4,*) Nt
      READ (4,*) (T(j),j=1,Nt)  ! Read list of temperatures
      
      Tmax=T(1)
      DO j=2, Nt
            IF(T(j).GT.Tmax) Tmax=T(j)      ! To find Tmax (maximum of Temperature) 
      ENDDO
      IF(Tmax.LT.3000.0d0) THEN
              Tmax=3000.0d0
      ENDIF      
      Emax=10.0d0*(0.6950356D+00)*Tmax      ! Emax = 10RT

      Hdiff = 0.0
      Nreac = 0
      Nprod = 0
      Nctst = 0
 
      READ (4,*) Ns
c
c      Number of species
c
      DO i = 1 , Ns
         Warning(i)=0
         nso(i) = 0
         READ (4,*) REPROD(i) , MOLNAME(i) , DelH(i)   ! reac/prod, name, Enthalpy(0 K)
         CALL ucase ( REPROD(i) )   ! convert to upper case

         IF ( REPROD(i) .EQ. CTST ) THEN
           BACKSPACE(UNIT=4)
           READ(4,*,IOSTAT=istat) REPROD(i) ,MOLNAME(i) ,
     &           DelH(i) , vimag, VVr  ! the usual + imaginary frequency and barrier height for reverse rxn (in Eunits)
           CALL ucase ( REPROD(i) )   ! convert to upper case
           IF ( vimag .LT. 0.0 ) istat = +1
           IF (   VVr .LT. 0.0 ) istat = +1
           IF ( istat .ne. 0 ) then
             WRITE(*,*)  '*** FATAL ERROR in CTST input data ***'
             WRITE(*,*)  '    This line should contain:'
             WRITE(*,*)  '    ''CTST'' , MOLNAME , DelH , vimag, VVr'
             WRITE(*,*)  ' *** vimag and VVr should be positive*** '
             WRITE(*,*)  ' '
             WRITE(3,*)  '*** FATAL ERROR in CTST input data ***'
             WRITE(3,*)  '    This line should contain:'
             WRITE(3,*)  '    ''CTST'' , MOLNAME , DelH , vimag, VVr'
             WRITE(3,*)  ' *** vimag and VVr should be positive*** '
            STOP
           ENDIF

         ENDIF
         
         READ (4,*) Formula                            ! empirical chemical formula
         CALL ChemFormula(i,Formula)
c
c   empirical formula
c
         AMU(i) = 0.0d+00
         Hfxn298(i) = 0.0d+00
         DO ne = 1 , Nelement(i)                                    ! Empirical formula
            CALL element( ATYPE(i,ne) , amass , Hfxn , Sfxn)
            AMU(i) = AMU(i) + Natom(i,ne)*amass                     ! Atomic masses
            Hfxn298(i) = Hfxn298(i) + Natom(i,ne)*Hfxn*1000.d+00    ! [H(298.15) - H(0)]/Rgas (enthalpy function)
            Sfel(i) = Sfel(i) + Natom(i,ne)*Sfxn                    ! Sf/Rgas entropy of elements at 298.15 K
         END DO
         Hfxn298(i) = Hfxn298(i)*Rgas
         Sfel(i) = Sfel(i)*Rgas
c
c   Comment lines in data file
c
         ADUM = '!'
         line = 0
         DO WHILE ( ADUM .EQ. '!' )
            READ (4, *) ADUM
            IF ( ADUM .EQ. '!' ) THEN
               line = line + 1
               BACKSPACE (4)
               READ (4,99114) TITLELINE( i , line )
            ENDIF
         END DO
         nlines( i ) = line
         BACKSPACE (4)

99114 FORMAT (A150)

c
c   symmetry number, optical isomers, number of electronic levels
c               READ (4,99114) TITLELINE( i , nlines(i)+1 )
c               WRITE(*,*) TITLELINE( i , nlines(i)+1 )
         READ (4,*) Sym(i) , Sopt(i) , Nele(i) 

         IF ( Nele(i) .GE. 1 ) THEN
           DO ki = 1 , Nele(i)
            READ (4,*) Elev(i,ki) , gele(i,ki)
           END DO
         ENDIF
c
c   number of internal degrees of freedom (rots + vibs) to be read,
c      and two keywords
         READ (4,*) N(i) , keyword(1) , keyword(2)
         CALL ucase ( keyword(1) )   ! convert to upper case
         CALL ucase ( keyword(2) )   ! convert to upper case
         
         IF ( N(i) .GT. Maxvibs ) THEN
            Write(*,*) '*** FATAL: number of internal DOF > Maxvibs ***'
            write(*,*) '   (suggest increase Maxvibs in declare.inc)'
            write(*,*) ' '
            STOP
         ENDIF
c     KEYWORDS:
c          VHAR = 'HAR' for harmonic vibs, 
c       or VHAR = 'OBS' or 'FUN' FOR observed (0-1 transitions)
c
c      for qro and rot types, including hindered rotor types EXCEPT hrd:
c          VROT = 'AMUA' for moment of inertia in units of amu*Ang^2
c       or      = 'GMCM' for moment of inertia units of gram*cm^2
c       or      = 'CM-1' for rotational constant in units of cm^-1
c       or      = 'MHZ' for rotational constant in units of MHz
c       or      = 'GHZ' for rotational constant in units of GHz
c         
         DO ki = 1 , 2
           IF     ( keyword(ki) .EQ. 'HAR' ) THEN
             VHAR = 'HAR'
           ELSEIF ( keyword(ki) .EQ. 'OBS' .OR. 
     &              keyword(ki) .EQ. 'FUN' ) THEN
             VHAR = 'OBS'
           ELSE
             VROT = keyword(ki)
             dum = 1.0
             CALL rotunits( VROT , dum )   ! check validity of keyword VROT: if not valid, rotunits will kill execution
           ENDIF
         END DO
c
c               Electronic level energies and degeneracies (for Nele levels); note that
c               lowest electronic energy level should be at Elev = 0 cm-1.
c
            Nflag = 0
         IF ( REPROD(i).EQ.REAC ) THEN
            a1=a1+1
            mass(a1)=AMU(i)
            Nreac = Nreac + 1
            Hdiff = Hdiff - DIV*DelH(i)       ! DIV=1000 for Hdiff expressed in cal/mole or J/mole; DIV=1 for cm-1
            Nflag = 1
         ELSEIF ( REPROD(i).EQ.PROD ) THEN
            Nprod = Nprod + 1
            Hdiff = Hdiff + DIV*DelH(i)       ! DIV=1000 for Hdiff in cal/mole or J/mole; DIV=1 for cm-1
            Nflag = 1
         ELSEIF ( REPROD(i).EQ.CTST) THEN
            Nctst = Nctst + 1
            Hdiff = Hdiff + DIV*DelH(i)       ! DIV=1000 for Hdiff in cal/mole or J/mole; DIV=1 for cm-1
            Nflag = 1
         ELSEIF ( REPROD(i).EQ.NON ) THEN
            Nnone = Nnone + 1
         ELSE 
            WRITE(3,*) '***FATAL: unrecognized reac/prod type: ',
     &          REPROD(i)
            write(*,*) '***FATAL: unrecognized reac/prod type: ',
     &          REPROD(i)
            stop 
         ENDIF

      IF ( Nctst .GT. 1 ) THEN
         WRITE(3,*) 'Too many transition states; only is one allowed!'
         write(*,*) 'Too many transition states; only is one allowed!'
         STOP 'Too many transition states; only is one allowed!'
      ENDIF

      IF ( Nctst .EQ. 1 .AND. Nprod.GT.0 ) THEN
         WRITE(3,*) '"prod" not allowed in Canonical TST calculation!'
         write(*,*) '"prod" not allowed in Canonical TST calculation!'
         STOP '"prod" not allowed in Canonical TST calculation!'
      ENDIF

c
c     Tunneling parameters: define and set units for forward and reverse barrier heights
c         Vf = forward barrrier in cm-1
c         Vr = reverse barrrier in cm-1
c         Hdiff = rxn enthalpy difference = Vf CTST calculation
c     Convert to units of cm-1
c
      IF ( Eunits .EQ. 'KCAL' ) THEN
        Vf = 349.76d-03*Hdiff                    ! Hdiff expressed in cal (not kcal)
        Vr = 349.76d+00*VVr                      ! convert reverse barrier to cm-1 units
      ELSEIF ( Eunits .EQ. 'KJOU' ) THEN
        Vf = 83.5946d-03*Hdiff                   ! Hdiff expressed in J (not kJ)
        Vr = 83.5946d+00*VVr                     ! convert reverse barrier to cm-1 units
      ELSEIF ( Eunits .EQ. 'CM-1' ) THEN
        Vf = Hdiff                               ! Hdiff expressed in cm-1
        Vr = VVr                                 ! reverse barrier in cm-1 units
      ENDIF

         zzpe(i) = 0.0
         ntop = 0             ! for counting symmetric tops (only one is allowed)
         IF ( N(i).GT.0 ) THEN
            DO II = 1 , N(i)
              READ (4,*) MODE(i,II) , IDOF(i,II)
              CALL ucase ( IDOF(i,II) )   ! convert to upper case

              IF ( IDOF(i,II) .EQ. RSO ) THEN  		! SPIN-ORBIT + ROTATION TREATMENT
                nso(i) = 1
                IF ( Nele(i) .NE. 0 ) then
                  WRITE(3,*)'*** FATAL ***'
                  WRITE(3,*)'  RSO requires Nele = 0 (Line 8 in input)'
                  WRITE(*,*)'*** FATAL ***'
                  WRITE(*,*)'  RSO requires Nele = 0 (Line 8 in input)'
                  STOP 
                ENDIF
                BACKSPACE (UNIT=4)
                READ (4,*) MODE(i,II) , IDOF(i,II) ,  
     &            MULT(i), LAMBDA(i), ASO(i), BSO(i), D(i), NRSO(i)
                CALL ucase ( IDOF(i,II) )   ! convert to upper case

              ELSEIF ( IDOF(i,II).EQ.VIB ) THEN                     ! VIBRATION
               BACKSPACE (UNIT=4)
               READ (4,*) MODE(i,II) , IDOF(i,II) , WE(i,II) , ANH(i,II)
     &                    , NG(i,II)
               CALL ucase ( IDOF(i,II) )   ! convert to upper case
                if(NG(i,II).eq.10.AND.REPROD(i).EQ.PROD) then
                   freq=WE(i,II)
                   NG(i,II)=1
                endif
                  W(i,II) = WE(i,II)
                  IF ( VHAR .EQ. HAR ) THEN
                   WOBS(i,II) = WE(i,II)+2.*ANH(i,II)           ! WOBS = observed freq
                  ELSEIF ( VHAR .EQ. OBS ) THEN
                   WOBS(i,II) = WE(i,II)
                   WE(i,II)= WE(i,II)-2.*ANH(i,II)        ! Convert WE to harmonic frequency
                  ENDIF
                  CALL qmorse(300.d+00,WE(i,II),ANH(i,II),NG(i,II),
     &                    qq,Cvib,Svib,Hvib,zap)
                  zzpe(i) = zzpe(i) + zap

              ELSEIF ( IDOF(i,II).EQ.BOX ) THEN                ! PARTICLE IN A BOX 
               BACKSPACE (UNIT=4)
               READ (4,*) MODE(i,II) , IDOF(i,II) , WE(i,II) , ANH(i,II)
     &                    , NG(i,II)
               CALL ucase ( IDOF(i,II) )   ! convert to upper case
                  W(i,II) = WE(i,II)
                  IF ( VHAR .EQ. HAR ) THEN
                   WOBS(i,II) = 3.0d+00*WE(i,II)                ! WOBS = observed freq (note that ANH is a placeholder: not used)
                  ELSEIF ( VHAR .EQ. OBS ) THEN
                   WOBS(i,II) = WE(i,II)
                   WE(i,II)= WE(i,II)/3.0d+00             ! Convert WE to harmonic frequency
                  ENDIF
                 zap = W(i,II)
                 zzpe(i) = zzpe(i) + zap

              ELSEIF ( IDOF(i,II).EQ.ROT ) THEN                ! CLASSICAL ROTOR
               BACKSPACE (UNIT=4)
               READ (4,*) MODE(i,II) , IDOF(i,II) , WE(i,II) , ANH(i,II)
     &                    , NG(i,II)
               CALL ucase ( IDOF(i,II) )   ! convert to upper case
     
       IF ( NG(i,II) .GT. 2 ) THEN
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
                  call rotunits( VROT , WE(i,II) )                ! unit conversion
                  AMOM(i,II) = WE(i,II)

              ELSEIF ( IDOF(i,II).EQ.QRO ) THEN                ! QUANTUM ROTOR
               BACKSPACE (UNIT=4)
               READ (4,*) MODE(i,II) , IDOF(i,II) , WE(i,II) , ANH(i,II)
     &                    , NG(i,II)
               CALL ucase ( IDOF(i,II) )   ! convert to upper case
                  call rotunits( VROT , WE(i,II) )                ! unit conversion
                  AMOM(i,II) = WE(i,II)
             
              ELSEIF ( IDOF(i,II).EQ.TOP ) THEN               ! Symmetric Top
               BACKSPACE (UNIT=4)
               READ (4,*) MODE(i,II) , IDOF(i,II) , WE(i,II) , ANH(i,II)
     &                    , NG(i,II)
               CALL ucase ( IDOF(i,II) )   ! convert to upper case
                    ntop = ntop + 1     ! for counting symmetric tops
                    IF ( ntop .GT. 1) THEN
                    write (*,*) '*** FATAL: only one "top" allowed ***'
                    write (3,*) '*** FATAL: only one "top" allowed ***'
                    STOP
                    ENDIF
                  call rotunits( VROT , WE(i,II) )                ! unit conversion BJ
                  AMOM(i,II) = WE(i,II)
                  call rotunits( VROT , ANH(i,II) )               ! unit conversion BK

              ELSEIF (IDOF(i,II).EQ. GOR) THEN                 ! GORIN
               BACKSPACE (UNIT=4)
               READ (4,*) MODE(i,II) , IDOF(i,II) , WE(i,II) , ANH(i,II)
     &                    , NG(i,II)
               CALL ucase ( IDOF(i,II) )   ! convert to upper case
                  call rotunits( VROT , WE(i,II) )                ! unit conversion
                  AMOM(i,II) = WE(i,II)
                  Gorin=Gorin+1
              ELSEIF ( IDOF(i,II).EQ.FIT) THEN                ! HINDRANCE
               BACKSPACE (UNIT=4)
               READ (4,*) MODE(i,II) , IDOF(i,II) , WE(i,II) , ANH(i,II)
     &                    , NG(i,II)
               CALL ucase ( IDOF(i,II) )   ! convert to upper case
                  call rotunits( VROT , WE(i,II) )                ! unit conversion
                  AMOM(i,II) = WE(i,II)
                  hindrance=hindrance+1
                  Weorig(hindrance)= WE(i,II)

              ELSEIF ( IDOF(i,II).EQ.HRA ) THEN                ! SYMMETRICAL HINDERED ROTOR 
               BACKSPACE (UNIT=4)
               READ (4,*) MODE(i,II) , IDOF(i,II) , WE(i,II) , ANH(i,II)
     &                    , NG(i,II)
               CALL ucase ( IDOF(i,II) )   ! convert to upper case
                  IF ( NG(i,II) .LT. 0 ) THEN
                      READ(4,*) NSIG(i,II)
                      NG(i,II) = -NG(i,II)
                    ELSE
                      NSIG(i,II) = NG(i,II)
                  ENDIF
                  W(i,II) = WE(i,II)                    ! small vibration frequency
                  AMOM(i,II) = ANH(i,II)                ! moment of inertia
                  call rotunits( VROT , AMOM(i,II) )                ! unit conversion to AMU*Ang^2
                  B(i,II) = 16.85763D+00/AMOM(i,II)     ! rotational constant (cm-1)
                  VV(i,II) = ((W(i,II)/NG(i,II))**2)/B(i,II)
                  NSV(i,II)=NSIG(i,II)
c                 to find ZPE:
                  CALL SHRLEV(Emax,EVh,B(i,II),VV(i,II),
     &            NSIG(i,II),1,IMAX(i,II),zap,'Vhrd1','Bhrd1',
     &            0.0d0,0.0d0,Vo,Vm)
                  zpe(i,II) = zap
                  Vmax(i,II) = Vo
                  Vmin(i,II) = Vm 
                  DO LL=1, IMAX(i,II)
                    EV(i,II,LL)=EVh(LL)
                  ENDDO
                  zzpe(i) = zzpe(i) + zap
                  IDOF(i,II) = HIN

              ELSEIF ( IDOF(i,II).EQ.HRB ) THEN                ! SYMMETRICAL HINDERED ROTOR 
               BACKSPACE (UNIT=4)
               READ (4,*) MODE(i,II) , IDOF(i,II) , WE(i,II) , ANH(i,II)
     &                    , NG(i,II)
               CALL ucase ( IDOF(i,II) )   ! convert to upper case
                  IF ( NG(i,II) .LT. 0 ) THEN
                      READ(4,*) NSIG(i,II)
                      NG(i,II) = -NG(i,II)
                    ELSE
                      NSIG(i,II) = NG(i,II)
                  ENDIF
                  W(i,II) = WE(i,II)                    ! small vibration frequency
                  VV(i,II) = ANH(i,II)                  ! hindrance barrier (cm-1)
                  IF (VV(i,II) .LT. 1.0d-05) STOP
     &            'Rotation barrier too small. Hit RETURN to terminate'
                  B(i,II) = ((W(i,II)/NG(i,II))**2)/VV(i,II)		! rotational constant (cm-1)
                  AMOM(i,II) = 16.85763D+00/B(i,II)     			! rotational constant (cm-1)
                  NSV(i,II)=NSIG(i,II)
                  CALL SHRLEV(Emax,EVh,B(i,II),VV(i,II),    ! To find zero-point energy
     &            NSIG(i,II),1,IMAX(i,II),zap,'Vhrd1','Bhrd1',
     &            0.0d0,0.0d0,Vo,Vm)
                  zpe(i,II) = zap
                  Vmax(i,II) = Vo
                  Vmin(i,II) = Vm
                  DO LL=1, IMAX(i,II)
                        EV(i,II,LL)=EVh(LL)
                  ENDDO
                  zzpe(i) = zzpe(i) + zap
                  IDOF(i,II) = HIN

              ELSEIF ( IDOF(i,II).EQ.HRC ) THEN                ! SYMMETRICAL HINDERED ROTOR 
               BACKSPACE (UNIT=4)
               READ (4,*) MODE(i,II) , IDOF(i,II) , WE(i,II) , ANH(i,II)
     &                    , NG(i,II)
               CALL ucase ( IDOF(i,II) )   ! convert to upper case
                  IF ( NG(i,II) .LT. 0) THEN
                      READ(4,*) NSIG(i,II)
                      NG(i,II) = -NG(i,II)
                    ELSE
                      NSIG(i,II) = NG(i,II)
                  ENDIF
                  AMOM(i,II) = WE(i,II)                 ! moment of inertia
                  call rotunits( VROT , AMOM(i,II) )                ! unit conversion to AMU*Ang^2
                  B(i,II) = 16.85763D+00/AMOM(i,II)     ! rotational constant (cm-1)
                  VV(i,II) = ANH(i,II)                  ! hindrance barrier (cm-1)
                  W(i,II) = NG(i,II)*SQRT(B(i,II)*VV(i,II))
                                                        ! small vibration frequency
                  NSV(i,II)=NSIG(i,II)
                  CALL SHRLEV(Emax,EVh,B(i,II),VV(i,II),    ! To find zero-point energy
     &            NSIG(i,II),1,IMAX(i,II),zap,'Vhrd1','Bhrd1',
     &            0.0d0,0.0d0,Vo,Vm)
                  zpe(i,II) = zap
                  Vmax(i,II) = Vo
                  Vmin(i,II) = Vm
                  DO LL=1, IMAX(i,II)
                        EV(i,II,LL)=EVh(LL)
                  ENDDO
                  zzpe(i) = zzpe(i) + zap
                  IDOF(i,II) = HIN

              ELSEIF ( IDOF(i,II).EQ.HRD ) THEN       ! GENERAL, UNSYMMETRICAL HINDERED ROTOR
               BACKSPACE (UNIT=4)
               READ (4,*) MODE(i,II) , IDOF(i,II) , WE(i,II) , ANH(i,II)
     &                    , NG(i,II)
               CALL ucase ( IDOF(i,II) )   ! convert to upper case
                  NSIG(i,II) = NG(i,II)
                  NCV(i,II)=INT(WE(i,II))
                  NCB(i,II)=INT(ANH(i,II))
                  READ(4,*) VHRD(i,II), NSV(i,II), Phav(i,II), 
     &                        (CV(i,II,LL), LL=1, NCV(i,II))
                  READ(4,*) BHRD(i,II), NSB(i,II), Phab(i,II), 
     &                        (CB(i,II,LL), LL=1, NCB(i,II))
                  CALL ucase ( VHRD(i,II) )   ! convert to upper case
                  CALL ucase ( BHRD(i,II) )   ! convert to upper case

                  DO LL=1, NCB(i,II)
                        CBB(LL)=CB(i,II,LL)
                  ENDDO
                  DO LL=1, NCV(i,II)
                        CVV(LL)=CV(i,II,LL)
                  ENDDO 
                  CALL UHRLEV(Emax,EVh,NCB(i,II),NCV(i,II),    ! To find zero-point energy
     &            CBB,CVV,NSV(i,II),NSB(i,II),IMAX(i,II),zap,
     &            VHRD(i,II),BHRD(i,II),Phav(i,II),Phab(i,II),Vo,Vm)
                  zpe(i,II) = zap
                  Vmax(i,II) = Vo
                  Vmin(i,II) = Vm
                  DO LL=1, IMAX(i,II)
                        EV(i,II,LL)=EVh(LL)
                  ENDDO
                  zzpe(i) = zzpe(i) + zap
                  IDOF(i,II) = HRD
               
              ELSEIF ( IDOF(i,II) .EQ. QVB ) THEN      ! ANHARMONIC VIBRATIONAL PARTITION FUNCTION, external file
               BACKSPACE (UNIT=4)
               READ (4,*) MODE(i,II) , IDOF(i,II) , WE(i,II) , ANH(i,II)
     &                    , NG(i,II)
               CALL ucase ( IDOF(i,II) )   ! convert to upper case
               OPEN (12,FILE=MolName(i)(1:lenstr(MolName(i)))//'.qvib',
     &                STATUS='OLD')
               NQVB(i) = 1
               READ (12,99030) dumcut                 ! Check to see if file format includes input data summary
               IF ( dumcut .EQ. cut2 ) THEN            ! Format starting with v.2008.1
                  dumcut = 'x'
                  DO WHILE (dumcut .NE. cut2)
                    READ (12,99030) dumcut
                  END DO
               ELSE
                  WRITE(3,*) '****FATAL: wrong format for qvib file: ',
     &                     IDOF(i,II)
                  WRITE(*,*) '****FATAL: wrong format for qvib file: ',
     &                     IDOF(i,II)
                  STOP
                ENDIF

               READ(12,99031) TITLE4(i,II)
               READ(12,*) Egrain1(i,II),Emax2(i,II), zpp(i),KEYWORD2(i)
               zzpe(i) = zzpe(i) + zpp(i)
               READ (12,*) Ntx                       ! *********** number of temperatures to be passed to THERMO
               IF ( Ntx .LE. maxNtx ) THEN
                 ELSE
                  WRITE(3,*) '****FATAL: too many Temps in qvib file: ',
     &                     IDOF(i,II)
                  WRITE(*,*) '****FATAL: too many Temps in qvib file: ',
     &                     IDOF(i,II)
                  STOP
               END IF
               READ(12,99032) DUMMY
               
               DO LL = 1 , Ntx
                 READ(12,*) LK, Tx(LL) , QT(i,LL) , Cx(i,LL),
     &                     Hx(i,LL) ,Sx(i,LL)                   ! i=species, LL=temperature
               END DO
               CLOSE (12)

              ELSEIF ( IDOF(i,II) .EQ. CRP ) THEN      ! SEMI-CLASSICAL TST, external file
               BACKSPACE (UNIT=4)
               READ (4,*) MODE(i,II) , IDOF(i,II) , WE(i,II) , ANH(i,II)
     &                    , NG(i,II)
               CALL ucase ( IDOF(i,II) )   ! convert to upper case
               IF (REPROD(i) .EQ. CTST) THEN
                  VVR = 0.0d+00
                  IF ( vimag .GE. 1.0 ) THEN
                     WRITE(*,*)
     &               '***FATAL: CRP type requires vimag < 1.0 ',
     &               'on CTST specification line (to prevent a second ',
     &               'tunneling correction); also, VVR is set to 0.0***'
                     WRITE(3,*)
     &               '***FATAL: CRP type requires vimag < 1.0 ',
     &               'on CTST specification line (to prevent a second ',
     &               'tunneling correction); also, VVR is set to 0.0***'
                     STOP
                  ENDIF
                ELSEIF (REPROD(i) .NE. CTST) THEN
                  WRITE(*,*) 
     &            '***FATAL: CRP type can only be used with CTST***'
                  WRITE(3,*) 
     &            '***FATAL: CRP type can only be used with CTST***'
                  STOP
               ENDIF
               NCRP = 1
               SCTST='YES'                            ! flag to control printing (or not) of Eckart tunneling data
               OPEN (12,FILE=MolName(i)(1:lenstr(MolName(i)))//'.qcrp',
     &                STATUS='OLD')
               READ (12,99030) dumcut                 ! Check to see if file format includes input data summary
               IF ( dumcut .EQ. cut2 ) THEN            ! Format starting with v.2008.1
                  dumcut = 'x'
                  DO WHILE (dumcut .NE. cut2)
                    READ (12,99030) dumcut
                  END DO
               ELSE
                  WRITE(3,*) '****FATAL: wrong format for qcrp file: ',
     &                     IDOF(i,II)
                  WRITE(*,*) '****FATAL: wrong format for qcrp file: ',
     &                     IDOF(i,II)
                  STOP
                ENDIF
               READ(12,99031) TITLE5(i,II)
               READ(12,*) Egrain1(i,II), Emax2(i,II), Vf, Vr, 
     &                     zpp(i), KEYWORD2(i) , VPTx(i)
               Emax1(i,II) = (Imax1-1.0)*Egrain1(i,II)
               zzpe(i) = zzpe(i) + zpp(i)
               READ (12,*) Ntx                       ! *********** number of temperatures to be passed to THERMO by SCTST
               IF ( Ntx .LE. maxNtx ) THEN
                 ELSE
                  WRITE(3,*) '****FATAL: too many Temps in qcrp file: ',
     &                     IDOF(i,II)
                  WRITE(*,*) '****FATAL: too many Temps in qcrp file: ',
     &                     IDOF(i,II)
                  STOP
               END IF
               READ(12,99032) DUMMY
               DO LL = 1 , Ntx                        ! Number of temperatures generated by ADENSUM and SCTST
                 READ(12,*) LK , Tx(LL) , QT(i,LL) , Cx(i,LL),
     &                     Hx(i,LL) ,Sx(i,LL)                   ! i=species, LL=temperature
               END DO
               CLOSE (12)

99030 FORMAT (A33)
99031 FORMAT(/,A150)
99032 FORMAT(/////,A150)

              ELSEIF ( IDOF(i,II) .NE. ECK ) THEN
                  WRITE(3,*) "**** FATAL: DOF type not recognized: ",
     &                   IDOF(i,II)        
                  write(*,*) "**** FATAL: DOF type not recognized: ",
     &                   IDOF(i,II)        
                  stop               
              ENDIF
            END DO !  II do loop
      ENDIF
 
      END DO   ! i species
  
      IF(hindrance.gt.0) READ (4,*) (Kexp(j),j=1,Nt)                !  Gorin 
      IF(Gorin.gt.0) then                                           !  Gorin 
         read(4,*) readpot   ! Potential
         read(4,*) freq      ! CM-1
         read(4,*) De        ! Eunit    Without ZPE
         read(4,*) re        ! Re, between centers of mass, A
         Deorig=De
         IF ( Eunits.EQ.'KCAL' ) THEN
          De=De*349.76
          ELSEIF ( Eunits.EQ.'KJOU' ) THEN
          De=De*83.59
         ENDIF

         a1=0
         do i=1,Ns
           IF(REPROD(i).EQ.REAC) then
            a1=a1+1
            mass(a1)=AMU(i)
           ENDIF 
         enddo 
         IF(a1.ne.2) then
           WRITE(3,*) "**** FATAL: Wrong number of reactants ****"
           write(*,*) "**** FATAL: Wrong number of reactants ****"
           stop
         ENDIF

         mu=(mass(1)*mass(2))/(mass(1)+mass(2))   
         beta=0.1217788*freq*SQRT(mu/De)
        
         readpot=readpot(1:LEN_TRIM(readpot))
         IF(readpot.eq."MORSE".OR.readpot.eq."morse") then
          POT="MORSE"
         elseif(readpot.eq."VARSHNI".OR.readpot.eq."varshni") then
          POT="VARSHNI"
          beta=0.5*(beta-(1/re))/re
         elseif (readpot.eq."sMORSE".OR.readpot.eq."smorse") then
          POT="sMORSE"
          read(4,*) c        ! read stiff Morse Parameter
         ELSE
          WRITE(3,*) "**** FATAL: Potential not recognized"
          write(*,*) "**** FATAL: Potential not recognized"
          stop
         ENDIF

      ENDIF  ! Gorin
        
      CLOSE(4)

      RETURN 
      end subroutine 
