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

      MODULE I_O_mod
      
!     The following procedures are contained in this module:
!
!     SUBROUTINE write_1      Copyright 2001 - 2021 John R. Barker
!     SUBROUTINE write_2      Copyright 2001 - 2021 John R. Barker
!     SUBROUTINE rxn_input    Copyright 2001 - 2021 John R. Barker

      IMPLICIT NONE
      SAVE

      CHARACTER(len=7)  :: AVERSION = '2023.1'
      CHARACTER(len=8)  :: ADATE = 'Apr 2023'
 
      REAL(8) :: Temp
      REAL(8) :: Tvib
      REAL(8) :: P(50) , PP(50) , Time , E , XRANST , AA , EE ,
     &        Einit , RR , CollFreq , EM, Rob_1 , SCALE , kcol

      CHARACTER(len=3) :: Punits
      CHARACTER(len=3) :: bar = 'BAR'
      CHARACTER(len=3) :: atm = 'ATM'
      CHARACTER(len=3) :: tor = 'TOR'
      CHARACTER(len=3) :: mcc = 'MCC'
      CHARACTER(len=4) :: POUT
      CHARACTER(len=4) :: Eunits
      CHARACTER(len=4) :: kcal = 'KCAL'
      CHARACTER(len=4) :: kjou = 'KJOU'
      CHARACTER(len=4) :: cmm1 = 'CM-1'
      CHARACTER(len=4) :: Rotatunits
      CHARACTER(len=4), DIMENSION(3) :: DUM 
      CHARACTER (len=10) :: dummy
      CHARACTER(len=8) :: KEYTEMP 
      CHARACTER(len=8), DIMENSION(5) :: KEYWORD
      CHARACTER(len=100) :: TITLE
      CHARACTER(len=10)  :: TS
      INTEGER :: II, III, Mdum, Mol, Icalc, MolInit, Idum0
      INTEGER :: IR 
      REAL(8) :: Rtrials


      CONTAINS

!---------------------------------------------------------------------------------------
      SUBROUTINE write_1
      
      USE declare_mod
      USE utili_mod
      USE subs_mod
      IMPLICIT NONE
      SAVE
      INTEGER:: i, ki
      
      DO ki = 7 , 11
        WRITE (ki,99001) AVERSION , ADATE , AVERSION , ADATE
        CALL DateTime(ki)
        WRITE (ki,99003) TITLE
        WRITE (ki,99044) Egrain1 , Emax1 , imax1 , Egrain2 , Emax2 , 
     &                   (Isize-imax1)
      END DO
      DO ki = 7 , 9
        WRITE (ki,99031) Temp , Tvib , NWells , NProds
        WRITE (ki,99036) (i,i=1,NWells+NProds)
        WRITE (ki,99035) (ADJUSTR(MolName(i)),i=1,NWells+NProds)
        WRITE (ki,99023) (HMol(i),i=1,NWells+NProds)
        WRITE (ki,99025) (HMol(i)/jtocm,i=1,NWells+NProds)
        WRITE (ki,99024) (HMol(i)/caltocm,i=1,NWells+NProds)
        WRITE (ki,99033) (Viblo(i),i=1,NWells)
        WRITE (ki,99045) (MolMom(i),i=1,NWells)
        WRITE (ki,99046) (Molsym(i),i=1,NWells)
        WRITE (ki,99047) (Molele(i),i=1,NWells)
        WRITE (ki,99048) (Molopt(i),i=1,NWells)
        WRITE (ki,99041) (Etherm(i,Temp),i=1,NWells)
        WRITE (ki,99009)
      END DO ! ki

      WRITE (KARY,99010) (ADJUSTR(MolName(i)),i=1,NWells)

      DO II = 1 , Isize
         IF ( II.LE.imax1 ) THEN
            E = (II-1)*Egrain1
         ELSE
            E = (II-imax1-1)*Egrain2
         ENDIF
         WRITE (KARY,99011) E , (exp(Dens(i,II)),i=1,NWells)
                                                        ! Densities stored as logs
       END DO ! II

      DO ki = 7 , 9
        WRITE(ki,99073) NWells
        IF ( ETKEY.EQ.newet ) THEN
          WRITE(ki,*) ETKEY
          WRITE(ki,*) 'BARKER''S ''NEW APPROACH'' TO ENERGY TRANSFER'
          WRITE(ki,*) 
     &       '[J.R. Barker, Int. J. Chem. Kinetics, 41, 748-763 (2009)]'
        ELSEIF ( ETKEY.EQ.oldet ) THEN
          WRITE(ki,*) ETKEY
          WRITE(ki,*) 'TRADITIONAL APPROACH TO ENERGY TRANSFER'
        ELSE
          WRITE(ki,*) ETKEY
          WRITE(ki,*) 'UNSPECIFIED ET TREATMENT: USE BARKER''S APPROACH'
          WRITE(ki,*) 
     &       '[J.R. Barker, Int. J. Chem. Kinetics, 41, 748-763 (2009)]'
        ENDIF

        WRITE (ki,99055)
        Mdum = 0
        WRITE (ki,99056) Mdum , SigM , EpsM , AmuM
      END DO

      DO Mol = 1, Nwells
        DO ki = 7 , 9
         WRITE (ki,99056) Mol , Sig(Mol) , Eps(Mol) , Amu , 
     &                    Sigma(Mol) , Epsil(Mol) ,  Mass(Mol) , 
     &                     klj(Mol) , kqm(Mol) , CollK(Mol)
        END DO
      END DO
c
c     Call PSTEP WITH E=-1.0 in order to store LI(Mol,i) in COMMON for print-out
c
c      DO ki = 7, 9
c        WRITE(ki,*) '(end of subroutine write_1)'
c      END DO       

! -----------------------------------------------------------------------       
! -----------------------------------------------------------------------       
!     FORMAT STATEMENTS
! -----------------------------------------------------------------------       

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
99003 FORMAT (A100)
99009 FORMAT (/'Units: 2D-moment : amu Ang^2',
     &         '       Thermal<E>: cm-1 (at Translational T)')
99010 FORMAT (/,5X,'Densities of States',/'     Energy',1x,150(1x,A10))
99011 FORMAT (1x,f10.1,150(1x,1pe10.3,1x,1pe10.3))
99023 FORMAT (' delH(cm-1):',150(1x,f10.1))
99024 FORMAT (' delH(kcal):',150(1x,f10.1))
99025 FORMAT ('   delH(kJ):',150(1x,f10.1))
99031 FORMAT (/,'TEMPERATURES:',/,
     &          '        Translational Temperature =',
     &        f8.1,1x,'K'/,    
     &          '  Initial Vibrational Temperature =',
     &        f8.1,1x,'K'///,I3,' WELLS and',I3,' PRODUCTS'/)
99033 FORMAT (' lowest_vib:',150(1x,f10.1))
99035 FORMAT (' Molec_Name:',150(1x,A10))
99036 FORMAT ('  Molec_No.:',150(1x,I10))
99041 FORMAT (' Thermal<E>:',150(1x,f10.1))
99044 FORMAT (/'PARAMETERS FOR DOUBLE-ARRAYS:'/
     &    '        Grain (cm-1)      Maximum Energy    Number of grains'
     &    /2(2F20.3,i20/))
99045 FORMAT ('  2D-moment:',150(1x,f10.1))
99046 FORMAT ('   ext_symm:',150(1x,I10))
99047 FORMAT ('Qel_partfxn:',150(1x,F10.2))
99048 FORMAT ('   opt_isom:',150(1x,I10))
99055 FORMAT (/,'Lennard-Jones Parameters & Pdown Models    ',
     &        '[rate constants k(LJ) and k(QM) units: cc/sec]',//,
     &'Mol  Sig(Angstrom)  Eps/k MolWt(amu)     NetSig   NetEps/k    Red
     &Mass     k(LJ)    k(QM)     Used')
99056 FORMAT (1x,I2,6(1X,f10.3),3(1x,1pe10.3))
99073 FORMAT (//,'COLLISION PARAMETERS FOR ',I2,' WELLS')

      END SUBROUTINE write_1

!---------------------------------------------------------------------------------------
      SUBROUTINE write_2
      
      USE declare_mod
      USE utili_mod
      USE subs_mod
      IMPLICIT NONE
      SAVE
      INTEGER:: i, j, k, ki
      
      REAL(8) :: xref
      REAL(8) :: kcol
      REAL(8) :: kx, ax, ae, k0
      REAL(8), DIMENSION(MaxWell,MaxChan) :: Kinf
      REAL(8), DIMENSION(MaxWell,MaxChan) :: Ainf
      REAL(8), DIMENSION(MaxWell,MaxChan) :: Einf
      REAL(8), DIMENSION(MaxWell,MaxChan) :: kosc 
      
      
      DO ki = 7, 9
        WRITE(ki,*) ' '
        WRITE(ki,*) 'ENERGY TRANSFER STEP-SIZE PARAMETERS'
      END DO

      DO 600 Mol = 1 , NWells
         Rob_1 = -10.0
         CALL PSTEP(Rob_1,Rob_1,Mol,Temp,XRANST,SCALE)            ! call to get LI(Mol) for writing output
         xref = Eref(Mol)
         IF ( ETKEY.Eq.oldet ) xref = 0.0
         DO ki = 7 , 9
           WRITE (ki,99057) Mol , ADJUSTR(MolName(Mol)),
     &      xref, (DC(Mol,i),i=1,8)
           WRITE (ki,99027) (LI(Mol,i),i=1,5)
         END DO
 600  CONTINUE
 
      WRITE (KARY,99030)
      WRITE (KARY,99012) (Mol,Mol,Mol=1,NWells)
      DO 700 II = 1 , Isize
         IF ( II.LE.imax1 ) THEN
            E = (II-1)*Egrain1
         ELSE
            E = (II-imax1-1)*Egrain2
         ENDIF
         WRITE (KARY,99011) E , 
     &        (CNORM(Mol,II),COLLUP(Mol,II),Mol=1,NWells)
 700  CONTINUE
 
c------------------------------------------------------------------------------------------
c   WRITE OUT RATE CONSTANT INFO TABLE
c------------------------------------------------------------------------------------------
!     COMPUTE THERMAL KINF, ETC.

      DO 900 j = 1 , NWells
         IF ( Nchan(j) .NE. 0 ) THEN
           DO 850 i = 1 , Nchan(j)            
             if(itun(j,i).eq.1)then
               CALL QKinftun(j,i,Temp,kx,ax,ae)    ! high-pressure rate constant w/ tunneling effects
               Kinf(j,i) = kx
               Ainf(j,i) = ax
               Einf(j,i) = ae
             else
               kcol = CollK(j)                                      ! collision rate constant previously selected (kLJ or kQM)
               CALL QKinf(j,i,Temp,kx,ax,ae,kcol,k0)
               Kinf(j,i) = kx
               Ainf(j,i) = ax
               Einf(j,i) = ae
               kosc(j,i) = k0
             end if
 850       CONTINUE
         END IF
 900  CONTINUE

      IF ( Nreac .GT. 0 ) THEN
      DO ki = 7 , 11                               ! Write to Units 7, 8, 9, 10, and 11
       WRITE (ki,99013)  Nreac
       WRITE (ki,99028) ((j,Jto(j,i),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99034) ((ADJUSTR(TSname(j,i)),i=1,Nchan(j)),
     &   j=1,NWells)
       WRITE (ki,99045) ((TSmom(j,i),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99046) ((TSsym(j,i),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99047) ((TSele(j,i),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99048) ((TSopt(j,i),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99026) ((Path(j,i),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99023) ((Hts(j,Jto(j,i)),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99025) ((Hts(j,Jto(j,i))/jtocm,i=1,Nchan(j)),j=1,
     &                   NWells)
       WRITE (ki,99024) ((Hts(j,Jto(j,i))/caltocm,i=1,Nchan(j)),j=1,
     &                   NWells)
       WRITE (ki,99021) ((Eo(j,i),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99022) ((Eo(j,i)/jtocm,i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99020) ((Eo(j,i)/caltocm,i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99029) ((Kinf(j,i),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99049) ((Ainf(j,i),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99050) ((Einf(j,i),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99051) ((Einf(j,i)*caltoj,i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99086) ((kosc(j,i),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99115) ((NCENT(j,i),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99114) ((Eor(j,i),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99052) ((Jread(j,i),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99088) ((vimag(j,i),i=1,Nchan(j)),j=1,NWells)
       WRITE (ki,99090) ((WARN1(j,i),i=1,Nchan(j)),j=1,NWells)
       IF ( nivrflag .EQ. 1 ) THEN
c        WRITE (ki,99080) ((iivr(j,i),i=1,Nchan(j)),j=1,NWells)
         WRITE (ki,99081) ((vivr(j,i),i=1,Nchan(j)),j=1,NWells)
         WRITE (ki,99087) ((vave(j,i),i=1,Nchan(j)),j=1,NWells)
         WRITE (ki,99082) ((pivr(j,i),i=1,Nchan(j)),j=1,NWells)
         WRITE (ki,99186) ((tivr(j,i),i=1,Nchan(j)),j=1,NWells)
         WRITE (ki,99083) ((civr(j,i,1),i=1,Nchan(j)),j=1,NWells)
         WRITE (ki,99084) ((civr(j,i,2),i=1,Nchan(j)),j=1,NWells)
         WRITE (ki,99085) ((civr(j,i,3),i=1,Nchan(j)),j=1,NWells)
       ENDIF
       IF ( Pressflag .GT. 0 ) THEN
         WRITE (ki,99089) ((Press(j,i),i=1,Nchan(j)),j=1,NWells)
       ENDIF
       WRITE (ki,99053)

99153  FORMAT(//,
     & 'SUPPLEMENTARY Rate Constants (unimol, or bimol w. Bath Gas)',/,
     &  'From    To    TSname  rxn_order  Afac(cm^3 s-1)  Ea/R  ',
     &  'ksup(T:', F5.0,  ')  QuenchE(cm-1)')
99154  FORMAT(I3, 3x, I3, 3x, A10, 2x, I2, 7x, 1PE10.3 ,2x, 0PF8.2, 1x, 
     &        1PE10.3, 2x, 0pF8.0 )
       IF ( nsux .GT. 0 ) THEN                                          ! if there are supplementary reactions
         WRITE (ki, 99153) Temp
         DO Mol = 1, NWells
           IF ( nsmax(Mol) .GT. 0 ) THEN
            DO i = 1, nsmax(Mol)
            WRITE (ki, 99154) Mol, sto( Mol,i ), SupName(Mol,i) ,  
     &       norder(Mol ,i), Afc(Mol, i), Bsr(Mol, i) , ksup( Mol,i ) , 
     &       dde( Mol,i )
            END DO   ! i
           END IF
         END DO   ! Mol
       END IF

      END DO  ! ki,  Write to Units 7, 8, 9, 10, and 11

c------------------------------------------------------------------------------------------
c   WRITE OUT RATE CONSTANTS IN MULTIWELL.ARRAY OUTPUT FILE
c------------------------------------------------------------------------------------------
      WRITE(KARY,*) '  '
      WRITE (KARY,99034) ((ADJUSTR(TSname(j,i)),i=1,Nchan(j)),
     &   j=1,NWells)
      WRITE (KARY,99114) ((Eor(j,i),i=1,Nchan(j)),j=1,NWells)
      WRITE (KARY,99028) ((j,Jto(j,i),i=1,Nchan(j)),j=1,NWells)

!     combine tunneling and nontunneling k(E) into KRate for output purposes

        do i=1,NWells
          do k=1,Nchan(i)
            do II=1,Isize+imax1
              if(II.le.imax1*2)then
                E=(II-1)*Egrain1-(imax1*Egrain1)
              else
                E=(II-imax1*2-1)*Egrain2
              end if
              if(itun(i,k).eq.1)then
                KRate(i,k,II)=TRate(i,k,II)
              else
                if(E.lt.0.0d0)then
                  KRate(i,k,II)= -5000.0d0     !-5000.0d0   ensures log[KRate(i,k,II)] is 0 below the rxn barrier  ???
                else
                  III=II-imax1                !III index ensures Rate array for the nontunneling case
                  KRate(i,k,II)=Rate(i,k,III) !is stored correctly in the KRate array
                end if
              end if
            end do
          end do
        end do

!     output all k(E)

      IF ( ntunflag .EQ. 1 ) THEN                ! IF TUNNELING, write k(E) from below critical energy Eo

        do II=1,Isize+imax1
          if(II.le.imax1*2)then
            E=(II-1)*Egrain1-(imax1*Egrain1)
          else
            E=(II-imax1*2-1)*Egrain2
          end if
          WRITE(KARY,99018)E,
     &         ((EXP(KRate(j,i,II)),i=1,Nchan(j)),j=1,NWells)
        end do

      ELSE                                       ! IF NO TUNNELING, write k(E) starting at critical energy Eo

        do II= imax1+1 , Isize+imax1
          if(II.le.imax1*2)then
            E=(II-1)*Egrain1-(imax1*Egrain1)
          else
            E=(II-imax1*2-1)*Egrain2
          end if
          WRITE(KARY,99018)E,
     &         ((EXP(KRate(j,i,II)),i=1,Nchan(j)),j=1,NWells)
        end do

      ENDIF
      
      ENDIF     ! End printing when Nreac .GT. 0

         WRITE (KSMM,99059) Rtrials
         IF ( Icalc.EQ.1 ) THEN
          WRITE (KSMM,99165) Temp , ADJUSTR(MolName(Molinit)) , 
     &                         Einit , IDUM0 , IDUM
         ELSEIF ( Icalc.EQ.2 ) THEN
          WRITE (KSMM,99166) Tvib , Temp , ADJUSTR(MolName(Molinit)) , 
     &                         Einit , IDUM0 , IDUM
         ELSEIF ( Icalc.EQ.3 ) THEN
          WRITE (KSMM,99167) Tvib , Temp , ADJUSTR(MolName(Molinit)) ,
     &                         MolName(IR) , Einit, IDUM0 , IDUM
         ELSEIF ( Icalc.EQ.4 ) THEN
          WRITE (KSMM,99168) Temp , ADJUSTR(MolName(Molinit)) ,
     &                         IDUM0 , IDUM
         ENDIF
 
         WRITE (KSMM,99063) (ADJUSTR(MolName(i)),i=1,NWells+NProds)
         WRITE (KSMM,99060) (HMol(i),i=1,NWells)
         WRITE (KSMM,99161) (i,i,i=1,NWells+NProds)

c      DO ki = 7, 9
c        WRITE(ki,*) '(end of subroutine write_2)'
c      END DO       
      
! -----------------------------------------------------------------------       
!     FORMAT STATEMENTS
! -----------------------------------------------------------------------       

99011 FORMAT (1x,f10.1,150(1x,1pe10.3,1x,1pe10.3))
99012 FORMAT (/'     Energy',150(:'   Norm(',I2.2,') UpProb(',I2.2,')'))
99013 FORMAT (/,5X,I3,' REACTION RATE CONSTANTS: Canonical k(infinity)')
99113 FORMAT ('  Eocent:',150(1x,f10.1),/)
99114 FORMAT ('Eocent/cm-1:',150(1x,f10.1))
99115 FORMAT ('Centr.Corr.:',150(1x,I10))
99016 FORMAT (' E-Eocent ',150(4X,I2.2,'...',I2.2))
99018 FORMAT (2x,f10.1,150(1x,1pe10.3))
99020 FORMAT ('   Eo(kcal):',150(1x,f10.2))
99021 FORMAT ('   Eo(cm-1):',150(1x,f10.1))
99022 FORMAT ('     Eo(kJ):',150(1x,f10.2))
99023 FORMAT (' delH(cm-1):',150(1x,f10.1))
99024 FORMAT (' delH(kcal):',150(1x,f10.1))
99025 FORMAT ('   delH(kJ):',150(1x,f10.1))
99026 FORMAT ('Rx_Path*Qel:',150(1x,f10.1))
99027 FORMAT (A100)
99028 FORMAT ('        ',4x,150(4X,I2.2,'-->',I2.2))
99029 FORMAT ('     k(inf):',150(1x,1pe10.2))
99030 FORMAT (/)
99034 FORMAT ('    TS_Name:',150(1x,A10))
99045 FORMAT ('  2D-moment:',150(1x,f10.1))
99046 FORMAT ('   ext_symm:',150(1x,I10))
99047 FORMAT ('Qel_partfxn:',150(1x,F10.2))
99048 FORMAT ('   opt_isom:',150(1x,I10))
99049 FORMAT ('     A(inf):',150(1x,1pe10.1))
99050 FORMAT (' E(inf)kcal:',150(1x,f10.1))
99051 FORMAT ('   E(inf)kJ:',150(1x,f10.1))
99052 FORMAT ('     Source:',150(1x,I10))
99053 FORMAT ('-----------',/,
     &'SOURCE: 0=> ILT method; 1=> sums of states file; 2=> k(E) file;',
     &' -2=> for reverse rxn; 3=> CRP file.',/,
     &'CENTR.CORR.: 0=> none; 1=> 1-D adiabatic rotor;',
     &' 2=> 2-D adiabatic rotor; -2=> legacy (not recommended).',/,
     &'WARNING: number of energy grains where Reactant Well has ',
     &'states, but Product Well has none to match [set k(E) = 0].')
99057 FORMAT (/,'Well number: ',I2,5x,'Name: ',A10,/,
     &      'Reference energy for collision frequency:',0pf8.1,' cm-1',
     &        /,'Collision parameters:',/,8(1X,1pe10.3))
99059 FORMAT (/'RESULTS FOR',1pe7.0,' TRIALS')
99060 FORMAT ('        deltaH(cm-1): ',75(5x,0pf10.1,18x))
99063 FORMAT ('                Name:    ',75(5x,A10,18x))
99161 FORMAT (3X ,'1/cc  Collisns(Mol=1)',
     &        75(: 3x,'Fract(',I2.2,') +/-Error  <Evib(',I2.2,')>'))
99165 FORMAT (/,/,'****Monoenergetic****'/'    Ttrans  = ',f6.1,
     &        ' K',/,'Initial Molecule: ',A10,/,'E(initial) =',
     &        0pf10.1,' cm-1',/,'Starting random number seed = ',
     &        I10,/,'Ending   random number seed = ',I10,/)
99166 FORMAT (/,/,'****Shifted Thermal (at Tvib)****'/'     Tvib   = ',
     &        f6.1,
     &        ' K',/,'     Ttrans = ',f6.1,' K',/,'Initial Molecule: ',
     &        A10,/,'E(shift) =',
     &        0pf10.1,' cm-1',/,'Starting random number seed = ',
     &        I10,/,'Ending   random number seed = ',I10,/)
99167 FORMAT (/,/,'****Chemical Activation (at Tvib)****'
     &        /'     Tvib   = ',
     &        f6.1,' K',/,'     Ttrans = ',f6.1,' K',/,
     &        'Excited Molecule: ', A10,/,
     &        '   Produced from: ', A10,/,
     &        '   E(shift) =', 0pf10.1,' cm-1'/,
     &        'Starting random number seed = ',I10,/,
     &        'Ending   random number seed = ',I10,/)
99168 FORMAT (/,/,'**** Energy Distribution from an EXTERNAL File ****',
     &        /,'     Ttrans = ',f6.1,' K',/,
     &        A10,/,'Starting random number seed = ',
     &        I10,/,'Ending   random number seed = ',I10,/)
99081 FORMAT ('IVR  v(IVR):',150(1x,f10.1),' cm-1')
99087 FORMAT ('IVR  v(ave):',150(1x,f10.1),' cm-1')
99082 FORMAT ('IVR  coll k:',150(1x,1pe10.3))
99186 FORMAT ('IVR  thresh:',150(1x,f10.1))
99083 FORMAT ('IVR civr(1):',150(1x,1pe10.3))
99084 FORMAT ('IVR civr(2):',150(1x,1pe10.3))
99085 FORMAT ('IVR civr(3):',150(1x,1pe10.3))
99086 FORMAT (' koSC(cc/s):',150(1x,1pe10.2))
99088 FORMAT ('TUN   vimag:',150(1x,0pf10.1))
99089 FORMAT ('React. [B] :',150(1x,1pe10.2))
99090 FORMAT (' WARNING!!!:',150(1x,I10))

      END SUBROUTINE write_2

!---------------------------------------------------------------------------------------
      SUBROUTINE rxn_input(Nforward)
!
!     Read parameters for multiwell reactions
!
c     Reaction parameters
c
c     Nforward= number of reaction parameters to be input (forward rxns, only)
c     Mol     = index of reactant molecule
c     ito     = Jto      = index of product molecule
c     TS     = TSname = Name of transition state (up to 10 characters)
c     RR     = TSmom  = moment of inertia (amu.ang^2, gram.cm^2, cm-1, MHz, GHz)
c     j     = TSsym  = external symmetry number for TS
c     qele     = TSele      = electronic partition function for TS
c     l     = TSopt      = number of optical isomers for TS
c     AA     = Afact     = A-factor for reaction (units: s-1)
c     EE     = Eo     = critical energy, relative to ZPE of reactant
c
c       KEYWORDS (5 are expected):
c
c         'NOCENT'   no centrifugal correction             [NCENT=  0]
c         'CENT1'    1-D adiabatic rotor centrifugal       [NCENT=  1]
c         'CENT2'    2-D adiabatic rotor centrifugal       [NCENT=  2]
c         'CENT3'    2-D adiabatic rotor centrifugal       [NCENT=  3]
c         'CENTX'    LEGACY centrifugal correction         [NCENT= -2] (not recommended)
c
c         'NOREV'    neglecting the reverse reaction        [Jrev = Jr]
c         'REV'      calculating reverse rxn rate; after calculation,
c                    Jrev=1 is changed to Jrev=2
c
c         'FAST'     neglecting limitations due to IVR
c         'SLOW'     including IVR limitations; a following line contains 
c                        'SLOW', pivr(Mol,i), vivr(Mol,i), vave(Mol,i), 
c                             tivr(Mol,i), civr(Mol,i,1), civr(Mol,i,2), 
c                             civr(Mol,i,3)
c
c         'NOTUN'    for neglecting tunneling
c         'TUN'      for including tunneling via unsymmetrical Eckart barrier;
c                       a following line contains: 'TUN', tunnel(Mol,i)
c
c          'ILT'     Inverse laplace transform method for k(E)
c          'SUM'     File containing sums of states (usually generated by Densum)
c          'CRP'     File containing CRP (usually generated by SCTST)
c          'RKE'     File containing k(E), generated externally
c          'REVRKE'  k(E) for reverse reaction calculated from k(E) file for forward rxn
c
c                Jread = 0 for ILT method
c                      = 1 for reading sums of states from external file
c                        *** filename: same name as TSame +'.dens'
c                        ***   e.g. 'TS-1.dens'
c                      = 2 for reading k(E) from external file
c                        *** filename: same name as TSame +'.rke'
c                        ***   e.g. 'TS-1.rke'
c                      = -2 for calculating reverse k(E) from external file for forward k(E)
c                        *** filename: same name as TSame +'.rke'
c                        ***   e.g. 'TS-1.rke'
c                      = 3 for reading CRP of states from external file CREATED BY SCTST or PARSCTST
c                        *** filename: same name as TSame +'.dens'
c                        ***   e.g. 'TS-1.dens'
c.......................................................................
c
c
c     Nchan(Mol) = number of rxn channels originating from a given molecule
c.......................................................................
      USE declare_mod
      USE utili_mod
      USE subs_mod
      IMPLICIT NONE
      SAVE
      INTEGER, INTENT(IN) :: Nforward
      REAL(8) :: Af, Bx, xd, qele
      INTEGER :: jbim, bt
      INTEGER :: i, j, l, ito, jivr, jr, jtun, KK, kkey, ndoubt
      INTEGER :: NSTOP      

      Nreac = Nforward      
      ntunflag = 0                     ! flag for tunneling
      nivrflag = 0                     ! flag for slow IVR
      Pressflag = 0                    ! flag for reactant [B]  (for microcanonical pseudo-first order rxn)

      DO 800 i = 1 , Nforward
 
         READ (KIN,*) Mol , ito , TS , RR , j , qele ,
     &                 l , AA , EE , (KEYWORD(II),II=1,5)                     ! *** LINE 14 ***
         DO II = 1 , 5
            call ucase( KEYWORD(II) )
         END DO
         WRITE (KOUT,99075) i , Mol , ito , TS , RR , j , qele , 
     &                            l , AA , EE , (KEYWORD(II),II=1,5)
     
         call rotunits( Rotatunits , RR , 'AMUA' )     ! convert rotational units to 'AMUA'

         Nreaindex(i)=Mol ! keep track of reactants and products so that
         Nproindex(i)=ito ! the initial distribution and Eor can be printed
                          ! correctly in multiwell.array

         Nchan(Mol) = Nchan(Mol) + 1
         Jto(Mol,Nchan(Mol)) = ito
         TSname(Mol,Nchan(Mol)) = TS
         TSmom(Mol,Nchan(Mol)) = RR
         TSsym(Mol,Nchan(Mol)) = j
         TSele(Mol,Nchan(Mol)) = qele
         TSopt(Mol,Nchan(Mol)) = l
         Afac(Mol,Nchan(Mol)) = AA
         IF(Afac(Mol,Nchan(Mol)).lt.0) THEN                 ! (for microcanonical pseudo-first order rxn)
           Press(Mol,Nchan(Mol))= -Afac(Mol,Nchan(Mol))
           Pressflag=Pressflag+1                            ! (for microcanonical pseudo-first order rxn)
          ELSE
           Press(Mol,Nchan(Mol))=1
         ENDIF
         IF ( Eunits.EQ.cmm1 ) THEN
            Eo(Mol,Nchan(Mol)) = EE                     ! cm-1 units
         ELSEIF ( Eunits.EQ.kcal ) THEN
            Eo(Mol,Nchan(Mol)) = EE*caltocm             ! Convert kcal/mole to cm-1
         ELSEIF ( Eunits.EQ.kjou ) THEN
            Eo(Mol,Nchan(Mol)) = EE*jtocm               ! Convert kJ/mole to cm-1
         ENDIF
         Hts(Mol,ito) = HMol(Mol) + Eo(Mol,Nchan(Mol))  ! Rxn barrier (or transition state) enthalpy same in both directions
c
c      KEYWORDS
c
c
c        DEFAULTS
c
         NCENT(Mol,Nchan(Mol)) = 2        ! CENT2 centrifugal correction
         Jr = 1                                       ! reversible reaction
         jtun = 0                               ! no tunneling
         jivr = 0                               ! fast ivr
         itun(Mol,Nchan(Mol)) = 0
         Jread(Mol,Nchan(Mol)) = 1        ! read from sum of states file (suffix .dens)

         kkey = 0

         DO 701 KK = 1 , 5
           call ucase( KEYWORD(KK) )

           IF     (KEYWORD(KK) .EQ. 'ILT') THEN   ! k(E) from Inverse Laplace Transform method
              Jread(Mol,Nchan(Mol)) = 0
              kkey = kkey + 1
           ELSEIF (KEYWORD(KK) .EQ. 'SUM') THEN   ! external sums of states file
              Jread(Mol,Nchan(Mol)) = 1
              kkey = kkey + 1
           ELSEIF (KEYWORD(KK) .EQ. 'RKE') THEN  ! external k(E) file
              Jread(Mol,Nchan(Mol)) = 2
              kkey = kkey + 1
           ELSEIF (KEYWORD(KK) .EQ. 'CRP') THEN   ! external CRP file (like a sums file) generated by SCTST
              Jread(Mol,Nchan(Mol)) = 3
              kkey = kkey + 1
           ELSEIF (KEYWORD(KK) .EQ. 'REVRKE') THEN  ! external k(E) file used to calculate REVERSE reaction
              Jread(Mol,Nchan(Mol)) = -2
              kkey = kkey + 1
           ELSEIF (KEYWORD(KK) .EQ. 'REV') THEN   ! Include reverse reaction
              Jr = 1
              kkey = kkey + 1
           ELSEIF (KEYWORD(KK) .EQ. 'NOREV') THEN
              Jr = 0
              kkey = kkey + 1
           ELSEIF (KEYWORD(KK) .EQ. 'TUN') THEN   ! Include tunneling
              ntunflag = 1
              jtun = 1
              kkey = kkey + 1
           ELSEIF (KEYWORD(KK) .EQ. 'NOTUN') THEN  ! No tunneling
              jtun = 0  
              kkey = kkey + 1
           ELSEIF (KEYWORD(KK) .EQ. 'SLOW') THEN   ! Include slow IVR
              nivrflag = 1
              jivr = 1
              kkey = kkey + 1
           ELSEIF (KEYWORD(KK) .EQ. 'FAST') THEN   ! Fast IVR
              jivr = 0
              kkey = kkey + 1
           ELSEIF (KEYWORD(KK) .EQ. 'NOCENT') THEN   ! No centrifugal correction
              NCENT(Mol,Nchan(Mol)) = 0
              kkey = kkey + 1
           ELSEIF (KEYWORD(KK) .EQ. 'CENTX') THEN   ! Include LEGACY centrifugal correction
              NCENT(Mol,Nchan(Mol)) = -2
              kkey = kkey + 1
           ELSEIF (KEYWORD(KK) .EQ. 'CENT1') THEN   ! Include 1-D ROTOR centrifugal correction
              NCENT(Mol,Nchan(Mol)) = 1
              kkey = kkey + 1
           ELSEIF (KEYWORD(KK) .EQ. 'CENT2') THEN   ! Include 2-D ROTOR centrifugal correction
              NCENT(Mol,Nchan(Mol)) = 2
              kkey = kkey + 1
           ELSEIF (KEYWORD(KK) .EQ. 'CENT3') THEN   ! Include 3-D ROTOR centrifugal correction
              NCENT(Mol,Nchan(Mol)) = 3
              kkey = kkey + 1
           ENDIF

 701     CONTINUE

           NSTOP = 0
           IF ( kkey .NE. 5 ) THEN
             write(6,*)'*** LESS THAN 5 REQUIRED KEYWORDS ***'
             write(KOUT,*)'*** LESS THAN 5 REQUIRED KEYWORDS ***'
             NSTOP = 1
           ENDIF
           IF ((Jread(Mol,Nchan(Mol)).eq.0).AND.(jtun.EQ.1))THEN
             write(6,*)'KEYWORDS: ILT CANNOT BE USED WITH TUN'
             write(KOUT,*)'KEYWORDS: ILT CANNOT BE USED WITH TUN'
             NSTOP = 1
           ENDIF
           IF ((Jread(Mol,Nchan(Mol)).eq.0).AND.
     &           (NCENT(Mol,Nchan(Mol)).NE.0))THEN
             write(6,*)'KEYWORDS: ILT CANNOT BE USED WITH CENT'
             write(KOUT,*)'KEYWORDS: ILT CANNOT BE USED WITH CENT'
             NSTOP = 1
           ENDIF
           IF ( ( NCENT(Mol,Nchan(Mol)) .NE. 0 ) .AND.
     &           (TSmom(Mol,Nchan(Mol)) .LE. 0.0))THEN
             write(6,*)'CENT KEYWORD REQUIRES TSmom >0'
             write(KOUT,*)'CENT KEYWORD REQUIRES TSmom >0'
             NSTOP = 1
           ENDIF
           IF((Jread(Mol,Nchan(Mol)).eq.2).AND.(jtun.EQ.1))THEN
             write(6,*)'KEYWORDS: RKE CANNOT BE USED WITH TUN'
             write(KOUT,*)'KEYWORDS: RKE CANNOT BE USED WITH TUN'
             NSTOP = 1
           ENDIF
           IF((Jread(Mol,Nchan(Mol)).eq.3).AND.(jtun.EQ.1))THEN
             write(6,*)'KEYWORDS: TUN CANNOT BE USED WITH CRP'
             write(6,*)'   BECAUSE CRP ALREADY INCLUDES TUNNELING'
             write(KOUT,*)'KEYWORDS: TUN CANNOT BE USED WITH CRP'
             write(KOUT,*)'   BECAUSE CRP ALREADY INCLUDES TUNNELING'
             NSTOP = 1
           ENDIF
           IF((NCENT(Mol,Nchan(Mol)).NE. 0).AND.(jtun.EQ.1))THEN
             write(6,*)'KEYWORDS: CENT CANNOT BE USED WITH TUN'
             write(KOUT,*)'KEYWORDS: CENT CANNOT BE USED WITH TUN'
             NSTOP = 1
           ENDIF
           IF((Jread(Mol,Nchan(Mol)).eq.-2).AND.(jtun.EQ.1))THEN
             write(6,*)'KEYWORDS: RKE CANNOT BE USED WITH TUN'
             write(KOUT,*)'KEYWORDS: RKE CANNOT BE USED WITH TUN'
             NSTOP = 1
           ENDIF

           IF ( NSTOP .NE. 0 ) THEN
             write(*,*)'***** FATAL ERROR(S) (NSTOP .NE. 0)'
             write(KOUT,*)'***** FATAL ERROR(S) (NSTOP .NE. 0)'
             STOP
           END IF
             
          IF (jivr .EQ. 1 ) THEN
             iivr(Mol,Nchan(Mol)) = 1
             READ(KIN,*) KEYTEMP, vivr(Mol,Nchan(Mol)), 
     &          vave(Mol,Nchan(Mol)) ,
     &          pivr( Mol , Nchan(Mol) ) , tivr( Mol , Nchan(Mol) ) ,
     &          ( civr( Mol , Nchan(Mol) , ii ) , ii = 1 , 3 )
          ENDIF

          IF (jtun .EQ. 1 ) THEN
             itun(Mol,Nchan(Mol)) = 1
             READ(KIN,*) KEYTEMP , vimag(Mol,Nchan(Mol))
          ENDIF

         IF ( Jr.EQ.1 ) THEN                    ! Jr=1 for reversible rxn
            Jrev(Mol,Nchan(Mol)) = 1
            IF ( ito .GT. NWells ) THEN
               Jrev(Mol,Nchan(Mol)) = 0         ! If 'ito' is an irreversible product molecule
            ELSE
               Hts(ito,Mol) = Hts(Mol,ito)      ! Rxn barrier (or transition state) enthalpy same in both directions
               Nreac = Nreac + 1                ! summing to get total number of rxns
               Nchan(ito) = Nchan(ito) + 1
               Jto(ito,Nchan(ito)) = Mol
               TSname(ito,Nchan(ito)) = TS
               TSmom(ito,Nchan(ito)) = RR
               Afac(ito,Nchan(ito)) = AA
               Eo(ito,Nchan(ito)) = Hts(ito,Mol) - HMol(ito)
               NCENT(ito,Nchan(ito)) = NCENT(Mol,Nchan(Mol))
               TSsym(ito,Nchan(ito)) = j
               TSele(ito,Nchan(ito)) = qele
               TSopt(ito,Nchan(ito)) = l
               Jrev(ito,Nchan(ito)) = 0
               iivr(ito,Nchan(ito)) = iivr(Mol,Nchan(Mol))
               vivr(ito,Nchan(ito)) = vivr(Mol,Nchan(Mol))
               vimag(ito,Nchan(ito)) = vimag(Mol,Nchan(Mol))
               iivr(ito,Nchan(ito)) = iivr(Mol,Nchan(Mol))
               pivr(ito,Nchan(ito)) = pivr(Mol,Nchan(Mol))
               tivr(ito,Nchan(ito)) = tivr(Mol,Nchan(Mol))
               vave(ito,Nchan(ito)) = vave(Mol,Nchan(Mol))
               civr(ito,Nchan(ito),1) = civr(Mol,Nchan(Mol),1)
               civr(ito,Nchan(ito),2) = civr(Mol,Nchan(Mol),2)
               civr(ito,Nchan(ito),3) = civr(Mol,Nchan(Mol),3)
               itun(ito,Nchan(ito)) = itun(Mol,Nchan(Mol))
              
               IF ( Jread(Mol,Nchan(Mol)).EQ.2 ) THEN    ! will open and read k(E) file (external file)
                  Jread(ito,Nchan(ito)) = -2             ! will open and read k(E) file for forward rxn and then calculate reverse k(E)
               ELSE
                  Jread(ito,Nchan(ito)) = Jread(Mol,Nchan(Mol))  ! i.e. if sums read from file for forward rxn, sums read from same file for reverse rxn.
               ENDIF
            ENDIF  ! ito
         ENDIF  ! Jr
 
 800  CONTINUE ! end of Nforward reactions
 
c------------------------------------------------------------------------------------------
c     SUPPLEMENTARY REACTIONS (IRREVERSIBLE)
c------------------------------------------------------------------------------------------
c     First (unimolecular) or second order (bimolecular)
c     BIMOLECULAR REACTIONS between Mol and bath gas,
c     Pseudo-first-order rate constant for reaction with bath gas is
c     the product:
c                 ksup(Mol,n)*[bath gas]
c
c     where ksup(Mol,n) = Afc(Mol,n) * exp( -Bsr(Mol,n)/Temp )
c
c     MaxSup      :     maximum allowed supplementary rxns for each Well (index Mol)
c     nsux  :     total number of supplementary reactions (all wells)
c     nsmax :     number of supplementary rxns  for a Well (index Mol) 
c     nsup        :     index number for supplentary rxn for a Well (index Mol)
c     SupName     :     name of supplementary reaction
c     sto         :     product index
c     norder      :     reaction order
c     Afc         :     Arrhenius A-factor (s-1 or cm3 s-1)
c     Bsr         :     E(activation)/R
c     dde         :     Quench energy (energy lost to environment)
c     ksup        :     Supplemntary rate constant = Afc*exp(-Bsr/T)
c
       DO Mol = 1, NWells
         nsmax(Mol) = 0                               ! initialize number of supplementary reactions for each Well
       END DO
       
       READ (KIN,*) dummy
       CALL UCASE ( dummy )
       IF ( dummy .EQ. 'MORERXN' ) THEN
         READ (KIN,*) nsux                                  ! total number of supplementary reactions to be read
         DO jbim = 1, nsux
           READ (KIN,*) Mol, bt , dummy, ndoubt, Af, Bx, xd       ! from, to, Name, rxn_order, Afactor, Ea/R , vib energy in Well(ito) lost to bath when quenching
           IF ( Mol.GT. NWells ) THEN                       ! can only react with Wells
             Write(*,*) 'FATAL: Supplementary rxns only from Wells'
             STOP
           ENDIF
           nsmax(Mol) = nsmax(Mol) + 1
           nsup(Mol,nsmax(Mol)) = nsmax(Mol)
           IF ( nsmax(Mol) .GT. MaxSup ) THEN                     ! maximum for each Mol: MaxSup
             Write(*,*) 'FATAL: Too many supplementary rxns, Mol:', Mol
             STOP
           ENDIF
           SupName( Mol,nsmax(Mol) ) = dummy                ! Rxn Name
           sto( Mol,nsmax(Mol) ) = bt                             ! index number of product (INTEGER)
           norder(Mol,nsmax(Mol) ) = ndoubt                 ! reaction order
           Afc( Mol,nsmax(Mol) ) = Af                       ! Afactor
           Bsr( Mol,nsmax(Mol) ) = Bx                       ! Ea/R activation energy
           dde( Mol,nsmax(Mol) ) = xd                       ! Quench energy (electronic lost from product to bath)
           ksup( Mol,nsmax(Mol) ) = Af*exp( -Bx/Temp )
         END DO      ! jbim
       ELSE
         BACKSPACE (KIN)
       ENDIF

99075 FORMAT(3(I2,1X),A10,1X,1PE10.2,1X,I2,1X,0PF10.2,1X,I2,
     &  2(1PE10.2,1X),5(A7,1x))
      
      END SUBROUTINE rxn_input

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
      END MODULE I_O_mod

      