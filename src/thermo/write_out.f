c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Copyright (C) 2011 Andrea Maranzana and John R. Barker
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

      Subroutine write_out

      include 'declare.inc'
      REAL(8) Dum, qqele, RK, xa, xb, xc
      DIMENSION qqele(10)
      INTEGER(4) LEN_TRIM
      SAVE

      DO Temp=1,Nt

      WRITE (3,99002) AVERSION , ADATE   ! Declared in "declare.inc"
c      WRITE (3,99902) AVERSION , ADATE
      
      CALL DateTime(3)

      
      DO  i=1,Ns
         
         WRITE (3,99014) MOLNAME(i),
     &           AMU(i) , ( ADJUSTR(ATYPE(i,ne)) , Natom(i,ne) , 
     &           ne=1,Nelement(i) )
         IF ( nlines(i) .GT. 0 ) THEN  
            DO ll = 1 , nlines(i)
               WRITE(3, 99224) TITLELINE( i , ll )
            END DO
         ENDIF

         WRITE (3,99022)  
         IF (nso(i) .GT. 0) write (3,99922) 			! if ITYPE=RSO was used for this species
         WRITE (3,99822) (Elev(i,ne) , gele(i,ne) , ne=1,Nele(i))
         IF (Warning(i).eq.1) then 
           WRITE(3,*) "**** WARNING: Used 10 energy levels, only. ****"
           WRITE(3,*) ' '
         ENDIF

         IF (REPROD(i).EQ.CTST .AND. NCRP.EQ.0 ) THEN
           WRITE (3,99050) vimag , Vf , Vr                       ! TUNNELING PARAMETERS
           IF ( vimag.LT.20. ) WRITE (3,99051)
           IF ( Vf.LT.1. ) WRITE (3,99051)
           IF ( Vr.LT.1. ) WRITE (3,99051)
         ENDIF
         
         IF ( N(i).GT.0 ) THEN
            DO 20 II = 1 , N(i)

               IF ( IDOF(i,II).EQ.VIB ) THEN                        ! VIBRATION
                  CALL qmorse(300.d+00,WE(i,II),ANH(i,II),NG(i,II),
     &                   qq,Cvib,Svib,Hvib,zap)
                  WRITE (3,99010) MODE(i,II) , WOBS(i,II) , We(i,II) ,
     &                      ANH(i,II) , zap , NG(i,II)
               ELSEIF ( IDOF(i,II).EQ.BOX ) THEN                    ! PARTICLE IN A BOX
                  zap = W(i,II)
                  WRITE (3,99110) MODE(i,II) , WOBS(i,II) , W(i,II) ,
     &                      zap , NG(i,II)
               ELSEIF ( IDOF(i,II).EQ.ROT ) THEN                ! CLASSICAL ROTOR
                  BB = 16.85763D+00/AMOM(i,II)                  ! rotational constant (cm-1)
                  WRITE (3,99011) MODE(i,II) , WE(i,II) , ANH(i,II) ,
     &                            NG(i,II) , BB
               ELSEIF ( IDOF(i,II).EQ.QRO ) THEN                ! QUANTUM ROTOR
                  BB = 16.85763D+00/AMOM(i,II)                  ! rotational constant (cm-1)
                  WRITE (3,99012) MODE(i,II) , WE(i,II) , ANH(i,II) ,
     &                            NG(i,II) , BB
               ELSEIF ( IDOF(i,II).EQ.RSO ) THEN                ! ROTATION + SPIN-ORBIT
                  WRITE (3,99070) MODE(i,II) , MULT(i) , LAMBDA(i) ,
     &                     ASO(i), BSO(i), D(i), NRSO(i)
              ELSEIF (IDOF(i,II).EQ.GOR) THEN                   ! GORIN
               AMOM(i,II)=I2DGorin(Temp)
               WE(i,II)=AMOM(i,II)
               BB = 16.85763D+00/AMOM(i,II)
               if (AMOM(i,II).le.10)  THEN
                  WRITE (3,99012) MODE(i,II) , WE(i,II) , ANH(i,II) ,
     &                            NG(i,II) , BB
               ELSE
                  WRITE (3,99011) MODE(i,II) , WE(i,II) , ANH(i,II) ,
     &                            NG(i,II) , BB
               ENDIF
               ELSEIF ( IDOF(i,II).EQ.FIT) THEN                 ! HINDRANCE
                  AMOM(i,II) = WE(i,II)*gamma(Temp)
                  BB = 16.85763D+00/AMOM(i,II)                  ! rotational constant (cm-1)
                  if (AMOM(i,II).le.10)  THEN
                  WRITE (3,99012) MODE(i,II) , AMOM(i,II) , ANH(i,II) ,
     &                            NG(i,II) , BB
                    ELSE
                  WRITE (3,99011) MODE(i,II) , AMOM(i,II), ANH(i,II) ,
     &                            NG(i,II) , BB
                  endif
              ELSEIF ( IDOF(i,II).EQ.HIN ) THEN                ! SYMMETRICAL HINDERED ROTOR 

                  WRITE (3,99001) MODE(i,II) , W(i,II) , AMOM(i,II) ,
     &            NG(i,II), NSIG(i,II), B(i,II), zpe(i,II), VV(i,II)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            ELSEIF ( IDOF(i,II).EQ.HRD ) THEN                  ! GENERAL, UNSYMMETRICAL HINDERED ROTOR

             WRITE (3,99040) MODE(i,II), NG(i,II), zpe(i,II)
             WRITE (3,99041) VHRD(i,II), NSV(i,II),Phav(i,II),
     &                   (CV(i,II,LL), LL=1, NCV(i,II))
             WRITE (3,99042) BHRD(i,II), NSB(i,II),Phab(i,II),
     &                   (CB(i,II,LL), LL=1, NCB(i,II))

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

               ELSEIF ( IDOF(i,II).EQ.TOP ) THEN               ! SYMMETRIC TOP
                  BJ = 16.85763D+00/WE(i,II)
                  BK = 16.85763D+00/ANH(i,II)
                  WRITE (3,99113) MODE(i,II) , WE(i,II) , 
     &                  BJ , ANH(i,II) , BK , NG(i,II)

               ELSEIF ( IDOF(i,II).EQ.QVB ) THEN               ! ANHARMONIC VIBRATIONAL PARTITION FUNCTION, external file
                 WRITE(3,99060) MODE(i,II), TITLE4(i,II)
                 WRITE(3,99062) Egrain1(i,II), Emax2(i,II), 
     &                   KEYWORD2(i)

               ELSEIF ( IDOF(i,II).EQ.CRP ) THEN               ! SEMI-CLASSICAL TST, external file
                 WRITE(3,99061) MODE(i,II), TITLE5(i,II)
                 WRITE(3,99063) Egrain1(i,II), Emax2(i,II), Vf, Vr, 
     &                   KEYWORD2(i), VPTx(i)

              ENDIF
20          ENDDO  ! II
       ENDIF
       ENDDO   ! i

      DDelH298 = 0.0d+00
      DDelG298 = 0.0d+00
      DO i = 1 , Ns
        IF ( REPROD(i).EQ.REAC ) THEN
           DDelH298 = DDelH298 - delH298(i)                      ! Delta H(298) for reaction
           DDelG298 = DDelG298 - delG298(i)                      ! Delta G(298) for reaction
        ELSEIF ( REPROD(i).EQ.PROD .OR. REPROD(i).EQ.CTST ) THEN
           DDelH298 = DDelH298 + delH298(i)                      ! Delta H(298) for reaction
           DDelG298 = DDelG298 + delG298(i)                      ! Delta G(298) for reaction
        ENDIF
      END DO

      IF ( Eunits.EQ.'KCAL' ) THEN
         WRITE (3,99004) UNITS 
      ELSEIF ( Eunits.EQ.'KJOU' ) THEN
         WRITE (3,99005) UNITS 
      ELSEIF ( Eunits.EQ.'CM-1' ) THEN
         WRITE (3,99055) UNITS 
      ENDIF

       WRITE (3,99016)

       DO i=1,Ns
         WRITE (3,99017) REPROD(i) , MOLNAME(i) ,
     &     AMU(i) , DelH(i) , delH298(i) , delG298(i) , zzpe(i) , 
     &     Sym(i) , Sopt(i)
       ENDDO  ! i

      IF ( Nflag .EQ. 1) THEN                   ! reaction thermo
         WRITE (3,99021)
         WRITE (3,99020) DDelH, DDelH298 , DDelG298 , Dzpe
      ENDIF

      IF (Nctst .EQ. 1 ) THEN                                       ! Table caption for Canonical TST rate constant
         WRITE (3,99106) UNITS, 1-Nreac, (REPROD(i), MOLNAME(i) ,i=1,Ns)
         WRITE (3,99103) (CAP,i=1,Ns)
      ELSE
         WRITE (3,99006) (REPROD(i), MOLNAME(i) ,i=1,Ns)             ! Table caption for Keq and thermodynamics output
         WRITE (3,99003) (CAP,i=1,Ns)
      ENDIF

      IF(hindrance.eq.0.AND.Gorin.eq.0) then            ! STANDARD OUTPUT: NOT GORIN OR FITTING
         DO j=1,Nt
           DO i=1,Ns
             CALL qelect(T(j),i,Nele(i),Elev,gele,qq,Cel,Sel,Hel)
             qqele(i) = qq
           END DO

           WRITE (3,99007) T(j) , T_AK(j) , T_AA(j) ,T_BB(j) ,T_DSR(j) ,
     &     T_DHR(j)/DIV , T_DCp(j) ,
     &     (T_DHR(j)-T(j)*T_DSR(j))/DIV , ( qqele(i) ,
     &     T_S(j,i) , T_Cp(j,i) , T_H(j,i)/DIV, i=1,Ns )
         ENDDO   ! end j

       CLOSE(3)
       STOP
      ENDIF                                                       ! end standard output

!                         GORIN and/or FITTING SECTIONS

       WRITE (3,99007) T(temp) , T_AK(temp) , T_AA(temp) ,T_BB(temp) ,
     & T_DSR(temp) , T_DHR(temp)/DIV , T_DCp(temp) ,
     & (T_DHR(temp)-T(temp)*T_DSR(temp))/DIV ,
     & (T_S(temp,i),T_Cp(temp,i),T_H(temp,i)/DIV, i=1,Ns)
         

        IF(Gorin.NE.0.OR.hindrance.ne.0) then
           if(hindrance.ne.0.AND.Gorin.eq.0) then
            write(*,*)
            write(3,*)
            write(*,*) "AUTOMATIC FITTING CALCULATION "
            write(3,*) "AUTOMATIC FITTING CALCULATION "
           endif

           if(hindrance.eq.0.AND.Gorin.ne.0) then
            write(*,*) 
            write(3,*) 
            write(*,*) "GORIN CALCULATION "
            write(3,*) "GORIN CALCULATION "
           endif

           if(hindrance.ne.0.AND.Gorin.ne.0) then
            write(*,*) 
            write(3,*) 
            write(*,*) "AUTOMATIC FITTING AND GORIN CALCULATION "
            write(3,*) "AUTOMATIC FITTING AND GORIN CALCULATION "
           endif

        write(*,*) 
        write(3,*)
        write(*,*) "T =",T(temp)
        write(3,*) "T =",T(temp) 
        ENDIF

        IF(Gorin.NE.0) then
         TypePOT=POT
         IF(POT.eq.'sMORSE') typePOT='Stiff MORSE'
         write(*,*) "   Potential:  ",TypePOT
         write(3,*) "   Potential:  ",TypePOT
         write(*,*) "   Frequency:  ", freq," cm-1"
         write(3,*) "   Frequency:  ", freq," cm-1"
!         write(*,*) "   2-D moment: ",I2d," amu A**2"
!         write(3,*) "   2-D moment: ",I2d," amu A**2"
         write(*,*) "   De Energy:  ", Deorig," ",Eunits
         write(3,*) "   De Energy:  ", Deorig," ",Eunits
         write(*,*) "   Equil. re:  ", re, " A"
         write(3,*) "   Equil. re:  ", re, " A"
         if(POT.eq.'VARSHNI') then
          write(*,*) "   Beta:       ", bet(temp), " A-2"
          write(3,*) "   Beta:       ", bet(temp), " A-2"
          else
          write(*,*) "   Beta:       ", bet(temp), " A-1"
          write(3,*) "   Beta:       ", bet(temp), " A-1"
          endif
         if(POT.eq.'sMORSE') then 
          write(*,*) "   stiff c:    ",c
          write(3,*) "   stiff c:    ",c
         endif 
         write(*,99023) I2dGorin(temp),rmax(temp)
         write(3,99023) I2dGorin(temp),rmax(temp)
         write(*,99024)  V(temp) , Eunits
         write(3,99024)  V(temp) , Eunits
         write(*,99025)  Veff(temp) , Eunits
         write(3,99025)  Veff(temp) , Eunits
        ENDIF

         IF(hindrance.ne.0) then
          write(*,*) "   Kexp       =  ",Kexp(temp)
          write(3,*) "   Kexp       =  ",Kexp(temp)
         IF(gamma(Temp).gt.1) then
         write(*,*) "**** WARNING: Gamma >1.0 (Maybe Kexp too big) ****"
         write(3,*) "**** WARNING: Gamma >1.0 (Maybe Kexp too big) ****"
         endif
        write(*,99026)  gamma(Temp), 1-gamma(Temp)**2                    ! Print Gamma and eta
        write(3,99026)  gamma(Temp), 1-gamma(Temp)**2                    ! Print Gamma and eta
        do i=1,hindrance
         write(*,99027) Weorig(i)*gamma(Temp)             !Print Optimized moment of inertia
         write(3,99027) WeOrig(i)*gamma(Temp)             !Print Optimized moment of inertia
        enddo
         ENDIF
       
       write(3,*) 
       write(3,*) 
       write(3,*) 

       ENDDO
       CLOSE(3)
       
       RETURN
 

99001 FORMAT (I5,2X,'Hind.Rot',2x,'Freq(har)=',F8.2,2x,'Mom=',F8.3,2x,
     &       'fold=',I2,2x,'symm=',I2,2x,'B=',F8.4,' cm-1',2x,
     &        'zpe(cm-1)=',F8.1,2x,'Uo=',F7.1,' cm-1')
99002 FORMAT (
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'//,
     &'                                  John R. Barker',/,
     &'        Thermo-',A7, '            University of Michigan',/,
     &'                                  Ann Arbor, MI 48109-2143',/,
     &'        ',A8,  '                  jrbarker@umich.edu',/,
     &'                                  (734) 763 6239',//,
     &'        http://clasp-research.engin.umich.edu/multiwell/'/)
99902 FORMAT (
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'//,
     &'Suggested Literature Citations:'//4x,
     &'a) J.R. Barker, T.L. Nguyen, J.F. Stanton, C. Aieta, M. Ceotto,',
     &/7x,'F. Gabas, T.J.D. Kumar, C.G.L.Li, L.L.Lohr, A. Maranzana,',
     &/7x,'N.F. Ortiz, J.M. Preses, J.M. Simmie, J.A. Sonk, and ',
     &'P.J. Stimac; ',
     &/7x,'MultiWell-',A7,' Software Suite; University of Michigan, ',
     &/7x,'Ann Arbor, Michigan, USA, ',A8,'; ',
     &/7x,'http://clasp-research.engin.umich.edu/multiwell/.',
     & //4x,
     &'  b) J.R. Barker, Int. J. Chem. Kinetics, 33, 232-45 (2001).',//,
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'/)
99003 FORMAT (8x,'T(K)',5x,'Kequil',5x,'A(T)',7x,'B(T)',4x,'DelS(rxn)',
     &       1x,'DelH(rxn)',1x,'DelCp(rxn)',1x,'DelG(rxn)',3x,6(A33,8x))
99103 FORMAT (8x,'T(K)',6x,'k(T)  ',4x,'A(T)',7x,'B(T)',4x,'DelS(rxn)',
     &       1x,'DelH(rxn)',1x,'DelCp(rxn)',1x,'DelG(rxn)',3x,6(A33,8x))
99004 FORMAT (//'*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  ',
     &'*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *'//, 
     &'THERMODYNAMIC QUANTITIES',10x,'Std. State = ',A12,/,3x,
     &        'Entropy & Cp units: cal/K/mole'/3x,
     &        'DelS         units: cal/K/mole'/3x,
     &        '[H(T)-H(0)]  units: kcal/mole'/3x,
     &        'DelH & DelG  units: kcal/mole'/)
99005 FORMAT (//'*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  ',
     &'*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *'//, 
     &'THERMODYNAMIC QUANTITIES',10x,'Std. State = ',A12,/,3x,
     &        'Entropy & Cp units: J/K/mole'/3x,
     &        'DelS         units: J/K/mole'/3x,
     &        '[H(T)-H(0)]  units: kJ/mole'/3x,
     &        'DelH & DelG  units: kJ/mole'/)
99055 FORMAT (//'*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  ',
     &'*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *'//, 
     &'THERMODYNAMIC QUANTITIES',10x,'Std. State = ',A12,/,3x,
     &        'Entropy & Cp units: CM-1/K'/3x,
     &        'DelS         units: CM-1/K'/3x,
     &        '[H(T)-H(0)]  units: CM-1'/3x,
     &        'DelH & DelG  units: CM-1'/)
99006 FORMAT (//21x,'Kequil = A*exp(B/T)',
     &        48x,6(:2x,'(',A4,') ', A32))
99106 FORMAT (//15x,'RATE k(T) = A*exp(B/T)',1x,
     &       '(',A12,')**(',I2,') * sec-1',
     &        22x,6(:2x,'(',A4,') ', A32))
99007 FORMAT (1x,f11.2, 2(1x,1pe10.3), 2x,0pf7.0,1x,0pf9.2,1x,0pf10.2, 
     &        1x,0pf10.2,1x,0pf10.2, 
     &        6(0pf8.2,1x,f8.2,1x,f8.2,1x,f8.2,6x)) 
99010 FORMAT (I5,2X,'Vibrator: Fund= ',F8.2,2x,
     &        'Harm= ',F8.2,2x,'Anh(cm-1)= ',F7.2,2X,
     &        'zpe(cm-1)= ',F8.1,2x,'   g= ',I2)
99110  FORMAT (I5,2X,'P-in-Box: Fund= ',F8.2,2x,
     &        'Harm= ',F8.2,2x,
     &        'zpe(cm-1)= ',F8.1,2x,'   g= ',I2)
99011 FORMAT (I5,2X,'Clas.Rot',2x,'AmuAng^2=',F8.2,2x,'Symm=',F4.0,2x,
     &        ' dim=',I2,2x,'B=',F8.4,' cm-1')
99012 FORMAT (I5,2X,'Quan.Rot',2x,'AmuAng^2 =',F8.2,2x,'Symm =',F4.0,2x,
     &        ' dim =',I2,2x,'B =',F8.4,' cm-1')
99113 FORMAT (I5,2X,'Symm-Top',2x,'2D-Moment =',F8.2,' amua, B2 =',F8.4,
     &       ' cm-1  ;  1D-Moment =',F8.2,' amua, B1 =',F8.4,' cm-1',2x,
     &       ';  Symm. No. = ',I1)
99014 FORMAT (//A32, 5X,'Molecular Wt. =', F8.3,
     &       '  ;  Empirical Formula:', 15(1x, A6 ,1x, I2 ,1x) ) 
99016 FORMAT ('Species Properties'/,
     &  54x,'----Enthalpy----',/,
     &  ' TYPE  NAME',32x,'Mol.Wt.    (O K)  (298.15 K) DelG(298)    ',
     &  'ZPE(cm-1) Symmetry  opt.isomers')
99017 FORMAT(1x,A4,2x,A32,1x,F9.3,2(2x,F8.2),2x,F8.2,
     & 6x,F9.1,2x,I5,I10)
99020 FORMAT (15x,'NET (prod - react):',
     &        16x,F9.2,1x,F9.2,1x,F9.2,4x,F11.1)
99021 FORMAT (52x,'------- ---------  --------      ---------')
99022 FORMAT (/,4X, 'Electronic Levels (degeneracies) ')
99822 FORMAT ( 10(10x,F10.2,1x,'cm-1 (',I2,')',/) )
99922 FORMAT(4x,'*** NOTE that RSO type d.o.f. (Spin-Orbit+rotation) ',
     &       'is invoked for this species ***',/, 
     &  4x,'*** CAREFUL: to avoid double-counting of rotation, do NOT',
     &  ' also use ROT, or QRO, for this species.')
99023 FORMAT ('>>> Gorin external 2-D moment =',F9.4,' amu A**2 <<<',/
     &        '    Rmax       = ',F8.3,' A')
99024 FORMAT ('    V(Rmax)    = ',F8.3, ' ',4A)
99025 FORMAT ('    Veff(Rmax) = ',F8.3, ' ',4A)
99026 FORMAT ('    Gamma      =  ',F8.4,/
     &        '    Eta        =  ',F8.4)
99027 FORMAT ('    Moment * Gamma =',F8.3,' amu A**2')

99040 FORMAT (I5,2X,'General HindRotor:',2x,'Rotor symm. =',I2,
     & 34x,'zpe(cm-1) =',F8.1 )
99041 FORMAT (8X,A5,1x,' ; symmV =',I2,2x,'; Phase(rad.) = ',
     & F8.4,2x,'; Coeff. =',20(F10.4,1x))

99042 FORMAT (8X,A5,1x,' ; symmM =',I2,2x,'; Phase(rad.) = ',
     & F8.4,2x,'; Coeff. =',20(F10.4,1x))

99050 FORMAT (4x,'TUNNELING PARAMETERS: Unsymmetrical Eckart Barrier',/
     &        5x,'Imaginary Freq.:',F10.2,'i cm-1',/
     &        5x,'    V(forward) :',F10.2,' cm-1',/
     &        5x,'    V(reverse) :',F10.2,' cm-1'  )
99051 FORMAT (4x,'TUNNELING IGNORED: V(forward) & V(reverse) must be',
     &           ' >0.01 cm-1, and Imag. Freq. must be >20i cm-1'/)
99060 FORMAT (I5,2X,'Anharmonic coupled vibrations: qvib from BDENS',
     &       1x,'external file'/,7x,A150)
99061 FORMAT (I5,2X,
     &          'Semi-classical TST: Cum. Rxn. Prob. (CRP) from SCTST',
     &       1x,'external file'/,7x,A150)
99062 FORMAT (8x,'Egrain1=',f10.2,3x,'Emax2=',F10.2,3x,'Trials: ',A6)
99063 FORMAT (8x,'Egrain1=',f10.2,3x,'Emax2=',f10.2,3x,'Vf=',F10.2,
     &        3x,'Vr=',f10.2,3x,'Trials: ',A6,3x,'Correction: ',A5,1x)
99070 FORMAT (I5,2X,'Spin-Orbit+Rotor: State (spin mult. =',I2,
     &        '; Lambda =',I2,'), Aso =',F9.3,'; B =',F7.3,'; D =',
     &          1PE10.3,' cm-1; rot sym =',I2)
99224 FORMAT( 3x, A150 )

      END subroutine

