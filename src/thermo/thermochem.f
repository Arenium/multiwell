      SUBROUTINE thermochem 

!  THERMOCHEMISTRY CALCULATION
      include 'declare.inc'
      REAL(8) Qelectron, Qtranslat, Qvibrat, Qrotat, Qhindrotat, Qtotal    ! partition fxns
      REAL(8) Celectron, Ctranslat, Cvibrat, Crotat, Chindrotat            ! heat capacities (Cp)
      REAL(8) Selectron, Stranslat, Svibrat, Srotat, Shindrotat, Sextsym   ! entropies
      REAL(8) Helectron, Htranslat, Hvibrat, Hrotat, Hhindrotat            ! enthalpy fxns
      REAL(8) CStop, HStop, SStop, QStop , delAMU
      REAL(8) tuneck , qtun
      REAL(8) zz , delz
      REAL(8) cc , en , hh
      INTEGER(4) Nprnt
      DIMENSION cc(Nmax,Maxvibs) , en(Nmax,Maxvibs) , hh(Nmax,Maxvibs)     ! Cp, S, and H for up to Nmax species and Maxvibs DOF; used only at 298.15 K

      SAVE

         Dzpe = 0.0d+00
         DDelH = 0.0d+00
         DDelH298 = 0.0d+00
         delAMU = 0.0d+00
         Nflag = 0     ! no reaction

      DO i = 1 , Ns
            IF ( REPROD(i).EQ.REAC ) THEN
               DDelH = DDelH - DelH(i)          ! Delta H(0) for reaction
               Dzpe = Dzpe - zzpe(i)            ! Delta zpe for reaction
               delAMU = delAMU - AMU(i)         ! Delta molar mass for reaction
               Nflag = 1                        ! reaction
            ELSEIF ( REPROD(i).EQ.PROD .OR. REPROD(i).EQ.CTST ) THEN
               DDelH = DDelH + DelH(i)          ! Delta H(0) for reaction
               Dzpe = Dzpe + zzpe(i)            ! Delta zpe for reaction
               delAMU = delAMU + AMU(i)         ! Delta molar mass for reaction
               Nflag = 1                        ! reaction
            ENDIF
      END DO
 
      IF ( ABS( delAMU ) .GT. 1.0d-004 ) THEN
        WRITE(*,*) '*** FATAL: Empirical formulas do not match ***'
        STOP '   Check to see that reac/prod type are labeled'
      ENDIF
            
c
c      Thermodynamic functions at T
c
       zz   = 1.d+00/T(Temp)
       delz = Tfrac*zz
       DO  jj = 1 , 3                        ! Three temperatures in order to get T-dependence
         TT = 1.0d+00/( zz + (jj-2)*delz )
         Qreac = 1.0d+00
         Qprod = 1.0d+00
         qtun = 1.0d+00
         DO i = 1 , Ns        ! species
            Qelectron = 1.d+00          ! partition fxns
            Qvibrat = 1.d+00
            Qrotat = 1.d+00
            Qhindrotat = 1.d+00
            Qrso = 1.0d+00
            Qtotal = 1.d+00
            Celectron = 0.d+00          ! Heat capacities (Cp)
            Ctranslat = 0.d+00
            Cvibrat = 0.d+00
            Crotat = 0.d+00
            Chindrotat = 0.d+00
            Crsot = 0.0d+00
            Sextsym = 0.d+00            ! Entropies
            Selectron = 0.d+00
            Stranslat = 0.d+00
            Svibrat = 0.d+00
            Srotat = 0.d+00
            Shindrotat = 0.d+00
            Srsot = 0.0d+00
            Helectron = 0.d+00          ! enthalpy fxns
            Htranslat = 0.d+00
            Hvibrat = 0.d+00
            Hrotat = 0.d+00
            Hhindrotat = 0.d+00
            Hrsot = 0.0d+00

            S(i) = Rgas*log( DBLE(Sopt(i)) / DBLE(Sym(i)) )
            Sextsym = S(i)
            en( i , N(i)+1 ) = Sextsym

            CALL qelect(TT,i,Nele(i),Elev,gele,q,Cel,Sel,Hel)   ! electronic partition function
            Qelectron = q
            H(i)  = Rgas*Hel*TT        ! H(T)-H(0)
            Cp(i) = Rgas*Cel
            S(i)  = S(i) + Rgas*Sel
            Celectron = Rgas*Cel
            Selectron = Rgas*Sel
            Helectron = Rgas*Hel
            cc( i , N(i)+2 ) = Celectron
            en( i , N(i)+2 ) = Selectron
            hh( i , N(i)+2 ) = Helectron

            CALL qtrans( TT, AMU(i), Sunits, Ctrans, Strans, Htrans,
     &                 Qtranslat )
            q = q*Qtranslat
            Ctranslat = Rgas*Ctrans    
            Stranslat = Rgas*Strans
            Htranslat = Rgas*Htrans
            Cp(i) = Cp(i) + Ctranslat
            S(i)  = S(i) + Stranslat
            H(i)  = H(i) + Htranslat
            cc( i , N(i)+3 ) = Ctranslat
            en( i , N(i)+3 ) = Stranslat
            hh( i , N(i)+3 ) = Htranslat

            IF ( N(i).NE.0 ) THEN                                   ! for non-monatomic species
               DO j = 1 , N(i)          ! degrees of freedom
                  IF ( IDOF(i,j).EQ.VIB ) THEN                          ! VIBRATION
                     CALL qmorse(TT,WE(i,j),ANH(i,j),NG(i,j),
     &                    qq,Cvib,Svib,Hvib,zap)
                     q = q*qq
                     Qvibrat=Qvibrat*qq
                     Cp(i) = Cp(i) + Cvib*Rgas
                     H(i) = H(i) + Rgas*Hvib*TT        ! H(T)-H(0)
                     S(i) = S(i) + Rgas*Svib
                     Cvibrat = Cvibrat + Cvib*Rgas
                     Svibrat = Svibrat + Rgas*Svib
                     Hvibrat = Hvibrat + Rgas*Hvib*TT        ! H(T)-H(0)
                     cc( i , j ) = Cvib*Rgas
                     en( i , j ) = Rgas*Svib
                     hh( i , j ) = Rgas*Hvib*TT
                  ELSEIF ( IDOF(i,j).EQ.BOX ) THEN                      ! PARTICLE IN A BOX
                     CALL Qboxmod(TT,WE(i,j),NG(i,j),
     &                     qq,Cpib,Spib,Hpib,zap)
                     q = q*qq
                     Cp(i) = Cp(i) + Cpib*Rgas
                     S(i) = S(i) +   Rgas*Spib
                     H(i) = H(i) +   Rgas*Hpib*TT        ! H(T)-H(0)
                     cc( i , j ) = Cpib*Rgas
                     en( i , j ) = Rgas*Spib
                     hh( i , j ) = Rgas*Hpib*TT
                  ELSEIF ( IDOF(i,j).EQ.GOR ) THEN   ! I2d Gorin         HINDERED GORIN MODEL
                     WE(i,j)=I2dGorin(Temp)
                      if (WE(i,j).gt.10) then 
                                                  ! CLASSICAL ROTOR
                       BB = 16.85763D+00/WE(i,j)  ! rotational constant (cm-1)
                       dim = NG(i,j)              ! rotor dimension
                       ss = INT(ANH(i,j))         ! rotor symmetry number (not included in external symmetry)
                       CALL Classrot(TT,BB,dim,ss,qq,Crot,Srot,Hrot)
                       q = q*qq
                       Qrotat = Qrotat*qq
                       Cp(i) = Cp(i) + Crot*Rgas
                       S(i) = S(i) +   Srot*Rgas
                       H(i) = H(i) +   Hrot*Rgas*TT        ! H(T)-H(0)
                       cc( i , j ) = Crot*Rgas
                       en( i , j ) = Srot*Rgas
                       hh( i , j ) = Hrot*Rgas*TT
                      else
                                                  ! QUANTUM ROTOR
                       BB = 16.85763D+00/WE(i,j)  ! rotational constant (cm-1)
                       dim = NG(i,j)
                                ! rotor dimension
                       ss = INT(ANH(i,j))
                                ! rotor symmetry number (not included in external symmetry)
                       CALL Quantrot(TT,BB,dim,ss,qq,Croq,Sroq,Hroq)
                       q = q*qq
                       Qrotat = Qrotat*qq
                       Cp(i) = Cp(i) + Croq*Rgas
                       S(i) = S(i) +   Sroq*Rgas
                       H(i) = H(i) +   Hroq*Rgas*TT        ! H(T)-H(0)
                       cc( i , j ) = Croq*Rgas
                       en( i , j ) = Sroq*Rgas
                       hh( i , j ) = Hroq*Rgas*TT
                      endif
                  ELSEIF ( IDOF(i,j).EQ.FIT ) THEN                      ! HINDERED GORIN MODEL FITTING
                       if (WE(i,j)*x.gt.10) then     ! x = provisional gamma (from subroutine hindrance_fit)
                                                 ! CLASSICAL ROTOR
                       BB = 16.85763D+00/(WE(i,j)*x)  ! rotational constant (cm-1)
                       dim = NG(i,j)              ! rotor dimension
                       ss = INT(ANH(i,j))         ! rotor symmetry number (not included in external symmetry)
                       CALL Quantrot(TT,BB,dim,ss,qq,Crot,Srot,Hrot)
                       q = q*qq
                       Qrotat = Qrotat*qq
                       Cp(i) = Cp(i) + Crot*Rgas
                       S(i) = S(i) +   Srot*Rgas
                       H(i) = H(i) +   Hrot*Rgas*TT        ! H(T)-H(0)
                       cc( i , j ) = Crot*Rgas
                       en( i , j ) = Srot*Rgas
                       hh( i , j ) = Hrot*Rgas*TT 
                      else
                                                ! QUANTUM ROTOR
                       BB = 16.85763D+00/(WE(i,j)*x)  ! rotational constant (cm-1); x = provisional gamma (from subroutine hindrance_fit)
                       dim = NG(i,j)
                                ! rotor dimension
                       ss = INT(ANH(i,j))
                                ! rotor symmetry number (not included in external symmetry)
                       CALL quantrot(TT,BB,dim,ss,qq,Croq,Sroq,Hroq)
                       q = q*qq
                       Qrotat = Qrotat*qq
                       Cp(i) = Cp(i) + Croq*Rgas
                       S(i) = S(i) +   Sroq*Rgas
                       H(i) = H(i) +   Hroq*Rgas*TT        ! H(T)-H(0)
                       cc( i , j ) = Croq*Rgas
                       en( i , j ) = Sroq*Rgas
                       hh( i , j ) = Hroq*Rgas*TT
                      endif
 
                  ELSEIF ( IDOF(i,j).EQ.ROT ) THEN
                                                                        ! CLASSICAL FREE ROTOR
                     BB = 16.85763D+00/WE(i,j)  ! rotational constant (cm-1)
                     dim = NG(i,j)                   ! rotor dimension
                     ss = INT(ANH(i,j))              ! rotor symmetry number (not included in external symmetry)
                     CALL classrot(TT,BB,dim,ss,qq,Crot,Srot,Hrot)
                     Qrotat = Qrotat*qq
                     q = q*qq
                     Cp(i) = Cp(i) + Crot*Rgas
                     H(i) = H(i) + Hrot*Rgas*TT        ! H(T)-H(0)
                     S(i) = S(i) + Srot*Rgas
                     Crotat = Crotat + Crot*Rgas
                     Srotat = Srotat + Srot*Rgas
                     Hrotat = Hrotat + Hrot*Rgas*TT        ! H(T)-H(0)
                     cc( i , j ) = Crot*Rgas
                     en( i , j ) = Srot*Rgas
                     hh( i , j ) = Hrot*Rgas*TT 

                  ELSEIF ( IDOF(i,j).EQ.QRO ) THEN
                                                                        ! QUANTUM FREE ROTOR
                     BB = 16.85763D+00/WE(i,j)  ! rotational constant (cm-1)
                     dim = NG(i,j)
                                ! rotor dimension
                     ss = INT(ANH(i,j))
                                ! rotor symmetry number (not included in external symmetry)
                     CALL quantrot(TT,BB,dim,ss,qq,Croq,Sroq,Hroq)
                     q = q*qq
                     Qrotat = Qrotat*qq
                     Cp(i) = Cp(i) + Croq*Rgas
                     H(i) = H(i) + Hroq*Rgas*TT        ! H(T)-H(0)
                     S(i) = S(i) + Sroq*Rgas
                     Crotat = Crotat + Croq*Rgas
                     Srotat = Srotat + Sroq*Rgas
                     Hrotat = Hrotat + Hroq*Rgas*TT        ! H(T)-H(0)
                     cc( i , j ) = Croq*Rgas
                     en( i , j ) = Sroq*Rgas
                     hh( i , j ) = Hroq*Rgas*TT

                  ELSEIF ((IDOF(i,j).EQ.HIN).OR.(IDOF(i,j).EQ.HRD)) THEN  ! HINDERED INTERNAL ROTATION

                  DO LL=1, IMAX(i,j)
                        EVh(LL)=EV(i,j,LL)
                  ENDDO
                     CALL qhinder(TT,EVh,IMAX(i,j),NSIG(i,j),
     &                        qq, Chin,Shin,Hhin,zap)
                     q = q*qq
                     Qhindrotat = Qhindrotat*qq
                     Cp(i) = Cp(i) + Chin*Rgas
                     H(i) = H(i) + Hhin*Rgas*TT        ! H(T)-H(0)
                     S(i) = S(i) + Shin*Rgas
                     Chindrotat = Chindrotat + Chin*Rgas
                     Shindrotat = Shindrotat + Shin*Rgas
                     Hhindrotat = Hhindrotat + Hhin*Rgas*TT        ! H(T)-H(0)
                     cc( i , j ) = Chin*Rgas
                     en( i , j ) = Shin*Rgas
                     hh( i , j ) = Hhin*Rgas*TT

                  ELSEIF ( IDOF(i,j).EQ.RSO ) THEN             		! ROTATION + SPIN-ORBIT
                     IF ( Temp .EQ. Nt+1 .AND. jj .EQ. 2 ) THEN 
                      WRITE (23,*) '   '
                      WRITE (23,*) MOLNAME(i), 'levels used at 298.15 K'
                      Nprnt = 1
                      ELSE
                      Nprnt = 0
                     ENDIF
                     CALL qdrso(TT,ASO(i),BSO(i),D(i),LAMBDA(i),
     &                    MULT(i),NRSO(i),qq,Crso,Srso,Hrso,Nprnt )
                     q = q*qq
                     Qrso = qq
                     Cp(i) = Cp(i) + Crso*Rgas
                     H(i) = H(i) + Hrso*Rgas*TT        ! H(T)-H(0)
                     S(i) = S(i) + Srso*Rgas
                     Crsot = Crsot + Crso*Rgas
                     Srsot = Srsot + Srso*Rgas
                     Hrsot = Hrsot + Hrso*Rgas*TT        ! H(T)-H(0)
                     cc( i , j ) = Crso*Rgas
                     en( i , j ) = Srso*Rgas
                     hh( i , j ) = Hrso*Rgas*TT

                  ELSEIF ( IDOF(i,j).EQ.TOP ) THEN                      ! SYMMETRIC TOP (QUANTUM)
                     BJ = 16.85763D+00/WE(i,j)    ! 2D Bx rotational constant (cm-1)
                     BK = 16.85763D+00/ANH(i,j)   ! 1D Bz rotational constant (cm-1)
                     ss = NG(i,j)                 ! rotor symmetry number (not included in external symmetry)
                     CALL Qsymtop(TT,BJ,BK,ss,qq,CStop,SStop,HStop)
                     q = q*qq
                     Qrotat = Qrotat*qq
                     Cp(i) = Cp(i) + CStop*Rgas
                     H(i) = H(i) + HStop*Rgas*TT        ! H(T)-H(0)
                     S(i) = S(i) + SStop*Rgas
                     Crotat = Crotat + CStop*Rgas
                     Srotat = Srotat + SStop*Rgas
                     Hrotat = Hrotat + HStop*Rgas*TT        ! H(T)-H(0)
                     cc( i , j ) = CStop*Rgas
                     en( i , j ) = SStop*Rgas
                     hh( i , j ) = HStop*Rgas*TT
                     
                  ELSEIF ( IDOF(i,j).EQ.QVB ) THEN                       ! ANHARMONIC VIBRATIONAL PARTITION FUNCTION, external file
                     CALL qvbint( i , TT , Qxx , Cxx , Hxx , Sxx )
                     q = q*Qxx                           ! i=species, Temp=index for temperature
                     Qvibrat=Qvibrat*Qxx
                     Cp(i) = Cp(i) + Cxx*Rgas
                     H(i) = H(i) + Hxx*Rgas*TT           ! H(T)-H(0)
                     S(i) = S(i) + Sxx*Rgas
                     Cvibrat = Cvibrat + Cxx*Rgas
                     Svibrat = Svibrat + Sxx*Rgas
                     Hvibrat = Hvibrat + Hxx*Rgas*TT           ! H(T)-H(0)
                     cc( i , j ) = Cxx*Rgas
                     en( i , j ) = Sxx*Rgas
                     hh( i , j ) = Hxx*Rgas*TT

                  ELSEIF ( IDOF(i,j).EQ.CRP ) THEN                       ! SEMI-CLASSICAL TST, external file
                     CALL qvbint( i , TT , Qxx , Cxx , Hxx , Sxx )
                     q = q*Qxx                           ! i=species, Temp=index for temperature
                     Cp(i) = Cp(i) + Cxx*Rgas
                     S(i) = S(i)   + Sxx*Rgas
                     H(i) = H(i)   + Hxx*Rgas*TT           ! H(T)-H(0)
                     cc( i , j ) = Cxx*Rgas
                     en( i , j ) = Sxx*Rgas
                     hh( i , j ) = Hxx*Rgas*TT
                  ENDIF

               ENDDO    ! DOF loop
            ENDIF

            IF ( REPROD(i) .EQ. REAC ) THEN   ! note that q is actual true total partition fxn.
               Qreac = Qreac*q
            ELSEIF ( REPROD(i).EQ.PROD .OR. REPROD(i).EQ.CTST ) THEN
               Qprod = Qprod*q
            ENDIF

            IF( Temp.LT.(Nt+1) .AND. jj.EQ.2 ) THEN                            ! Summary of contributions at each temperature

              WRITE (21,*) i, MOLNAME(i), TT, Qelectron, Qtranslat,      
     &               Qvibrat, Qrotat, Qhindrotat, Qrso
c              WRITE (21,9903) i, MOLNAME(i), TT, Qelectron, Qtranslat,      
c     &               Qvibrat, Qrotat, Qhindrotat, Qrso
              WRITE (22,9904) i, MOLNAME(i), TT, Sextsym,
     &              Celectron, Selectron, Helectron/DIV,
     &              Ctranslat, Stranslat, Htranslat/DIV,     
     &              Cvibrat, Svibrat , Hvibrat/DIV,
     &              Crotat, Srotat, Hrotat/DIV,
     &              Chindrotat, Shindrotat, Hhindrotat/DIV,
     &              Crsot, Srsot, Hrsot/DIV
            ENDIF

            IF ( Temp.EQ.(Nt+1) .AND. jj.EQ.2 ) THEN                   ! Detailed contributions at 298.15 K
               write (22,9905) i, MOLNAME(i), Cp(i), S(i) , H(i)/DIV,
     &         cc( i, N(i)+1 ), en( i, N(i)+1 ), hh( i, N(i)+1 )/DIV ,   ! external symm and optical isomers
     &         cc( i, N(i)+2 ), en( i, N(i)+2 ), hh( i, N(i)+2 )/DIV ,   ! electronic
     &         Ctranslat, Stranslat, Htranslat/DIV,                      ! translation
     &         Cvibrat, Svibrat , Hvibrat/DIV,                           ! vibrations
     &         Crotat, Srotat, Hrotat/DIV,                               ! free rotors
     &         Chindrotat, Shindrotat, Hhindrotat/DIV,                   ! hindered rotors
     &         Crsot, Srsot, Hrsot/DIV							   ! rotation + Spin-Orbit

              DO k = 1, N(i)
                write (22,9906) k, IDOF(i,k), cc(i,k) , en(i,k), 
     &                              hh(i,k)/DIV
              END DO  ! k; DOF loop
            ENDIF ! 298.15 K detailed contributions output

         ENDDO ! i; species loop
         
         IF(jj .EQ. 2) WRITE (21,*) ' '
         IF(jj .EQ. 2) WRITE (22,*) ' '
 
         DSR = 0.0D+00
         DHR = Hdiff                            ! Hdiff is delH(rxn) at 0 K
         DCp = 0.0d+00
         DO i = 1 , Ns
            IF ( REPROD(i).EQ.REAC ) THEN
               DSR = DSR - S(i)                 ! DelS for reaction at Temp
               DHR = DHR - H(i)                 ! DelH for reaction at Temp  [ H(i) is the enthalpy fxn at T ]
               DCp = DCp - Cp(i)                ! DelCp for reaction at Temp
            ELSEIF ( REPROD(i).EQ.PROD.OR.REPROD(i).EQ.CTST ) THEN
               DSR = DSR + S(i)                 ! Delta S for reaction
               DHR = DHR + H(i)                 ! Delta H for reaction
               DCp = DCp + Cp(i)                ! Delta Cp for reaction
            ENDIF
         ENDDO
 
         DG = (DSR - DHR/TT)/Rgas

         AK(jj) = exp(DG)                          ! Equilibrium constant

         IF ( Nctst .EQ. 1 ) THEN                  ! Canonical TST rate constant; Vf, Vr, and vimag in cm-1 units
           IF ( (Vf .GT. 0.01).AND.(vimag .GT. 20.).AND.
     &        (Vr .GT. 0.01) ) THEN
              qtun = tuneck( TT, vimag, Vf, Vr )   ! tunneling transmission coefficient
           ELSE
              qtun = 1.0d+00
           ENDIF
           AK(jj) = qtun*AK(jj)*2.0837d+10*TT      ! * kT/h
         ENDIF
         
       ENDDO ! end of jj finite difference temperature loop
 
      BB = LOG( AK(3)/AK(1) ) / ( 2.0d+00*delz )                   ! exponential factor
      AA = AK(2)*exp( -BB / T(Temp) )                              ! pre-exponential factor
 
9903  FORMAT ( 3x, I1, 1x, A32, 1x, F10.2, 6(1x, 1pe12.4)  )
9904  FORMAT ( 3x, I1, 1x, A32, 1x, F10.2, 19(1x, f10.3)  )
9905  FORMAT ( //'-------------------------------------------------',
     & '------------------------------------------------------',//,
     & 'DETAILED INFORMATION (298.15 K):      (',I1,') ',A32,//,
     &  35x,'Cp',14x,'S', 7x,'H(298.15) - H(0)',/,
     & 'TOTAL', 26x, f10.3, 2(5x, f10.3),/,
     & '(Ext)Symm & Opt.Isom Term      ',f10.3, 2(5x , f10.3),/,
     & 'Electronic                     ',f10.3, 2(5x , f10.3),/,
     & 'Translation                    ',f10.3, 2(5x , f10.3),/,
     & 'Vibrations                     ',f10.3, 2(5x , f10.3),/,
     & 'Free Internal & External Rotors',f10.3, 2(5x , f10.3),/,
     & 'Hindered Internal Rotors       ',f10.3, 2(5x , f10.3),/,
     & 'Rotation + Spin-Orbit          ',f10.3, 2(5x , f10.3),//,
     & 'Individual D.O.F. in input file',4x,'Cp',14x,'S', 7x,
     & 'H(298.15) - H(0)')
9906  FORMAT (10x, I3, 5x, A3, 10x, f10.3, 2(5x , f10.3) )

      end subroutine 
