      IMPLICIT NONE
      CHARACTER(len=7) AVERSION
      CHARACTER(len=8) ADATE
      PARAMETER ( AVERSION='2023', ADATE='Mar 2023' )

      INTEGER(4) Nmax , lenstr, a1, Temp , Maxvibs
      PARAMETER (Nmax=6)                                      ! maximum number of species
      PARAMETER (Maxvibs=1000)                                ! maximum number of internal DOF

      INTEGER(4) NCV, NCB, IMAX, LLL, LL, LK , NSV, NSB
      DIMENSION NCV(Nmax,Maxvibs), NCB(Nmax,Maxvibs), 
     &            IMAX(Nmax,Maxvibs),
     &            NSV(Nmax,Maxvibs), NSB(Nmax,Maxvibs)
      REAL(8) CV , CB ,  CBB , CVV , EV , EVh
      DIMENSION CV(Nmax,Maxvibs,20) , CB(Nmax,Maxvibs,20) , 
     &            CBB(20) , CVV(20) , EV(Nmax,Maxvibs,2000) , EVh(2000) 
      REAL(8) Emax , Phav , Phab
      DIMENSION Phav(Nmax,Maxvibs) , Phab(Nmax,Maxvibs)
      CHARACTER(len=3)   IDOF , VIB , ROT , QRO , RSO
      CHARACTER(len=4)   REPROD , REAC , PROD
      CHARACTER(len=32)  MOLNAME
      CHARACTER(len=33)  CAP
      CHARACTER(len=4)   Eunits , Sunits
      CHARACTER(len=12)  UNITS
      CHARACTER(len=3)   HIN , HRA , HRB , HRC , HRD , VHAR , 
     &                   OBS , HAR , BOX , QVB , ECK , CRP, SCTST 
      CHARACTER(len=4)   CTST , NON
      CHARACTER(len=150) TITLELINE , TITLE4 , TITLE5  
      CHARACTER(len=150) DUMMY
      CHARACTER(len=3)   GOR , FIT , TOP
      CHARACTER(len=5)   VHRD(Nmax,Maxvibs) , BHRD(Nmax,Maxvibs)

      CHARACTER(len=6) ATYPE , KEYWORD2
      DIMENSION ATYPE(Nmax,10), KEYWORD2(Nmax)
     
      CHARACTER(len=5) VPTx
      DIMENSION VPTx(Nmax)

      INTEGER(4) N(Nmax) , Ns , NG , MODE , i , Sym , Nele , Sopt , 
     &           Nreac , Nprod , II , j , Nt , k , ss , dim , NSIG, 
     &           gele , ne , Nflag , Nctst , Nnone , jj, hindrance, 
     &           gorin, Warning(6), vts , Nelement(Nmax) , 
     &           Natom(Nmax,10) , nlines , line , MULT , LAMBDA, nso ,
     &           NRSO

      REAL(8) WE , ANH , AMU , DelH , Cp , Sfac , ff ,T, Tm,
     &        Hvib , Hrot , Hroq , Svib , Srot , Sroq , S , H ,
     &        Rgas , B , AMOM , W , VV , Qhin , Hhin , Chin ,
     &        Shin , Croq , zpe , zap , zzpe, DDelH, DG , Qele,
     &        Elev , Sel , Hel , Cel , TT , TK , WOBS , 
     &        Hfxn298 , gamma(100) , delH298 , DDelH298 ,
     &        delG298 , DDelG298 , Sfxn , Sfel, Vmin, Vmax ,
     &        zpp , QT , ASO , BSO , D

      DIMENSION T(200), Tm(100),
     &      IDOF(Nmax,Maxvibs) , WE(Nmax,Maxvibs) , ANH(Nmax,Maxvibs) ,
     &      NG(Nmax,Maxvibs) , MODE(Nmax,Maxvibs) , REPROD(Nmax) ,
     &      NSIG(Nmax,Maxvibs), TITLELINE(Nmax,20) , nlines(Nmax) ,
     &      TITLE4(Nmax,Maxvibs), TITLE5(Nmax,Maxvibs) ,
     &      WOBS(Nmax,Maxvibs), A(Nmax,Maxvibs) ,
     &      Hfxn298(Nmax), delH298(Nmax), delG298(Nmax) , 
     &      Sfel(Nmax) , zpe(Nmax,Maxvibs) , QT(Nmax,200) , zpp(Nmax) ,
     &		  ASO(Nmax) , MULT(Nmax) , D(Nmax) , LAMBDA(Nmax) ,
     &      BSO(Nmax) ,  NRSO(Nmax) , nso(Nmax) 

      DIMENSION AMOM(Nmax,Maxvibs) , W(Nmax,Maxvibs) , 
     &          VV(Nmax,Maxvibs) ,
     &          B(Nmax,Maxvibs), Vmin(Nmax,Maxvibs), Vmax(Nmax,Maxvibs)

      DIMENSION AMU(Nmax) , MOLNAME(Nmax) , Sym(Nmax) ,
     &          Sopt(Nmax) , DelH(Nmax) , H(Nmax) , S(Nmax) ,
     &          Cp(Nmax) , zzpe(Nmax)

      DIMENSION  Nele(Nmax) , gele(Nmax,10) , Elev(Nmax,10)

      REAL(8) Ctrans, Strans, Htrans
      
      REAL(8) Qvib , Qbox , Qrot , q , Qreac , Qprod , Kequil ,
     &                 Hdiff , Cvib , Crot , qq , AA , BB , DSR , DHR ,
     &                 Qroq ,  DCp , Dzpe, Hpib, Cpib, Spib , AK ,
     &                 Tdel , Tfrac, BJ, BK , DIV ,
     &        Qrso , Crso , Hrso , Srso , Crsot , Hrsot , Srsot
     
      DIMENSION AK(3)

      PARAMETER (Tfrac=1.d-04)
      
      PARAMETER (VIB='VIB',ROT='ROT',QRO='QRO',REAC='REAC',PROD='PROD',
     &           CAP='Qelectr  Entropy  Cp  [H(T)-H(0)]', HIN='HIN',
     &           HRA='HRA',HRB='HRB',HRC='HRC',HRD='HRD',OBS='OBS',
     &           HAR='HAR',BOX='BOX', ECK='ECK', TOP='TOP', QVB='QVB', 
     &           CRP='CRP',NON='NONE', CTST='CTST', RSO='RSO',GOR='GOR', 
     &           FIT='FIT' )

      Character(len=7)  POT
      CHARACTER(len=15) TypePOT
      
      REAL(8) freq, I2d, De, mu,beta, re,rx,derVeff,der2Veff, c
      REAL(8) I2dGorin(100),rmax(100),Veff(100),V(100),bet(100)
      REAL(8) Kexp(100),count,Kcin,Kcinold
      REAL(8) I2d_Gorin, Weorig(10)
      REAL(8) a,aold,m, x,xold,xold2,start, Deorig
      REAL(8) T_AA(200),T_BB(200),T_DSR(200)
      REAL(8) T_DHR(200),T_DCp(200), T_AK(200)
      REAL(8) T_Sfac(200),T_S(200,6),T_Cp(200,6),T_H(200,6),V0
      REAL(8) vimag , Vf , Vr
      
      REAL(8) Egrain1 , Emax1, Emax2 , viblo
      DIMENSION Egrain1(Nmax,Maxvibs) , Emax1(Nmax,Maxvibs) ,
     &          Emax2(Nmax,Maxvibs)
      
c     External temperatures and variables from SCTST and ADENSUM; used for interpolation in THERMO
c
      INTEGER(4) maxNtx , Ntx
      PARAMETER (maxNtx=150)                                       ! 150 temperatures passed from BDENS or SCTST
      REAL(8) Tx, Cx , Hx , Sx , Qxx , Hxx , Cxx, Sxx
      DIMENSION Tx(maxNtx) , Cx(Nmax,maxNtx) , Hx(Nmax,maxNtx), 
     &          Sx(Nmax,maxNtx)
      
      INTEGER(4) Imax1, Imax2, NCRP, NQVB
      DIMENSION NQVB(Nmax)

!     Real numbers

      COMMON  DIV, WE, ANH, AMU , DelH, Cp, Sfac, T, Tm, 
     & Hvib, Hrot, Hroq, Svib, Srot, Sroq, S, H,
     & Rgas, B, AMOM, W, VV , Hhin , Chin ,
     & Shin , Croq , zap , zzpe , zpe, DDelH, DG,
     & Elev, Sel , Hel, Cel, TK, WOBS, A, gamma,
     & q , Qreac, Qprod, Kequil, Hfxn298, DDelH298 , delH298 ,
     & delG298 , DDelG298 , Sfel , 
     & Hdiff ,Cvib ,Crot, qq ,AA ,BB , 
     & DSR, DHR , DCp, Dzpe, Hpib, Cpib, Spib , AK,
     & freq, I2d, De,mu, beta, re, rx , c,
     & I2dGorin, rmax , Veff , V, bet,
     & Kexp , Kcin , I2d_Gorin, WEorig, x, Deorig,
     & T_AA , T_BB , T_DSR , T_DHR , T_DCp, T_AK,
     & T_Sfac , T_S , T_Cp , T_H , V0, IMAX, EV, CV, CB, Emax,
     & Phav, Phab, CBB, CVV, Vmin, Vmax ,
     & vimag, Vf, Vr, QT , Egrain1 , Emax1, Emax2 , viblo ,
     & Tx , Cx , Hx , Sx, 
     & ASO, BSO, D, Qrso, Crso, Hrso, Srso, Crsot, Hrsot, Srsot
  
!    Integers/int/

       COMMON N, Ns , NG, MODE, Sym, Nele, Sopt, Nt, ss, NSIG, gele, ne,
     & Nflag, Nctst, hindrance, gorin, Warning, vts, Temp, Nreac,
     & Natom , Nelement, NCV, NCB, NSV, NSB, Imax1, Imax2, NCRP, NQVB , 
     & Ntx , nlines , MULT , LAMBDA, NRSO , nso
       
!    Character constants

       COMMON/char/ TITLELINE , TITLE4 , TITLE5 ,
     & REPROD , MOLNAME , SCTST , KEYWORD2 ,
     & ATYPE , IDOF , Eunits , Sunits , UNITS , VHAR , POT , 
     & TypePOT, VHRD, BHRD, VPTx
     
 
       

