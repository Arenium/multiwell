
      SUBROUTINE GHRLEV(EV,NN,Nmax,Emax,B,NCV,NCF,CV,CF,
     &            Nsiv,Nsif,Vhr,Bhr,Phav,Phab,Vo,Vm)
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c GHRLEV: a code of general purposes to calculate eigenvalues of 1D-hindered internal rotation 
c Copyright (C) 2009 Lam T. Nguyen
c
c      DATE:      Feb. 13, 2009
c
c Lam T. Nguyen
c NGUYENLT@umich.edu
c University of Michigan
c Ann Arbor, MI 48109-2143
c
c or contact
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
c See the "ReadMe" file for a copy of the GNU General Public License,
c or contact:
c
c Free Software Foundation, Inc.
c 59 Temple Place - Suite 330
c Boston, MA 02111-1307, USA.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C      INPUT:      
C            1) Emax            : maximum energy in Master Equation
C            2) B            : average rotational constant, defined as:
C                  B = (Integral of B(x)dx from 0 to 2*PI) / (2*PI)
C            3) NCV, NCF      : number of elements in vectors of CV and CF, respectively
C            4) CV            : vector of coefficients with NCV elements in torsional 
C                  potential energy function, which is expressed as = CVo/2 + Sum of CVi*(1-cos(i*x*Nsig))/2
C            5) CF            : vector of coefficients with NCF elements in rotational constant,
C                  which is given by a Fourier series = CFo + Sum of CFi*cos(i*x*Nsig)
C            6) Nsiv, Nsif      : symmetry numbers for tortional potential E function and rotational constant, respectively.
C            7) NN            : maximum number of elements in vector EV
C            8) Vhr            : model for torsional potential energy function
C            9) Bhr            : model for rotational constant
C            10) Phav, Phab      : phases of torsional angles in Vhr and Bhr, respectively.         
C
C      OUTPUT:      1) EV(Nmax)      : vector of eigenvalues with Nmax elements
C            2) Nmax            : number of elements of vector EV
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      IMPLICIT NONE
      REAL(8) EV, Phav, Phab, Vo, Vm, CV, CF, Emax,b,Vtp
      INTEGER(4) NN,NCV,NCF,Nsiv,Nsif,NX,Nmax
      PARAMETER (NX=501)
C
C.....Declare arrays.
      DIMENSION EV(NN), CV(NCV), CF(NCF)
      CHARACTER(len=5) Vhr, Bhr
      SAVE

      IF((Vhr.EQ.'VHRD2').OR.(Vhr.EQ.'Vhrd2').OR.
     &      (Vhr.EQ.'VHRD3').OR.(Vhr.EQ.'Vhrd3')) THEN
        CALL CALVo(Vo,Vm,CV,NCV,Nsiv,Vhr,Phav)
      CALL ODQHR(EV,NN,NCV,NCF,CV,CF,Nsiv,Nsif,Vhr,Bhr,Phav,Phab,Vm)
      ELSE
        CALL CALVo(Vo,Vm,CV,NCV,Nsiv,Vhr,Phav)
      Vm=0.0d0
        CALL ODQHR(EV,NN,NCV,NCF,CV,CF,Nsiv,Nsif,Vhr,Bhr,Phav,Phab,Vm)
      ENDIF
      Nmax=NX

      Vtp=EV(Nmax)

!      write(*,999) "Maximum energy (cm-1) = ", Emax
!      write(*,999) "The highest eigenvalue (cm-1) = ", Vtp
!      write(*,999) "Vo = Vmax (cm-1) = ", Vo
!      write(*,999) "Vmin (cm-1) = ", Vm
!999      FORMAT(1X,A33,1X,F10.1)

      IF(Vtp.LT.Vo) THEN
            write(*,*) "EIGENVALUES from HinRotor are not reliable !!!"
            write(*,*) "You have to increase number of grid points"
            write(*,*) "OR Harmonic approximation is good enough"
            STOP
      ENDIF
 
      IF(Emax.LT.Vo*10) Emax=Vo*10 
      CALL CALEI(EV,NN,Nmax,Emax,B)

      RETURN 
      END


C****************************************************************************************************************
C      NAME:      Subroutine ODQHR 
C      USE:      To obtain a vector of quantum eigenvalues of 1D-hindered internal rotor, in which  
C            both torsional potential energy function and moment inertia are explicitly treated as 
C            functions of internal rotational angle, i.e. non-rigid rotor, based on Meyer's algorithm:
C            R. Meyer J. Chem. Phys. 52 (1970) 2053-2059.  
C
C      AUTHORS:      Lam T. Nguyen and John R. Barker 
C
C      DATE:      Jan. 26, 2009
C
C****************************************************************************************************************
C     
C       NX  : Number of grid points
C        XA  : Position on a grid.
C        AR  : Hamiltonian Matrix, i.e. Meyer's symmetrical matrix
C        ER  : Vector of eigenvalues 
C        ZR  : Eigenvectors (X,Y); where
C                                      X : Wavefunction.
C                                      Y : Energy level.       
C        NPR      : 0 for eigenvalue only and 1 for both eigenvalue and eigenvector 
C        RMIN     : Starting point of grid
C        RMAX     : End point of grid
C        ZL       : Grid length
C        DX       : Grid spacings
C
C****************************************************************************************************************
C
C      INPUT:      1) NCV, NCF      = Number of coefficients derived from fitting rotational energy profile and moment inertia
C            2) CV(NCV)      = Vector of coefficients from rotational energy profile, with NCV elements
C            3) CF(NCF)      = Vector of coefficients from moment inertia (rotational constant), with NCF elements
C            4) Nsig            = rotational symmetry number
C      OUTPUT:      1) ER             = Vector of eigenvalues of 1D-hindered internal rotor, with NX=501 elements 
C
C*****************************************************************************************************************

      SUBROUTINE ODQHR(ER,NN,NCV,NCF,CV,CF,Nsiv,Nsif,Vhr,Bhr,
     &            Phav,Phab,Vm)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
      PARAMETER (NX=501)
C
C.....Declare arrays.
      DIMENSION AR(NX,NX),ER(NN),ZR(NX,NX),XA(NX) 
      DIMENSION D(NX,NX), F(NX), V(NX), DT(NX,NX),
     & T(NX,NX), FT(NX,NX)
      DIMENSION  CF(NCF), CV(NCV), FV1(NX), FV2(NX)
        CHARACTER Vhr*5, Bhr*5

C.....      DIMENSION WORK(NX*3)

C
C.....Variable data input.

      DATA NPR/0/

C.....Now compute Hamiltonian matrix:

      CALL CMM(AR,NCV,NCF,CV,CF,Nsiv,Nsif,Vhr,Bhr,Phav,Phab,Vm)
                              
C.....Now call eigenvalue solver using RS subroutine of EISPACK.

      CALL rs(NX,NX,AR,ER,NPR,ZR,FV1,FV2,IERR) 

C...... Alterative option: using subroutine DSYEV from LAPACK 
C......      CALL DSYEV( 'V', 'U', NX, AR, NX, ER, WORK, NX*3, IERR )
C......

      RETURN 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C.....      CMM stands for Constructing Meyer's Matrix, AR, for 1D-hindered internal rotor using Meyer's method, 
C.....      R. MEYER J. Chem. Phys. 52 (1970) 2053.
C.....      Both rotational energy profile and moment inertia are properly treated as functions of 
C.....      internal rotation angle, i.e. non-rigid rotor
C.....      INPUT:      1) NCF, NCV = Number of coefficients drived from fitting rotational energy profile and moment inertia
C.....            2) CV(NCV) = Vector of coefficients from rotational energy profile, with NCV elements
C.....            3) CF(NCF) = Vector of coefficients from moment inertia, with NCF elements
C.....            4) Nsig = Rotational symmetry number, it can be either one or two or three or whatever 
C.....      OUTPUT: 1) AR(NX,NX) = Meyer's symmetrical matrix, with a size of NX--number of grid points   


      SUBROUTINE CMM(AR,NCV,NCF,CV,CF,Nsiv,Nsif,Vhr,Bhr,Phav,Phab,Vm)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
      PARAMETER (NX=501)
C
C.....Declare arrays.
      DIMENSION AR(NX,NX),XA(NX) 
      DIMENSION D(NX,NX), F(NX), V(NX), DT(NX,NX),
     & T(NX,NX), FT(NX,NX)
      DIMENSION  CF(NCF), CV(NCV)
      CHARACTER Vhr*5, Bhr*5

C
C.....Variable data input.

      PI=acos(-1.D0)
C
C.....Set up grid
      RMIN=0.0D0
      RMAX=PI*2
      ZL=(RMAX-RMIN)
      DX=ZL/NX

C.....Now compute Meyer's symmetrical matrix:

      DO I=1, NX
            XA(I)=RMIN+DX*I
            CALL VX(V(I),XA(I)+Phav,CV,NCV,Nsiv,Vhr,Vm)
            CALL FX(F(I),XA(I)+Phab,CF,NCF,Nsif,Bhr)
            DO J=1, NX
                  ImJ=I-J
                  IF(ImJ.EQ.0) THEN
                        D(I,J)=0.0D0
                  ELSE
                        tp1=PI*ImJ/dble(NX)
                        tp2=(-1.0D0)**ImJ
                        D(I,J)=tp2/2.0D0/dsin(tp1)
                  ENDIF
                  DT(I,J)=-D(I,J)
                  FT(I,J)=F(I)*DELTA(I,J)
            ENDDO
      ENDDO

      DO I=1, NX
            DO J=1, NX
                  SU=0.0D0
                  DO K=1, NX
                        SU=SU + FT(I,K)*D(K,J)
                  ENDDO
                  T(I,J)=SU
            ENDDO
      ENDDO

      DO I=1, NX
            DO J=1, NX
                  SU=V(I)*DELTA(I,J)
                  DO K=1, NX
                        SU=SU + DT(I,K)*T(K,J)
                  ENDDO
                  AR(I,J)=SU
            ENDDO
      ENDDO      
      RETURN 
      END
C..........
C..........
      SUBROUTINE VX(V,X,CV,N,Nsig,Vhr,Vm)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
      DIMENSION CV(N)
      CHARACTER Vhr*5
      DOUBLE PRECISION Vm

      IF((Vhr.EQ.'Vhrd1').OR.(Vhr.EQ.'VHRD1')) THEN
            V=0.0d0
            DO I=1, N
                  V=V+CV(I)*(1.0d0-dcos(X*Nsig*I))/2.0d0
            ENDDO
      ELSEIF((Vhr.EQ.'Vhrd2').OR.(Vhr.EQ.'VHRD2')) THEN
            V=CV(1)
            DO I=2, N
                  V=V+CV(I)*dcos(X*Nsig*(I-1))
            ENDDO
        ELSEIF((Vhr.EQ.'Vhrd3').OR.(Vhr.EQ.'VHRD3')) THEN
            V=CV(1)
            Nc=(N+1)/2
            DO I=2, N
                  IF(I.LE.Nc) THEN
                        V=V+CV(I)*dcos(X*Nsig*(I-1))
                  ELSE
                        V=V+CV(I)*dsin(X*Nsig*(I-Nc))
                  ENDIF
            ENDDO
      ELSE
            write(*,*) "ERROR at INPUT for torsional function = Vhrd"
            STOP
      ENDIF

      V=V-Vm

      RETURN
      END
C..........
C..........
      SUBROUTINE FX(F,X,CF,N,Nsig,Bhr)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
      DIMENSION CF(N)
      PARAMETER (CON=16.85763D0)
      CHARACTER Bhr*5

      IF((Bhr.EQ.'Bhrd1').OR.(Bhr.EQ.'BHRD1')) THEN
            F=CF(1)
            DO I=2, N
                  F=F+CF(I)*dcos(X*Nsig*(I-1))
            ENDDO
        ELSEIF((Bhr.EQ.'Ihrd1').OR.(Bhr.EQ.'IHRD1')) THEN
                F=CF(1)
                DO I=2, N
                        F=F+CF(I)*dcos(X*Nsig*(I-1))
                ENDDO
            F=CON/F
      ELSE
                write(*,*) "ERROR at INPUT for Bhrd or Ihrd"
            STOP
      ENDIF

      RETURN
      END
C.........
C.........      
      Double Precision Function DELTA(I,J)
      INTEGER I, J

      IF(I.EQ.J) THEN
            DELTA=1.0D0
      ELSE
            DELTA=0.0D0
      ENDIF
      
      RETURN
      END
C.........
C.........

CCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CALEI(EV,NN,Nmax,Emax,B)
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c CALEI: a code to compute/select eigenvalues of torsional motion based on Emax 
c Copyright (C) 2009 Lam T. Nguyen and John R. Barker
c
c John R. Barker
c jrbarker@umich.edu
c University of Michigan
c Ann Arbor, MI 48109-2143
c (734) 763 6239
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C      INPUT:      1) EV(Nmax)      : a vector of eigenvalues with Nmax elements
C            2) Nmax            : number of elements of vector EV
C            3) Emax            : maximum energy in Master Equation
C            4) B            : average rotational constant, defined as:
C                  B = (Integral of B(x)dx from 0 to 2*PI) / (2*PI)
C      OUTPUT:      1) EV(Nmax)
C            2) Nmax      
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
      DIMENSION EV(NN)

      I=Nmax

C.......IF The largest eigenvalue is larger than the maximum energy, 
C.......THEN searching for an eigenvalue, which is <= Emax
11      E=EV(I)
      IF(E.GT.Emax) THEN
            I=I-1
            goto 11
      ENDIF
      Nmax=I

C.......Looking for higher-lying eigenvalues with double degeneracy
c111      NE1=NINT(EV(I))
c      NE2=NINT(EV(I-1))
c      IF(NE1.GT.NE2) THEN
c            I=I-1
c            goto 111
c      ENDIF

C.......Switching from a hindered rotor to a Pitzer rotor
C.......Pitzer rotor is defined as E = Eo + B*(JP**2)
C.......Two variables Eo and JP will be computed as below
c      tp=EV(I-1) - EV(I-2)
c      JP=NINT((tp/B - 1.0d0)/2.0d0)
c      Eo=EV(I-2)-B*JP*JP

c1111      EV(I+1)=Eo + B*(JP+2)*(JP+2)
c      EV(I+2)=EV(I+1)
c      IF(EV(I+2).LE.Emax) THEN
c            I=I+2
c            JP=JP+1
c            goto 1111
c      ENDIF
c      Nmax=I

      RETURN
 
      END

CCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CALVo(Vo,Vm,CV,NCV,Nsiv,Vhr,Pha)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
      DIMENSION CV(NCV)
      CHARACTER Vhr*5

      PI=acos(-1.0d0)

        RMIN=0.0D0
        RMAX=PI*2
        ZL=(RMAX-RMIN)
        DX=ZL/3600.0d0
        CALL VX(Vo,0.0d0+Pha,CV,NCV,Nsiv,Vhr,0.0d0)
      Vm=Vo
        DO I=1, 3600 
                XA = RMIN + DX*I
                CALL VX(Vtp,XA+Pha,CV,NCV,Nsiv,Vhr,0.0d0) 
            IF(Vtp.GT.Vo) Vo=Vtp
                IF(Vtp.LT.Vm) Vm=Vtp
      ENDDO
      RETURN
      END

CCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC

