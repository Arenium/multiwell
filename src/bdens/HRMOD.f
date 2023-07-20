      MODULE HRMOD
!
!    Module for sharing data between HINDERED ROTOR subroutines
!
c Copyright (C) 2017 John R. Barker
c
c John R. Barker
c University of Michigan
c Ann Arbor, MI 48109-2143
c (734) 763 6239
c jrbarker@umich.edu
c
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License (version 2)
c as published by the Free Software Foundation.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
c GNU General Public License for more details:
c
c     Free Software Foundation, Inc.
c     59 Temple Place - Suite 330
c     Boston, MA 02111-1307, USA.
c
      IMPLICIT NONE
      SAVE
 
      INTEGER, PARAMETER :: Levels = 20001          ! maximum no. of computed energy levels 
      
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      CONTAINS

      SUBROUTINE SHRLEV(T,DELE,JMAX,B,V,Nsiv,Nsif,IMAX,zpe,
     &          Vhr,Bhr,Phav,Phab)
c SHRLEV: a code to calculate eigenvalues of symmetrical, rigid 1D-hindered internal rotation
c Copyright (C) 2009 Lam T. Nguyen and John R. Barker
c
c	DATE:	Feb. 13, 2009
c
c John R. Barker
c jrbarker@umich.edu
c University of Michigan
c Ann Arbor, MI 48109-2143
c (734) 763 6239

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JMAX                      ! number of energy bins
      REAL(8), DIMENSION(JMAX), INTENT(INOUT) :: T     ! Beyer-Swinehart array for density of states
      REAL(8), INTENT(IN) :: DELE                      ! energy grain size
      REAL(8), INTENT(IN) :: B                         ! average rotational constant, defined as: B = (Integral of B(x)dx from 0 to 2*PI) / (2*PI)
      REAL(8), INTENT(IN) :: V                         ! hindrance barrier height (from minimum)
      INTEGER, INTENT(IN) :: Nsiv                      ! symmetry number for torsional potential function
      INTEGER, INTENT(IN) :: Nsif                      ! symmetry number for torsional rotational constant (Bhr)
      INTEGER, INTENT(IN) :: IMAX                      ! grid points and maximum no. of elements in EV array
      REAL(8), INTENT(OUT) :: zpe                      ! zero point energy
      CHARACTER(LEN=5), INTENT(IN) :: Vhr              ! model for torsional potential energy function (Vhr)
      CHARACTER(LEN=5), INTENT(IN) :: Bhr              ! model for rotational constant
      REAL(8), INTENT(IN) :: Phav                      ! phase of torsional angle in Vhr
      REAL(8), INTENT(IN) :: Phab                      ! phase of torsional angle in Bhr
      INTEGER :: NN                               ! number of elements in vector EV where EV(NN) <= Emax
      REAL(8) :: Emax                             ! maximum energy in Master Equation (corresponding to JMAX and DELE)
      REAL(8), DIMENSION(IMAX+1) :: EV            ! eigenvalues with Imax elements
      INTEGER :: NCV                              ! number of elements in CV
      INTEGER, PARAMETER :: No=1                  ! index number
      REAL(8), DIMENSION(No) :: CV                ! coefficients for Vhr fourier series 
      INTEGER :: NCF                              ! number of elements in CF
      REAL(8), DIMENSION(No) :: CF                ! coefficients for rotational constant fourier series
      REAL(8), DIMENSION(JMAX) :: AT              ! Beyer-Swinehart array for density of states
	 INTEGER, DIMENSION(IMAX) :: IR              ! array of bin numbers for energy levels
	 INTEGER :: I , J, L, LL                     ! index numbers
	 REAL(8) :: tp                               ! temporary variable
      
        Emax=(JMAX-1)*DELE + 1000.              ! *****raise max energy to account for later subtraction of zpe******

        CV(1)=V
        CF(1)=B

        CALL GHRLEV(EV,NN,IMAX,Emax,B,No,No,CV,CF,Nsiv,Nsif,
     &          Vhr,Bhr,Phav,Phab)

        zpe=EV(1)

c*********************************************************************
c     The following lines are for output of hindered rotor eigenvalues
c*********************************************************************
        write(12,*) 'Eigenvalues for Separable 1-D Symm Hindered Rotor'
        write(12,*) '          [shrlev]'
        write(12,*) '    i     E(i)            E(i)-zpe'

        DO I=1, IMAX
          tp=EV(I)-zpe
          write(12, 99) I, EV(I), tp
        ENDDO
99    FORMAT (3x, I5, 1x, ES12.5 , 5x, ES12.5 )
c*********************************************************************

      DO I=1, IMAX
		tp=EV(I)-zpe
		IR(I)=NINT( tp/DELE )
      ENDDO 

      DO J=1, JMAX                 ! initialize AT(i)
		AT(J)=0.0d0
      ENDDO 

      DO J=1, IMAX
		LL=IR(J)
		DO L=1, JMAX - LL
			AT(L+LL) = AT(L+LL) + T(L)
		ENDDO
      ENDDO
      DO J=1, JMAX
		T(J)=AT(J)/Nsiv
      ENDDO 

      RETURN
      END SUBROUTINE SHRLEV

!-----------------------------------------------------------------------
        SUBROUTINE UHRLEV(T,DELE,JMAX,NCB,NCV,CB,CV,Nsiv,Nsif,IMAX,
     &          zpe,Vhr,Bhr,Phav,Phab,Nsig)
c UHRLEV: a code to calculate eigenvalues of unsymmetrical, non-rigid 1D-hindered internal rotation
c Copyright (C) 2009 Lam T. Nguyen and John R. Barker
c
c	DATE:	Feb. 13, 2009
c
c Contact:
c John R. Barker
c jrbarker@umich.edu
c University of Michigan
c Ann Arbor, MI 48109-2143
c (734) 763 6239
c
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JMAX                      ! number of energy bins
      REAL(8), DIMENSION(JMAX), INTENT(INOUT) :: T     ! Beyer-Swinehart array for density of states
      REAL(8), INTENT(IN) :: DELE                      ! energy grain size
      INTEGER, INTENT(IN) :: NCB                       ! number of elements in CF
      INTEGER, INTENT(IN) :: NCV                       ! number of elements in CV
      REAL(8), DIMENSION(NCB), INTENT(IN) :: CB        ! coefficients for rotational constant fourier series
      REAL(8), DIMENSION(NCV), INTENT(IN) :: CV        ! coefficients for Vhr fourier series 
      INTEGER, INTENT(IN) :: Nsiv                      ! symmetry number for torsional potential function
      INTEGER, INTENT(IN) :: Nsif                      ! symmetry number for torsional rotational constant (Bhr)
      INTEGER, INTENT(IN) :: IMAX                      ! grid points and maximum no. of elements in EV array
      REAL(8), INTENT(OUT) :: zpe                      ! zero point energy
      CHARACTER(LEN=5), INTENT(IN) :: Vhr              ! model for torsional potential energy function (Vhr)
      CHARACTER(LEN=5), INTENT(IN) :: Bhr              ! model for rotational constant
      REAL(8), INTENT(IN) :: Phav                      ! phase of torsional angle in Vhr
      REAL(8), INTENT(IN) :: Phab                      ! phase of torsional angle in Bhr
      INTEGER, INTENT(IN) :: Nsig                      ! symmetry number

      INTEGER :: NN                               ! number of elements in vector EV where EV(NN) <= Emax
      REAL(8) :: Emax                             ! maximum energy in Master Equation (corresponding to JMAX and DELE)
      REAL(8), DIMENSION(IMAX) :: EV              ! eigenvalues with Imax elements
      REAL(8), DIMENSION(JMAX) :: AT              ! Beyer-Swinehart array for density of states
	 INTEGER, DIMENSION(IMAX) :: IR              ! array of bin numbers for energy levels
	 INTEGER :: I , J, L, LL                     ! index numbers
	 REAL(8) :: tp                               ! temporary variable
      REAL(8) :: B                         ! average rotational constant, defined as: B = (Integral of B(x)dx from 0 to 2*PI) / (2*PI)
      REAL(8) :: V                         ! hindrance barrier height (from minimum)
      REAL(8) :: PI

        PI=acos(-1.0d0)
        Emax=(JMAX-1)*DELE + 1000.              ! *****raise max energy to account for later subtraction of zpe******

        IF((Bhr.EQ.'Bhrd2').OR.(Bhr.EQ.'BHRD2')) THEN   ! Average rotational constant
                B=0.0d0
                DO I=1, Nsif
                  B=B + (CB(I)/I)*((PI*2)**(I-1))
                ENDDO
        ELSEIF((Bhr.EQ.'Ihrd2').OR.(Bhr.EQ.'IHRD2')) THEN       ! Average moment of inertia
                B=0.0d0
                DO I=1, Nsif
                   B=B + (CB(I)/I)*((PI*2)**(I-1))
                ENDDO
                B=16.85763d0 / B        ! convert from I to B
        ELSEIF((Bhr.EQ.'Ihrd1').OR.(Bhr.EQ.'IHRD1')) THEN       ! Average moment of inertia
                B=16.85763d0/CB(1)      ! convert from I to B
        ELSEIF((Bhr.EQ.'Bhrd1').OR.(Bhr.EQ.'BHRD1')) THEN       ! Average rotational constant
                B=CB(1)
        ELSE
                write(*,*) "ERROR at INPUT for Bhr or Ihr"
        ENDIF

        CALL GHRLEV(EV,NN,IMAX,Emax,B,NCV,NCB,CV,CB,Nsiv,Nsif,
     &          Vhr,Bhr,Phav,Phab)

	zpe=EV(1)

c*********************************************************************
c     The following lines are for output of hindered rotor eigenvalues
c*********************************************************************
        write(12,*) 'Eigenvalues for Separable 1-D Gen. Hindered Rotor'
        write(12,*) '          [uhrlev]'
        write(12,*) '    i     E(i)            E(i)-zpe'

        DO I=1, IMAX
          write(12, 99) I, EV(I), (EV(I)-zpe)
        ENDDO
99    FORMAT (3x, I5, 1x, ES12.5 , 5x, ES12.5 )
c*********************************************************************

	DO I=1, IMAX
		tp=EV(I)-zpe
		IR(I)= NINT( tp/DELE )
C......		write(*,*) I, tp
	ENDDO 

      DO J=1, JMAX                 ! initialize AT(i)
		AT(J)=0.0d0
      ENDDO 

	DO J=1, IMAX
		LL=IR(J)
		DO L=1, JMAX - LL
			AT(L+LL) = AT(L+LL) + T(L)
		ENDDO
	ENDDO
	DO J=1, JMAX
		T(J)=AT(J)/Nsig
		AT(J)=0.0d0
	ENDDO 

	RETURN
	END SUBROUTINE UHRLEV

!-----------------------------------------------------------------------
	 SUBROUTINE GHRLEV(EV,NN,IMAX,Emax,B,NCV,NCF,CV,CF,
     &		Nsiv,Nsif,Vhr,Bhr,Phav,Phab)
!
! GHRLEV: a code to calculate eigenvalues for a general 1D-hindered internal rotation 
! Copyright (C) 2009 Lam T. Nguyen
!
!	DATE:	Feb. 13, 2009
!
! Lam T. Nguyen          or contact:    John R. Barker
! NGUYENLT@umich.edu                    jrbarker@umich.edu
! University of Michigan                University of Michigan
! Ann Arbor, MI 48109-2143              Ann Arbor, MI 48109-2143
!                                       (734) 763 6239
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: NN                  ! number of elements in vector EV where EV(NN) <= Emax
      REAL(8), INTENT(IN) :: Emax                 ! maximum energy in Master Equation
      REAL(8), DIMENSION(IMAX), INTENT(OUT) :: EV ! eigenvalues with Imax elements
      INTEGER, INTENT(IN) :: IMAX                 ! grid points and maximum no. of elements in EV array
      REAL(8), INTENT(IN) :: B                    ! average rotational constant, defined as: B = (Integral of B(x)dx from 0 to 2*PI) / (2*PI)
      CHARACTER(LEN=5), INTENT(IN) :: Vhr         ! model for torsional potential energy function (Vhr)
      INTEGER, INTENT(IN) :: NCV                  ! number of elements in CV
      INTEGER, INTENT(IN) :: Nsiv                 ! symmetry number for torsional potential function
      REAL(8), INTENT(IN) :: Phav                 ! phase of torsional angle in Vhr
      REAL(8), DIMENSION(NCV), INTENT(IN) :: CV   ! coefficients for Vhr fourier series 
      CHARACTER(LEN=5), INTENT(IN) :: Bhr         ! model for rotational constant
      INTEGER, INTENT(IN) :: NCF                  ! number of elements in CF
      INTEGER, INTENT(IN) :: Nsif                 ! symmetry number for torsional rotational constant (Bhr)
      REAL(8), INTENT(IN) :: Phab                 ! phase of torsional angle in Bhr
      REAL(8), DIMENSION(NCF), INTENT(IN) :: CF   ! coefficients for rotational constant fourier series

      REAL(8) :: Vo                               ! Vmax of V(hr)
      REAL(8) :: Vm                               ! Vmin of V(hr)
      REAL(8) :: zpe                              ! zero point energy
!
!.....CALVo used to find Vmax (Vo) and Vmin (Vm) of torsional potential energy function
	 CALL CALVo(Vo,Vm,CV,NCV,Nsiv,Vhr,Phav)

	 CALL ODQHR(EV,IMAX,NCV,NCF,CV,CF,Nsiv,Nsif,Vhr,Bhr,Phav,Phab,Vm)
	 CALL CALEI(EV,NN,IMAX,Emax,B)

        zpe = EV(1)
        write(*,9991) Vo, Vm, EV(1)
9991    FORMAT(11X,'Torsional potential: ',1X,'Vmax (cm-1) = ',F8.1,1X,
     &  '; Vmin (cm-1) = ',F8.1,1X,'; zpe (cm-1) = ',F5.1)
     
!*********************************************************************
!     The following lines are for output of hindered rotor eigenvalues
!*********************************************************************
        write(12,*) '   '
        write(12,9991) Vo, Vm, zpe

!        DO I=1, IMAX+1
!          write(12, 99) I, EV(I), (EV(I)-zpe)
!        ENDDO
!99      FORMAT (3x, I5, 1x, F10.3 , 1x, F10.3 )
!*********************************************************************
	 RETURN 
      END SUBROUTINE GHRLEV

!-----------------------------------------------------------------------
	 SUBROUTINE ODQHR(ER,NX,NCV,NCF,CV,CF,Nsiv,Nsif,
     &	Vhr,Bhr,Phav,Phab,Vm)
!	USE:	To obtain a vector of quantum eigenvalues of 1D-hindered internal rotor, in which  
!		both torsional potential energy function and moment inertia are explicitly treated as 
!		functions of internal rotational angle, i.e. non-rigid rotor, based on Meyer's algorithm:
!		R. Meyer J. Chem. Phys. 52 (1970) 2053-2059.  
!
!	AUTHORS:	Lam T. Nguyen and John R. Barker 
!
!	DATE:	Jan. 26, 2009
!     
!	 NX  : Number of grid points
!        XA  : Position on a grid.
!        AR  : Hamiltonian Matrix, i.e. Meyer's symmetrical matrix
!        ER  : Vector of eigenvalues 
!        ZR  : Eigenvectors (X,Y); where
!                                      X : Wavefunction.
!                                      Y : Energy level.       
!        NPR      : 0 for eigenvalue only and 1 for both eigenvalue and eigenvector 
!        RMIN     : Starting point of grid
!        RMAX     : End point of grid
!        ZL       : Grid length
!        DX       : Grid spacings
!
!	INPUT:	1) NCV, NCF	= Number of coefficients derived from fitting rotational energy profile and moment inertia
!		2) CV(NCV)	= Vector of coefficients from rotational energy profile, with NCV elements
!		3) CF(NCF)	= Vector of coefficients from moment inertia (rotational constant), with NCF elements
!		4) Nsig		= rotational symmetry number
!	OUTPUT:	1) ER 		= Vector of eigenvalues of 1D-hindered internal rotor, with NX=501 elements 
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX                   ! grid points and maximum no. of elements in EV array
      REAL(8), DIMENSION(NX), INTENT(OUT) :: ER   ! eigenvalues with Imax elements
      INTEGER, INTENT(IN) :: NCV                  ! number of elements in CV
      INTEGER, INTENT(IN) :: NCF                  ! number of elements in CF
      REAL(8), DIMENSION(NCV), INTENT(IN) :: CV   ! coefficients for Vhr fourier series 
      REAL(8), DIMENSION(NCF), INTENT(IN) :: CF   ! coefficients for rotational constant fourier series
      INTEGER, INTENT(IN) :: Nsiv                 ! symmetry number for torsional potential function
      INTEGER, INTENT(IN) :: Nsif                 ! symmetry number for torsional rotational constant (Bhr)
      CHARACTER(LEN=5), INTENT(IN) :: Vhr         ! model for torsional potential energy function (Vhr)
      CHARACTER(LEN=5), INTENT(IN) :: Bhr         ! model for rotational constant
      REAL(8), INTENT(IN) :: Phav                 ! phase of torsional angle in Vhr
      REAL(8), INTENT(IN) :: Phab                 ! phase of torsional angle in Bhr
      REAL(8), INTENT(IN) :: Vm
      INTEGER, PARAMETER :: NPR = 0
      REAL(8), DIMENSION(NX,NX) :: ZR
      REAL(8), DIMENSION(NX,NX) :: AR
      REAL(8), DIMENSION(NX) :: XA 
      REAL(8), DIMENSION(NX,NX) :: D
      REAL(8), DIMENSION(NX) :: F
      REAL(8), DIMENSION(NX) :: V
      REAL(8), DIMENSION(NX,NX) :: DT
      REAL(8), DIMENSION(NX,NX) :: T
      REAL(8), DIMENSION(NX,NX) :: FT
	 REAL(8), DIMENSION(NX) :: FV1
	 REAL(8), DIMENSION(NX) :: FV2
	 INTEGER :: IERR
!.....Now compute Hamiltonian matrix:

	 CALL CMM(NX,AR,NCV,NCF,CV,CF,Nsiv,Nsif,Vhr,Bhr,Phav,Phab,Vm)
					
!.....Now call eigenvalue solver using RS subroutine of EISPACK.

	 CALL rs(NX,NX,AR,ER,NPR,ZR,FV1,FV2,IERR) 

!...... Alterative option: using subroutine DSYEV from LAPACK 
!......	CALL DSYEV( 'V', 'U', NX, AR, NX, ER, WORK, NX*3, IERR )

	 RETURN 
      END SUBROUTINE ODQHR

!-----------------------------------------------------------------------
	SUBROUTINE CMM(NX,AR,NCV,NCF,CV,CF,Nsiv,Nsif,Vhr,Bhr,Phav,Phab,Vm)
!
! CMM stands for Constructing Meyer's Matrix, AR, for 1D-hindered internal rotor using Meyer's method, 
! R. MEYER J. Chem. Phys. 52 (1970) 2053.
! Both rotational energy profile and moment inertia are properly treated as functions of 
! internal rotation angle, i.e. non-rigid rotor
! INPUT:
! 1) NCF, NCV = Number of coefficients drived from fitting rotational energy profile and moment inertia
! 2) CV(NCV) = Vector of coefficients from rotational energy profile, with NCV elements
! 3) CF(NCF) = Vector of coefficients from moment inertia, with NCF elements
! 4) Nsig = Rotational symmetry number, it can be either one or two or three or whatever 
! OUTPUT: 1) AR(NX,NX) = Meyer's symmetrical matrix, with a size of NX--number of grid points   

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX                   ! grid points and maximum no. of elements in EV array
      REAL(8), DIMENSION(NX,NX), INTENT(OUT) :: AR
      INTEGER, INTENT(IN) :: NCV                  ! number of elements in CV
      INTEGER, INTENT(IN) :: NCF                  ! number of elements in CF
      REAL(8), DIMENSION(NCV), INTENT(IN) :: CV   ! coefficients for Vhr fourier series 
      REAL(8), DIMENSION(NCF), INTENT(IN) :: CF   ! coefficients for rotational constant fourier series
      INTEGER, INTENT(IN) :: Nsiv                 ! symmetry number for torsional potential function
      INTEGER, INTENT(IN) :: Nsif                 ! symmetry number for torsional rotational constant (Bhr)
      CHARACTER(LEN=5), INTENT(IN) :: Vhr         ! model for torsional potential energy function (Vhr)
      CHARACTER(LEN=5), INTENT(IN) :: Bhr         ! model for rotational constant
      REAL(8), INTENT(IN) :: Phav                 ! phase of torsional angle in Vhr
      REAL(8), INTENT(IN) :: Phab                 ! phase of torsional angle in Bhr
      REAL(8), INTENT(IN) :: Vm

      REAL(8), DIMENSION(NX,NX) :: ZR
      REAL(8), DIMENSION(NX) :: XA 
      REAL(8), DIMENSION(NX,NX) :: D
      REAL(8), DIMENSION(NX,NX) :: DT
      REAL(8), DIMENSION(NX) :: F
      REAL(8), DIMENSION(NX) :: V
      REAL(8), DIMENSION(NX,NX) :: T
      REAL(8), DIMENSION(NX,NX) :: FT
	 REAL(8), DIMENSION(NX) :: FV1
	 REAL(8), DIMENSION(NX) :: FV2
	 REAL(8) :: PI
	 REAL(8) :: RMIN
	 REAL(8) :: RMAX
	 REAL(8) :: ZL
	 REAL(8) :: DX
	 REAL(8) :: TP1, TP2, SU 
	 INTEGER :: I, J, K, ImJ

	 PI=acos(-1.D0)
!
!.....Set up grid
	RMIN=0.0D0
	RMAX=PI*2
	ZL=(RMAX-RMIN)
	DX=ZL/NX

!.....Now compute Meyer's symmetrical matrix:

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
      END SUBROUTINE CMM

!-----------------------------------------------------------------------
	 REAL(8) FUNCTION DELTA(I,J)
!     Kronecker delta function
!	 
	   INTEGER, INTENT(IN) :: I, J
	   IF(I.EQ.J) THEN
		DELTA=1.0D0
	   ELSE
		DELTA=0.0D0
	   ENDIF
	 END FUNCTION DELTA

!-----------------------------------------------------------------------
	SUBROUTINE VX(V,X,CV,N,Nsig,Vhr,Vm)
!    Computes potential energy
!
	 IMPLICIT NONE
	 REAL(8), INTENT(OUT) :: V                   ! potential
	 REAL(8), INTENT(IN) :: X                    ! dihedral
      INTEGER, INTENT(IN) :: N                    ! number of elements in CV
      REAL(8), DIMENSION(N), INTENT(IN) :: CV     ! coefficients for Vhr fourier series 
      INTEGER, INTENT(IN) :: Nsig                 ! symmetry number for potential
      CHARACTER(LEN=5), INTENT(IN) :: Vhr         ! model for torsional potential energy function (Vhr)
      REAL(8) :: Vm
      INTEGER :: I, Nc

	IF( Vhr .EQ. 'VHRD1' ) THEN
		V=0.0d0
		DO I=1, N
			V=V+CV(I)*(1.0d0-dcos(X*Nsig*I))/2.0d0
		ENDDO
	ELSEIF( Vhr .EQ. 'VHRD2' ) THEN
		V=CV(1)
		DO I=2, N
			V=V+CV(I)*dcos(X*Nsig*(I-1))
		ENDDO
        ELSEIF( Vhr .EQ. 'VHRD3' ) THEN
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
		write(*,*) "ERROR at INPUT for torsion function: Vhrd= ",Vhr
		STOP
	ENDIF

	V=V-Vm

	RETURN
	END SUBROUTINE VX

!-----------------------------------------------------------------------
	 SUBROUTINE FX(F,X,CF,N,Nsig,Bhr)
	 IMPLICIT NONE
      REAL(8), INTENT(OUT) :: F                   ! rotational constant
      REAL(8), INTENT(IN) :: X                    ! dihedral
      INTEGER, INTENT(IN) :: N                    ! number of elements in CF
      REAL(8), DIMENSION(N), INTENT(IN) :: CF     ! coefficients for rotational constant fourier series
      INTEGER, INTENT(IN) :: Nsig                 ! symmetry number for rotational constant
      CHARACTER(LEN=5), INTENT(IN) :: Bhr         ! model for rotational constant
	 REAL(8), PARAMETER:: CON = 16.85763D0
      INTEGER :: I, Nc

	 IF( Bhr .EQ. 'BHRD1' ) THEN
		F=CF(1)
		DO I=2, N
			F=F+CF(I)*dcos(X*Nsig*(I-1))
		ENDDO
        ELSEIF( Bhr .EQ. 'IHRD1' ) THEN
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
	 END SUBROUTINE FX

!-----------------------------------------------------------------------
	SUBROUTINE CALEI(EV,NN,Imax,Emax,B)
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! CALEI: a code to compute/select eigenvalues of torsional motion based on Emax 
! Copyright (C) 2009 Lam T. Nguyen and John R. Barker
!
! John R. Barker
! jrbarker@umich.edu
! University of Michigan
! Ann Arbor, MI 48109-2143
! (734) 763 6239
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!	INPUT:	1) EV(Nmax)	: a vector of eigenvalues with Nmax elements
!		2) Nmax		: number of elements of vector EV
!		3) Emax		: maximum energy in Master Equation
!		4) B		: average rotational constant, defined as:
!			B = (Integral of B(x)dx from 0 to 2*PI) / (2*PI)
!	OUTPUT:	1) EV(Nmax)
!		2) Nmax	

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Imax                 ! maximum no. of elements in EV array
      REAL(8), DIMENSION(Imax), INTENT(OUT) :: EV ! eigenvalues with Nmax elements
      INTEGER, INTENT(OUT) :: NN                  ! number of elements in vector EV where EV(NN) <= Emax
      REAL(8), INTENT(IN) :: Emax                 ! maximum energy in Master Equation
      REAL(8), INTENT(IN) :: B                    ! average rotational constant, defined as: B = (Integral of B(x)dx from 0 to 2*PI) / (2*PI)
      INTEGER :: I
      REAL :: E

! IF The largest eigenvalue is larger than the maximum energy, 
! SEARCH for an eigenvalue, which is <= Emax

	 I = Imax
	 E = EV(I)
	 DO WHILE (E.GT.Emax)
		I=I-1
		E = EV(I)
	 END DO
	 NN = I

	 RETURN
      END SUBROUTINE CALEI

!-----------------------------------------------------------------------
        SUBROUTINE CALVo(Vo,Vm,CV,NCV,Nsiv,Vhr,Pha)
! CALVo is used to compute Vmax and Vmin 
! of the torsional potential energy function 
!

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

!-----------------------------------------------------------------------

      END MODULE HRMOD


