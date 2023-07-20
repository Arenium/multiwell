      MODULE vibcalcs
!
! Copyright (C) 2009, 2017 John R. Barker, Thanh Lam Nguyen, and 
!                          Collin G. L. Li
!
! Authors: Thanh Lam Nguyen, Collin G. L. Li, and John R. Barker
!          nguyenlt@umich.edu
!          September 2009, August 2017
!
! John R. Barker
! jrbarker@umich.edu
! University of Michigan
! Ann Arbor, MI 48109-2143
! (734) 763 6239
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License (version 2)
! as published by the Free Software Foundation.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! See the GNU General Public License:
!    Free Software Foundation, Inc.
!    59 Temple Place - Suite 330
!    Boston, MA 02111-1307, USA.
!      
      IMPLICIT NONE
      SAVE
      
      INTEGER, PARAMETER :: nsm = 100             ! maximum no. of vibrational modes
      INTEGER :: ns                               ! number of vibrational modes
      REAL(8), DIMENSION(nsm) :: wa               ! harmonic frequencies (cm-1) with zero of energy at the potential minimum
      REAL(8), DIMENSION((nsm)) :: w0             ! frequencies (cm-1) with zero of energy at zpe
      REAL(8), DIMENSION(nsm) :: wf               ! fundamental frequencies (cm-1; for 0-1 transitions with all vj = 0)
      REAL(8), DIMENSION(nsm,nsm) :: xa           ! X(i,j) anharmonicities
      REAL(8), DIMENSION(nsm,nsm,nsm) :: ya       ! Y(i,j,k) anharmonicities
      REAL(8), DIMENSION(nsm,nsm,nsm,nsm) :: za   ! Z(i,j) anharmonicities
      REAL(8) :: zpe                              ! vibrational zero point energy (cm-1) 
      REAL(8) :: ave                              ! average harmonic <wa> (cm-1) 
      REAL(8) :: av0                              ! average frequency <w0> (cm-1) 
      REAL(8) :: avf                              ! average fundamental <wf> (cm-1) 
      REAL(8) :: Viblo                            ! smallest harmonic frequency (cm-1) 
      REAL(8) :: Sfac                             ! Factorial ns!
      REAL(8) :: Prod                             ! product of harmonic frequencies: PROD (wa(i))
      REAL(8) :: Meansq                           ! mean square harmonic frequency: SUM( wa(i)**2 )/ns
      REAL(8) :: Sqmean                           ! square of mean harmonic freq:  [ SUM( wa(i) )/ns ]**2
      INTEGER, DIMENSION(nsm) :: nv , nvold       ! vibrational quantum numbers (current and old)
      
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
      CONTAINS

!----------------------------------------------------------------------
      SUBROUTINE READ_WXYZ( KIN , NY, NZ, WW ) 
!
!      Read vibrational frequencies and anharmonicities
!
!       ns  = number of vibrational frequencies
!       wa(i): harmonic frequencies with zero of energy at the potential minimum
!          Evib = SUMi( wai*(vi+1/2) ) + SUMi( SUMj( Xij*(vi+1/2)*(vj+1/2) ) )
!       w0(i): frequencies with zero of energy at zpe
!          Evib = SUMi( w0i*vi ) + SUMi( SUMj( Xij*vi*vj ) )
!       wf(i) = fundamental frequencies (for 0-1 transitions with all vj = 0)
!      
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: KIN               ! Unit no. for reading input
      INTEGER(4), INTENT(IN) :: NY                ! no. of elements in Y anharmonicity matrix to be read in
      INTEGER(4), INTENT(IN) :: NZ                ! no. of elements in Z anharmonicity matrix to be read in
      CHARACTER(len=2), INTENT(IN) :: WW          ! keyword ('we', 'w0', or 'wf') for the three types of frequencies
      INTEGER(4) :: i , j , k , l  , I1 , J1 , K1 , L1
      CHARACTER(LEN=5) :: KEYWORD
!
      IF ( WW .EQ. 'WE' ) THEN
        READ(KIN,*)     (wa(i), i=1,ns)
      ELSEIF ( WW .EQ. 'W0') THEN
        READ(KIN,*)     (w0(i), i=1,ns)
      ELSEIF ( WW.EQ.'WF' ) THEN
        READ(KIN,*)     (wf(i), i=1,ns)
      ENDIF
!                                                        X matrix (vibrational anharmonicity)
      READ(KIN,9011) KEYWORD
      CALL ucase ( KEYWORD )
      IF ( KEYWORD .EQ. 'UPPER')  THEN
        DO j = 1, ns
          READ(KIN,*) (xa(j,i), i=j, ns)
          IF ( j .LT. ns) THEN
            DO k = j+1 , ns
              xa(k,j) = xa(j,k)
            END DO
          ENDIF
        ENDDO
      ELSEIF ( KEYWORD .EQ. 'LOWER') THEN
        DO j = 1, ns
          READ(KIN,*) (xa(j,i), i=1,j)
          IF ( j .GT. 1) THEN
            DO k = 1 , j-1
              xa(k,j) = xa(j,k)
            END DO
          ENDIF
        ENDDO
      ENDIF
!                                                        Y matrix (vibrational anharmonicity)
       DO i=1, ns
         DO j=1, ns
             DO k=1, ns
                ya(i,j,k)=0.0d0
             ENDDO
         ENDDO
       ENDDO

        DO i=1, NY
         READ(KIN,*) I1, J1, K1, ya(I1,J1,K1)
         ya(K1,J1,I1)=ya(I1,J1,K1)
         ya(K1,I1,J1)=ya(I1,J1,K1)
         ya(J1,I1,K1)=ya(I1,J1,K1)
         ya(J1,K1,I1)=ya(I1,J1,K1)
         ya(I1,K1,J1)=ya(I1,J1,K1)
        ENDDO 
!                                                        Z matrix (vibrational anharmonicity)
        DO i=1, ns
          DO j=1, ns
            DO k=1, ns    
              DO l=1, ns
                za(i,j,k,l)=0.0d0
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        DO i=1, NZ
         READ(KIN,*) I1, J1, K1, L1, za(I1,J1,K1,L1)
         za(I1,J1,L1,K1)=za(I1,J1,K1,L1)
         za(I1,K1,J1,L1)=za(I1,J1,K1,L1)
         za(I1,K1,L1,J1)=za(I1,J1,K1,L1)
         za(I1,L1,K1,J1)=za(I1,J1,K1,L1)
         za(I1,L1,J1,K1)=za(I1,J1,K1,L1)
         za(J1,I1,K1,L1)=za(I1,J1,K1,L1)
         za(J1,I1,L1,K1)=za(I1,J1,K1,L1)
         za(J1,K1,I1,L1)=za(I1,J1,K1,L1)
         za(J1,K1,L1,I1)=za(I1,J1,K1,L1)
         za(J1,L1,I1,K1)=za(I1,J1,K1,L1)
         za(J1,L1,K1,I1)=za(I1,J1,K1,L1)
         za(K1,I1,J1,L1)=za(I1,J1,K1,L1)
         za(K1,I1,L1,J1)=za(I1,J1,K1,L1)
         za(K1,J1,I1,L1)=za(I1,J1,K1,L1)
         za(K1,J1,L1,I1)=za(I1,J1,K1,L1)
         za(K1,L1,I1,J1)=za(I1,J1,K1,L1)
         za(K1,L1,J1,I1)=za(I1,J1,K1,L1)
         za(L1,I1,J1,K1)=za(I1,J1,K1,L1)
         za(L1,I1,K1,J1)=za(I1,J1,K1,L1)
         za(L1,J1,I1,K1)=za(I1,J1,K1,L1)
         za(L1,J1,K1,I1)=za(I1,J1,K1,L1)
         za(L1,K1,I1,J1)=za(I1,J1,K1,L1)
         za(L1,K1,J1,I1)=za(I1,J1,K1,L1)
        ENDDO

9011  FORMAT(A5)

      RETURN
      END SUBROUTINE READ_WXYZ
      
!-----------------------------------------------------------------------------------------
      SUBROUTINE convib ( WW ) 
!
!     From input as wa, w0, or wf frequencies, convert to the other forms:
!       wa(i): harmonic frequencies with zero of energy at the potential minimum
!          Evib = SUMi( wai*(vi+1/2) ) + SUMi( SUMj( Xij*(vi+1/2)*(vj+1/2) ) )
!       w0(i): frequencies with zero of energy at zpe
!          Evib = SUMi( w0i*vi ) + SUMi( SUMj( Xij*vi*vj ) )
!       wf(i) = fundamental frequencies (for 0-1 transitions with all vj = 0)
!
!     zpe = zero point energy, based on We and X matrix
!
!     Average Frequences:
!       ave = average wa(i)
!       av0 = average w0(i)
!       avf = average wf(i)
!
!     [Herzberg, Infrared and Raman Spectra (D. van Nostrand Co., 1945), p. 206ff]
!     Eq. numbers from Herzberg
      
      IMPLICIT NONE
      INTEGER i , j
      CHARACTER(len=2) :: WW            ! keyword ('we', 'w0', or 'wf') for the three types of frequencies
      
      IF ( WW .EQ. 'WE' ) THEN

      DO i = 1 , ns                               ! start conversion
        w0(i) = wa(i) + xa(i,i)                   ! Eq. (II,273)
        wf(i) = wa(i) + 2.0d+00*xa(i,i)           ! fundamental
          DO j = 1, ns
            IF ( i .NE. j ) THEN
            w0(i) = w0(i) + 0.5d+00*xa(i,j)       ! Eq. (II,273)
            wf(i) = wf(i) + 0.5d+00*xa(i,j)       ! fundamental
            ENDIF
          END DO
      END DO                                      ! end conversion

      ELSEIF ( WW .EQ. 'W0' ) THEN

      DO i = 1 , ns                               ! start conversion
        wa(i) = w0(i) - xa(i,i)                   ! Eq. (II,273)
          DO j = 1, ns
            IF ( i .NE. j ) THEN
            wa(i) = wa(i) - 0.5d+00*xa(i,j)       ! Eq. (II,273)
            ENDIF
          END DO
      END DO                                      ! end conversion
      DO i = 1 , ns                               ! start conversion
        wf(i) = wa(i) + 2.0d+00*xa(i,i)           ! fundamental
          DO j = 1, ns
            IF ( i .NE. j ) THEN
            wf(i) = wf(i) + 0.5d+00*xa(i,j)       ! fundamental
            ENDIF
          END DO
      END DO                                      ! end conversion

      ELSEIF ( WW .EQ. 'WF' ) THEN

      DO i = 1 , ns                               ! start conversion
        wa(i) = wf(i) - 2.d+00*xa(i,i)
          DO j = 1, ns
            IF ( i .NE. j ) THEN
            wa(i) = wa(i) - 0.5d+00*xa(i,j)       ! Eq. (II,273)
            ENDIF
          END DO
      END DO                                      ! end conversion
      DO i = 1 , ns                               ! start conversion
        w0(i) = wa(i) + xa(i,i)                   ! Eq. (II,273)
          DO j = 1, ns
            IF ( i .NE. j ) THEN
            w0(i) = w0(i) + 0.5d+00*xa(i,j)       ! Eq. (II,273)
            ENDIF
          END DO
      END DO                                      ! end conversion

      ELSE
        WRITE(*,*) 'FATAL: Vib data type not recognized: ', WW
        
      ENDIF
      
      zpe = 0.0d+00
      ave = 0.0d+00
      av0 = 0.0d+00
      avf = 0.0d+00
      Viblo = 10000.
      Sfac = 1.d+00
      Prod = 1.d+00
      Meansq = 0.d+00
      DO i = 1 , ns                          ! start ZPE
        zpe = zpe + 0.5d+00*wa(i)            ! Eq. (II,267)
        ave = ave + wa(i)                    ! summed frequency
        av0 = av0 + w0(i)                    ! summed frequency
        avf = avf + wf(i)                    ! summed frequency
        Sfac = Sfac*I                        ! factorial ns!
        Prod = Prod*wa(i)                    ! product of harmonic frequencies
        Meansq = Meansq + ( wa(i)**2 )/ns    ! mean squared harmonic freq
        Viblo = MIN( Viblo , wa(i) )         ! find lowest frequency
        DO j = 1, i
          zpe = zpe + 0.25d+00*xa(j,i)       ! Eq. (II,267)
        END DO
      END DO                                 ! end 

      ave = ave/ns                           ! average frequency
      av0 = av0/ns                           ! average frequency
      avf = avf/ns                           ! average frequency
      Sqmean = (2.d+0*zpe/ns)**2             ! squared mean frequency

      RETURN
      END SUBROUTINE convib
      
!-----------------------------------------------------------------------------------------
      REAL(8) FUNCTION feswitch( Emax,DELE,LTMODE,MVAL,KEYWORD)
! Computes energy for switching from direct count to Wang-Landau algorithm
!
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: Emax                 ! maximum energy of calculation
      REAL(8), INTENT(IN) ::  DELE                ! energy grain size
      REAL(8), INTENT(IN) :: MVAL                 ! MVAL = switch point energy if manually set
      CHARACTER(6), INTENT(IN):: LTMODE           ! key word for how switch energy is determined: "PRAC" (practical), "AUTO" (automated), or 'MAN" (manual)
      CHARACTER(6), INTENT(IN):: KEYWORD          ! designates no. of Wang-Landau stochastic trials; assumed to be upper case 
      REAL(8) :: Time, Gx, C, Beta 
      Real(8) :: a1, alpha, b1, Be                ! Tadensum = b1*(Emax - x)^Be, Tnewloops = a1*Gx^alpha
      Real(8) :: Eswitch                          ! switch point energy
      REAL(8) :: dgxsolver, gxsolver
      EXTERNAL gxsolver, dgxsolver
!
      Time = MVAL
      Gx = sqrt((Time-0.830d0)*(1/9.350d0)*10.0d0**10)
      Beta = ((ns-1.0d0)/ns)*(Meansq/Sqmean)
!
!       G(X) SOLVER  --Practical limit
!
!     C = (s!w(s)*G)^(1/s) = (Sfac*Prod*Gx)^(1/ns) 

      IF ( LTMODE.EQ.'PRAC' ) THEN
        C = (Sfac*Prod)
!        C = DExp((1.0d0/ns)*DLOG(Sfac*Prod*Gx))      !Recursive Method
        feswitch = gxsolver( C, zpe, DELE, Beta, ns, Gx, Emax )    
        IF(feswitch .GT. Emax)THEN
	      feswitch = Emax
        ENDIF
        feswitch = INT( feswitch/DELE )*DELE
        write(*,99052) LTMODE, feswitch
!      ENDIF
!
!       dG(x) SOLVER  --Optimum limit
!
!     C = a1*alpha*ns/(b1*Be*(Sfac*Prod)^alpha)
!
      ELSEIF ( LTMODE.EQ.'AUTO' ) THEN
        a1 = 1.0d-5                                         !empirical results  1.616E-5 is the calculated result
        b1 = 0.07544d-00*ns-0.4317d-00                              !empirical results
        alpha = 0.9813d-00
        Be = 1.1d-00										
        C = a1*alpha*ns/(b1*Be*(Sfac*Prod)**alpha)

        feswitch = dgxsolver( Emax,C,zpe,Be,ns,alpha,b1,Beta,KEYWORD )
        
        feswitch = INT( feswitch/DELE )*DELE                ! makes eswitch an integer multiple of the grain size
        write(*,99052) LTMODE, feswitch

      ELSEIF (LTMODE.EQ.'MAN'.OR.LTMODE.EQ.'Man'.OR.LTMODE.EQ.'man')THEN
        feswitch = INT( MVAL/DELE )*DELE                    ! makes eswitch an integer multiple of the grain size
        write(*,99052) LTMODE, feswitch

      ELSE
        write(*,*) 'FATAL: unknown Eswitch KEYWORD: ', LTMODE
        write(3,*) 'FATAL: unknown Eswitch KEYWORD: ', LTMODE
        write(4,*) 'FATAL: unknown Eswitch KEYWORD: ', LTMODE
        write(5,*) 'FATAL: unknown Eswitch KEYWORD: ', LTMODE
        STOP
      ENDIF

      RETURN 
99052 FORMAT ( A6, ': Eswitch = ', F7.1 )
      END FUNCTION feswitch
      
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------      
      
      END MODULE vibcalcs


