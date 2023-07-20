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

      MODULE utili_mod
      
!     The following procedures are contained in this module:
!
!     SUBROUTINE NTERP       Copyright 2001 - 2021 John R. Barker
!     SUBROUTINE NTERP2      Copyright 2001 - 2021 John R. Barker, Philip Stimac
!     FUNCTION fxnofe        Copyright 2001 - 2021 John R. Barker
!     FUNCTION fxnexp        Copyright 2001 - 2021 John R. Barker
!     FUNCTION rxnofe        Copyright 2001 - 2021 John R. Barker
!     FUNCTION XLINT         Copyright 2001 - 2021 John R. Barker
!     SUBROUTINE DateTime    Copyright 2001 - 2021 John R. Barker
!     FUNCTION RAN1          Copyright 1992 American Institute of Physics
!     function lenstr         
!     SUBROUTINE ucase


      CONTAINS
      
!  ---------------------------------------------------------------------
      SUBROUTINE NTERP(E,IE,X)
!
!       For Interpolating double arrays
!
C       Calculate nearest element in Double Array: upper bound= Isize+Imax1
c
c       'Double arrays' have two sections:
c              segment 1 consists of equally spaced (Egrain1) values ranging
c                     from E=0 to Emax1
c              segment 2 consists of equally spaced values from E=0 to Emax2;
c                     the spacing depends on the number of arrays elements
c                     remaining: Isize-Emax1/Egrain1-1
c
c       E        = energy (units: cm-1)
c       Egrain1  = energy grain of lower segments in 'double arrays' (units: cm-1)
c       Egrain2  = energy grain of upper segments in 'double arrays' (units: cm-1)
c       Emax1    = maximum energy of 1st segment of double arrays (units: cm-1)
c       Emax2    = maximum energy of 2nd segment of double arrays (units: cm-1)
c       IE       = nearest element index
c       X        = fraction of element difference from IE (units: dimensionless)
 
      USE declare_mod, ONLY: Imax1, Emax1, Egrain1, Egrain2, Isize
      IMPLICIT NONE

      REAL(8), INTENT (IN)  :: E
      REAL(8), INTENT (OUT) :: X
      INTEGER, INTENT (OUT) :: IE
      
        IF ( E .LE. Emax1 ) THEN
          IE = 1 + NINT(E/Egrain1)       ! Nearest integer
          IF (IE .LT. 2) IE = 2
          if (IE .GT. (Imax1-1) ) IE = Imax1-1
          X = E/Egrain1 + 1 - IE
        else  
          IE = Imax1 + 1 + NINT(E/Egrain2)       ! Nearest integer
          if( IE .GE. Isize ) IE = Isize - 1
          if( IE .LT. Imax1 + 2 ) IE = Imax1 + 2
          X = E/Egrain2 + 1 + Imax1 - IE
        end if

      RETURN
      END SUBROUTINE NTERP
      
      
!  ---------------------------------------------------------------------
      SUBROUTINE NTERP2(E,IE,X)
!
!     FOR ARRAYS USED IN ECKART AND TUNNELING ROUTINES
!
C     Calculate nearest element in Double Array: upper bound= Isize+Imax1
c
c     'Double arrays' have two sections:
c           segment 1 consists of equally spaced (Egrain1) values ranging
c                 from E=0 to Emax1
c           segment 2 consists of equally spaced values from E=0 to Emax2;
c                 the spacing depends on the number of arrays elements
c                 remaining: Isize-Emax1/Egrain1-1
c
c     E     = energy (units: cm-1)
c     Egrain1     = energy grain of lower segments in 'double arrays' (units: cm-1)
c     Egrain2     = energy grain of upper segments in 'double arrays' (units: cm-1)
c     Emax1 = maximum energy of 1st segment of double arrays (units: cm-1)
c     Emax2 = maximum energy of 2nd segment of double arrays (units: cm-1)
c     IE    = nearest element index
c     X     = fraction of element difference from IE (units: dimensionless)
 
      USE declare_mod

      IMPLICIT NONE
      REAL(8), INTENT (IN) :: E
      REAL(8), INTENT (OUT) :: X
      INTEGER, INTENT (OUT) :: IE
      INTEGER :: IUP                ! upper bound
      
      IUP = Isize + Imax1

! code modified to correctly interpolate k(E) array with Isize+imax1
! elements; imax1 of which are below the barrier


        if(E.le.Emax1)then
          IE=1+NINT(E/Egrain1)+Imax1
          if(IE.lt.2) IE=2
          if(IE.gt.(2*Imax1-1)) IE=2*Imax1-1
          X=E/Egrain1-(IE-1.0)+Imax1
        else  
         IE = 2*imax1 + 1 + NINT(E/Egrain2)       ! Nearest integer
          if(IE.ge. IUP) IE = IUP - 1
          if(IE.le.(2*Imax1+1)) IE=2*Imax1+2
         X = E/Egrain2 - (IE-2*imax1-1.0)
      end if
      RETURN
      END SUBROUTINE NTERP2

!  ---------------------------------------------------------------------
      REAL(8) FUNCTION fxnofe(E , Mol , ARRAY)
 
C     for a function of E, calculate nearest element in Double Array of size Isize
c
c     'Double arrays' have two sections:
c           segment 1 consists of equally spaced (Egrn1) values ranging
c                 from EN=0 to Emax1
c           segment 2 consists of equally spaced values from EN=0 to Emax2;
c                 the spacing depends on the number of arrays elements
c                 remaining: Ize-Emax1/Egrn1-1
c
c     E     = energy (units: cm-1)
c     Egrain1     = energy grain of lower segments in 'double arrays' (units: cm-1)
c     Egrain2     = energy grain of upper segments in 'double arrays' (units: cm-1)
c     Emax1 = maximum energy of 1st segment of double arrays (units: cm-1)
c     Emax2 = maximum energy of 2nd segment of double arrays (units: cm-1)
c     i     = nearest element index
c     X     = EN - EN(i)

      USE declare_mod, ONLY: Imax1,Emax1,Egrain1,Egrain2

      IMPLICIT NONE
      REAL(8), INTENT(IN)  :: E                       ! Energy
      INTEGER, INTENT(IN)  :: Mol                     ! index for molecule
      REAL(8), INTENT(IN), DIMENSION(:,:) :: ARRAY
      REAL(8) :: X, Y, y1, y2, y3
      INTEGER :: IE
      INTEGER :: IUP          ! upper energy bound of ARRAY
      
      IUP = UBOUND(ARRAY,2)
 
      IF ( E .LE. Emax1 ) THEN
         IE = 1 + NINT(E/Egrain1)               ! Nearest integer
         Y = ARRAY(Mol,IE)                      ! no interpolation
      ELSE
         IE = Imax1 + 1 + NINT(E/Egrain2)       ! Nearest integer
         IF ( IE .LE. Imax1+1 ) IE = Imax1 + 2
         IF ( IE .GE. IUP ) IE = IUP - 1
         X = E/Egrain2 - (IE-imax1-1.0)         ! displacement from IE
         Y1 = ARRAY(Mol,IE-1)
         Y2 = ARRAY(Mol,IE)
         Y3 = ARRAY(Mol,IE+1)
         IF ( X .LE. 0.0d+00) THEN
              Y = Y2 + X*(Y2-Y1)
           ELSEIF ( X .GT. 0.0d+00) THEN
              Y = Y2 + X*(Y3-Y2)
         ENDIF
      ENDIF

      fxnofe = Y
       
      RETURN
      END FUNCTION fxnofe
      
!  ---------------------------------------------------------------------
      REAL(8) FUNCTION fxnexp(E , Mol , ARRAY)
 
C     for a function of E stored as a log( Double Array ) of size Isize
c
c     'Double arrays' have two sections:
c           segment 1 consists of equally spaced (Egrn1) values ranging
c                 from EN=0 to Emax1
c           segment 2 consists of equally spaced values from EN=0 to Emax2;
c                 the spacing depends on the number of arrays elements
c                 remaining: Ize-Emax1/Egrn1-1
c
c     E     = energy (units: cm-1)
c     Egrain1     = energy grain of lower segments in 'double arrays' (units: cm-1)
c     Egrain2     = energy grain of upper segments in 'double arrays' (units: cm-1)
c     Emax1 = maximum energy of 1st segment of double arrays (units: cm-1)
c     Emax2 = maximum energy of 2nd segment of double arrays (units: cm-1)
c     i     = nearest element index
c     X     = EN - EN(i)

      USE declare_mod, ONLY: Imax1,Emax1,Egrain1,Egrain2

      IMPLICIT NONE
      REAL(8), INTENT(IN)  :: E                       ! Energy
      INTEGER, INTENT(IN)  :: Mol                     ! index for molecule
      REAL(8), INTENT(IN), DIMENSION(:,:) :: ARRAY
      REAL(8) :: X, Y, y1, y2, y3
      INTEGER :: IE
      INTEGER :: IUP          ! upper energy bound of ARRAY
      
      IUP = UBOUND(ARRAY,2)
 
      IF ( E .LE. Emax1 ) THEN
         IE = 1 + NINT(E/Egrain1)               ! Nearest integer
         Y = ARRAY(Mol,IE)                      ! no interpolation
      ELSE
         IE = Imax1 + 1 + NINT(E/Egrain2)       ! Nearest integer
         IF ( IE .LE. Imax1+1 ) IE = Imax1 + 2
         IF ( IE .GE. IUP ) IE = IUP - 1
         X = E/Egrain2 - (IE-imax1-1.0)         ! displacement from IE
         Y1 = ARRAY(Mol,IE-1)
         Y2 = ARRAY(Mol,IE)
         Y3 = ARRAY(Mol,IE+1)
         IF ( X .LE. 0.0d+00) THEN
              Y = Y2 + X*(Y2-Y1)
           ELSEIF ( X .GT. 0.0d+00) THEN
              Y = Y2 + X*(Y3-Y2)
         ENDIF
      ENDIF

      fxnexp = exp( Y )
       
      RETURN
      END FUNCTION fxnexp

!  ---------------------------------------------------------------------
      REAL(8) FUNCTION rxnofe(E,Mol,nc,ARRAY)
 
c     For interpolating k(Mol,nc,E), which is stored as log
c
C     a function of energy E, calculate nearest element in Double Array of size Isize
c
c     'Double arrays' have two sections:
c           segment 1 consists of equally spaced (Egrain1) values ranging
c                 from E=0 to Emax1
c           segment 2 consists of equally spaced values from E=0 to Emax2;
c                 the spacing depends on the number of arrays elements
c                 remaining: Isize-Emax1/Egrain1-1
c
c     E     = energy (units: cm-1)
c     Egrain1     = energy grain of lower segments in 'double arrays' (units: cm-1)
c     Egrain2     = energy grain of upper segments in 'double arrays' (units: cm-1)
c     Emax1 = maximum energy of 1st segment of double arrays (units: cm-1)
c     Emax2 = maximum energy of 2nd segment of double arrays (units: cm-1)
c     i     = nearest element index
c     X     = E - E(i)

      USE declare_mod, ONLY: Imax1,Emax1, Egrain1,Egrain2

      IMPLICIT NONE
      REAL(8), INTENT(IN)  :: E                       ! Energy
      INTEGER, INTENT(IN)  :: Mol                     ! index for molecule
      INTEGER, INTENT(IN)  :: nc                      ! index for rxn channel
      REAL(8), INTENT(IN), DIMENSION(:,:,:) :: ARRAY  ! rate constant array (stored as log)
      REAL(8) :: X, Y, Ei , slope
      INTEGER :: i
      INTEGER :: IUP          ! upper energy bound of ARRAY
      
      IUP = UBOUND(ARRAY,3)

      IF ( E .LE. Emax1 ) THEN
         i = 1 + NINT(E/Egrain1)                ! Nearest integer
         Y = ARRAY(Mol,nc,i)                    ! no interpolation
      ELSE
         i = Imax1 + 1 + NINT(E/Egrain2)       ! Nearest integer
         IF ( i.GE.IUP ) i = IUP - 1
         IF ( i.LE.Imax1+1 ) i = Imax1 + 2
         Ei = ( i - 1 - Imax1 )*Egrain2
         X = E - Ei
         IF ( X .GE. 0.0 ) THEN
            slope = ( ARRAY(Mol,nc,i+1) - ARRAY(Mol,nc,i) )/Egrain2
            Y = ARRAY(Mol,nc,i) + X*slope
           ELSE
            slope = ( ARRAY(Mol,nc,i) - ARRAY(Mol,nc,i-1) )/Egrain2
            Y = ARRAY(Mol,nc,i) + X*slope
         ENDIF
      ENDIF

      rxnofe = exp( Y )       ! rate constants stored as log
       
      RETURN
      END FUNCTION rxnofe
      
!  ---------------------------------------------------------------------
      REAL(8) FUNCTION XLINT(Y1,Y2,Y3,X)
C     (previously named SUBROUTINE QUDINT)
c
C     LINEAR INTERPOLATION FROM EQUALLY-SPACED DATA, RETURNING
C     VALUE OF Y AT VALUE OF X:
C
C       Y1 AT X = -1
C       Y2 AT X =  0
C       Y3 AT X = +1
c
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: Y1, Y2, Y3 , X
      REAL(8) :: Y
 
      IF ( X .LE. 0.0d+00) THEN
           Y = Y2 + X*(Y2-Y1)
      ELSEIF ( X .GT. 0.0d+00) THEN
           Y = Y2 + X*(Y3-Y2)
      ENDIF
      
      XLINT = Y
 
      RETURN

      END FUNCTION XLINT
      
!  ---------------------------------------------------------------------
      SUBROUTINE DateTime(KUNIT)
c
c      Date & time subroutine
c
      IMPLICIT NONE
      INTEGER(4) :: KUNIT
      REAL(8) :: duration
      CHARACTER(len=10) date, time, zone
      INTEGER values(8)

      call date_and_time(date,time,zone,values)
      WRITE (KUNIT,99001) values(1), values(2), values(3), values(5), 
     &                    values(6), values(7), values(8)

      CALL CPU_TIME(duration)
      WRITE (KUNIT,99002) duration, duration/3600.
      
      RETURN
      
99001 FORMAT (/'....................................................',/,
     &     10x, I4,'/',I2.2,'/',I2.2,3x,I2.2,':',I2.2,':',
     &          I2.2,'.',I3.3,'   (Local Time)',/
     &     10x, 'year/mm/dd   hr:mi:second'   )
 
99002 FORMAT (17x, 'CPU :   ',f10.3,' s (',f7.3,' hr)',/
     &        '....................................................')

      RETURN
      END SUBROUTINE DateTime
      
!  ---------------------------------------------------------------------

      REAL(8) FUNCTION RAN1(idum)
C     "Minimal" random number generator of Park and Miller with Bays-Durham shuffle
C     and added safeguards.  Returns a uniform random deviate between 0.0 and 1.0
C     (exclusive of the endpoint values).  Call with idum a negative integer to
C     initialize; thereafter do not alter idum between successive deviates in a
C     sequence.  RNMX should approximate the largest floating value that is less
C     than 1.
C
C     W. H. Press and S. A. Teukolsky, Computers in Physics, 6(5), 522-4 (1992).

      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: idum
      
      INTEGER, PARAMETER :: IA = 16807
      INTEGER, PARAMETER :: IM = 2147483647
      INTEGER, PARAMETER :: IQ = 127773
      INTEGER, PARAMETER :: IR = 2836
      INTEGER, PARAMETER :: NTAB = 32
      INTEGER, PARAMETER :: NDIV = 1+INT( (IM-1)/NTAB )
      REAL(8), PARAMETER :: AM = 1.0D+00/IM
      REAL(8), PARAMETER :: EPS = 2.3D-16
      REAL(8), PARAMETER :: RNMX = 1.0D+00-EPS
      REAL(8) :: T
      INTEGER :: iy = 0
      INTEGER, DIMENSION(NTAB) :: iv = 0
      INTEGER :: j , k 

      IF ( idum.LE.0 .OR. iy.EQ.0 ) THEN
                                        ! Initialize
         idum = max(-idum,1)            ! Be sure to prevent idum=0
         DO 50 j = NTAB + 8 , 1 , -1    ! Load shuffle table after 8 warm-ups
            k = idum/IQ
            idum = IA*(idum-k*IQ) - IR*k
            IF ( idum.LT.0 ) idum = idum + IM
            IF ( j.LE.NTAB ) iv(j) = idum
 50      CONTINUE
         iy = iv(1)
      ENDIF
      k = idum/IQ                       ! Start here when not initializing
      idum = IA*(idum-k*IQ) - IR*k      ! Compute idum=mod(IA*idum,IM) without
      IF ( idum.LT.0 ) idum = idum + IM !    overflows by Schrage's Method
      j = 1 + iy/NDIV                   ! Will be in range 1:NTAB
      iy = iv(j)                        ! Output previously stored value
      iv(j) = idum                           ! Refill shuffle table
      T = AM*iy
      RAN1 = min(T,RNMX)        ! Because users don't expect endpoint values
c11   format(d90.50)
c     write(6,11)RAN1
      RETURN
      END FUNCTION RAN1

!  ---------------------------------------------------------------------
      integer function lenstr(string)
c     length of string

      implicit      none
c
      character(len=*), INTENT(IN) :: string   !input string
      character(len=1) :: null     !null character
c
      integer::     i,       !string index
     2              iend,    !candidate for string length
     3              j,       !index used for stepping through string
     4              length   !defined length of input
c
      length=len(string)
c
c     Look for the existence of nulls; if there
c     are any, the first one marks one character past the
c     end of the defined string. If there are no nulls,
c     pick up the final nonspace character, and call that
c     the defined length of the string.
c
      null=char(0)
c
c     Look for nulls
c
      iend=index(string(1:length),null)
c
c     No nulls, look for the location of the last nonspace.
c
      if(iend.eq.0)then
          do i=1,length
              j=length+1-i
              if(string(j:j).ne.' ')then
                  lenstr=j
                  return
              endif
          enddo
c
c         all spaces
c
          lenstr=0
          return
      else
c
c         first null-1
c
          lenstr=iend-1
          return
      endif
      end function lenstr
      
!  ---------------------------------------------------------------------
      SUBROUTINE ucase( string )
!
!  Shifts an arbitrary character string to UPPER CASE on any processor.
!
! From Stephen J. Chapman, "Fortran 90/95 for Scientist and Engineers, 
! Second Edition", McGraw Hill Higher Education, Boston, 2004), 
! pp. 434-435.
!
      CHARACTER(len=*), INTENT(INOUT) :: string
      INTEGER :: i, length

      length = LEN ( string )
   
      DO i = 1 , length
         IF ( LGE(string(i:i),'a') .AND. LLE(string(i:i),'z') ) THEN
            string(i:i) = ACHAR ( IACHAR ( string(i:i) ) - 32 )
         ENDIF
      END DO
      
      END SUBROUTINE ucase

!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------

      END MODULE utili_mod