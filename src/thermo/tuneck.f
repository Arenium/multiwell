c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Thermo: a code for thermochemical calculations.
c Copyright (C) 2010 John R. Barker
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
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION tuneck(T,vimag,Vf,Vr)
c ----------------------------------------------------------------------
c      Tunneling via unsymetrical Eckart barrier
c
c      T      = temperature
c      vimag  = magnitude of imaginary frequency (cm-1)
c      Vf     = Barrier height, forward direction (cm-1)
c      Vr     = Barrier height, reverse direction (cm-1)
c
c      C      = heat capacity at constant pressure (units of R)
c      S      = entropy (units of R)
c      H      = [H¡(T)-H¡(0)]/RT
c
c      tuneck = transmission coefficient
c    
c    Unsymmetrical Eckart potential: formulae from 
c       W. Forst, Unimolecular Reactions - A Concise Introduction 
c                 (Cambridge University Press, 2003), p. 172ff.
c ----------------------------------------------------------------------
!
!	2020-6-14 Bug reported by and fixed by Gabriel R. Da Silva and 
!              Milad Marimani (Univ. of Melbourne, Australia)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8) tuneck
      REAL(8) Q , T , vimag , Vf , Vr , dE , E
      REAL(8) V1 , V2 , aa , bb , cc , dd , arg , TERM , uu , term1 , 
     &        term2 , TEST
      REAL(8) p , A , B, Estart
      REAL(8) hor , PI , TWOPI
      PARAMETER (hor = 1.4387752D+00)
      SAVE 
      
      PI = 4.0D+00*DATAN(1.0D+00)
c      PI = 3.14159265358979d+00
      TWOPI = 2.d+00*PI
      
        V2 = Vr
        V1 = Vf

      IF ( V2 .GT. V1 ) THEN
        Estart = 0.0d+00
       ELSE 
        Estart = V1-V2
      ENDIF

      A = ( V1 - V2 )
      B = ( SQRT(V1) + SQRT(V2) )**2
      cc = ( A-B )*( A+B )/( 2.d+00*vimag*SQRT(B**3) )
      arg = 4.d+00*B*cc**2 - 1.d+00
      dd = PI*SQRT( ABS(arg) )
       
      dE = T/(hor*100.d+00)               ! units of cm-1

      tuneck = 0.0d+00
      Q = 0.d+00
      p = 0.d+00
      E = Estart                         					! units of cm-1
      DO WHILE ( p .LT. 0.9999d+00 )
         E = E + dE                      					! Translational energy, units of cm-1
         aa = TWOPI*cc*sqrt(E)
         bb = TWOPI*cc*SQRT( E - A )
         uu = aa + bb
         IF ( arg .GT. 0.d+00 ) THEN
c           p = 2.d+00*SINH(aa)*SINH(bb)/( COSH(uu) + COSH(dd) )          ! Forst (2003) Eq. 6.6
           term1 = -exp(-2.d+00*aa) - exp(-2.d+00*bb) + exp(-2.d+00*uu)  ! Forst (2003) Eq. 6.8
           term2 = exp(-2.d+00*uu) + exp(-(uu-dd)) + exp(-(uu+dd))       ! Forst (2003) Eq. 6.8
           p = (1.d+00 + term1)/(1.d+00 + term2)                         ! Forst (2003) Eq. 6.8
          ELSE
           p = 2.d+00*SINH(aa)*SINH(bb)/( COSH(aa+bb) + COS(dd) )        ! Forst (2003) Eq. 6.7
         ENDIF

         TERM = dE*p*exp(-E*hor/T)*hor/T
         Q = Q + TERM

       END DO
       
       p = 1.0d+00			           ! assume p=1.00 at higher energies until Emax
       TEST = 1.0
       DO WHILE ( TEST .GT. 1.0d-08 )
           E = E + dE                      ! units of cm-1
           TERM = dE*p*exp(-E*hor/T)*hor/T
           Q = Q + TERM
           TEST = TERM/Q
       END DO
       
       tuneck = Q*EXP(Vf*hor/T)   ! Partition function for Transmission coefficient
      
      RETURN
 
      END FUNCTION
