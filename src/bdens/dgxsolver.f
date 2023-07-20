!    
!    LICENSE NOTICE
!
!    bdens: sums and densities of coupled anharmonic vibrations
!    Copyright (C) 2015 Collin G. L. Li, Thanh Lam Nguyen, and John R. Barker
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License: <http://www.gnu.org/licenses/>.
!
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c                                                                 
c    PROGRAM bdens
c    by Collin G. L. Li, Thanh Lam Nguyen, and John R. Barker    
c                                                                  
c  ***Direct Count and Wang-Landau Algorithm for Densities of States***
c
c    Contact:
c    John R. Barker   (email: jrbarker@umich.edu)
c    Department of Atmospheric, Oceanic, and Space Sciences
c    College of Engineering
c    University of Michigan
c    Ann Arbor, MI 48109
c
c    http://aoss.engin.umich.edu/multiwell
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-----------------------------------------------------------------------------------------
      FUNCTION dgxsolver(Emax,C,zpe,Be,s,alpha, b1,Beta, KEYWORD)
!
!     dgxsolver: a function that solves G(x) for x while considering the effect of Wang-Landou programs
!		 First, it decides the accuracy of Wang-Landou program.  Then, it tried to find the optimal
!		 energy point combining Wang-Landou program and G(x) using linear search approach. Let x 
! 		 start from 1 and go to infinity with step size 1.  For each x, calculate results comparing, 
!		 until we arrive at a turnover point.
!
      IMPLICIT NONE
      REAL(8) :: dgxsolver
      Real(8), INTENT(IN) :: Emax, C, zpe, Be, alpha, b1, Beta
      CHARACTER(6), INTENT(IN) :: KEYWORD              ! designates no. of Wang-Landau stochastic trials; assumed to be upper case
      Real(8) :: Ep, En, Er1, Er2                      !Ep: energy previous, En:energy next
      Integer :: s, x
      Real(8) :: ax, kindex
      
      IF    ( KEYWORD .EQ. 'EXTRA' ) THEN
        kindex = 1.d+1
      ELSEIF ( KEYWORD .EQ. 'BEST' ) THEN
        kindex = 1.0d0
      ELSEIF ( KEYWORD .EQ. 'BETTER' ) THEN
        kindex = 1.d-1
      ELSEIF ( KEYWORD .EQ. 'GOOD' ) THEN
        kindex = 1.d-2
      ELSEIF ( KEYWORD .EQ. 'FAIR' ) THEN
        kindex = 1.d-3
      ELSE
        kindex = 1.d0              !default set as 'best'
      ENDIF
      
      x = 1
      ax = 1.d0 - Beta*(LOG(5.0d0)-LOG(x/zpe))/19.d0
      Ep = C*(x+zpe*ax)**(s*alpha-1.d0)*(1.d0+zpe*(Beta/(19.d0*x)))
      En = kindex*(Emax - x)**(Be-1.0d0)
      Er1 = ABS(Ep-En)
      
      x = x+1
      ax = 1.d0 - Beta*(LOG(5.0d0)-LOG(x/zpe))/19.d0      
      Ep = C*(x+zpe*ax)**(s*alpha-1.d0)*(1.d0+zpe*(Beta/(19.d0*x)))
      En = kindex*(Emax - x)**(Be-1.0d0)
      Er2 = ABS(Ep-En)
      DO WHILE(Er1.GT.Er2 .AND. x.LT.Emax)
        x = x+1
        Er1 = Er2
        ax = 1.d0 - Beta*(LOG(5.0d0)-LOG(x/zpe))/19.d0
        Ep = C*(x+zpe*ax)**(s*alpha-1.d0)*(1.d0+zpe*(Beta/(19.d0*x)))
        En = kindex*(Emax - x)**(Be-1.0d0)  
        Er2 = ABS(Ep-En)
      ENDDO            
      dgxsolver = x
      
      RETURN
      END FUNCTION dgxsolver
      
