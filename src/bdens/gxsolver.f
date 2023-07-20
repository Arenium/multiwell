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
      REAL(8) FUNCTION gxsolver(C,zpe,DELE,Beta,ns,Gx,Emax)
!
!     gxsolver: a subroutine that solves G(x) for x.
!     		Two algorithms: one is recursive, like using excel work sheet.  
!		Use an initial value x, repeatedly calculating G(x), until x become stable.
!		The other one is a linear approach. Let x start from 1 and go to infinity with step size 1
!		For each x, calculate results comparing with G(x), until we arrive at a turnover point.
!  
      IMPLICIT NONE
      Real(8) C, zpe, DELE, Beta, Gx, Emax
!      Real(8) Ep, En        !Energy prev and Energy next for recursive way
      Real(8) Gp, Gn         !G(x) previous and G(x) next
      Real(8) Gr1, Gr2
      Real(8) ax, x
      Integer ns
      
      x = 1
      ax = 1.d0 - Beta*(LOG(5.0d0)-LOG(x/zpe))/19.d0
      Gp = (x + ax*zpe)**(ns+0.0d0)/C
      Gr1 = Gx - Gp
      
      x = x+1
      ax = 1.d0 - Beta*(LOG(5.0d0)-LOG(x/zpe))/19.d0
      Gn = (x + ax*zpe)**(ns+0.0d0)/C
      Gr2 = Gx - Gn
      DO WHILE((Gr1*Gr2).GT.0 .AND. x.LT.Emax)
        Gr1 = Gr2
        x = x+1
        ax = 1.d0 - Beta*(LOG(5.0d0)-LOG(x/zpe))/19.d0
        Gn = (x + ax*zpe)**(ns+0.0d0)/C
        Gr2 = Gx - Gn  
      END DO      
!      Recursive way
!      Ep = C - zpe
!      En = C - zpe + Beta*zpe*(LOG(5.0)-LOG(Ep/zpe))/19
      
!      DO WHILE(ABS(INT((Ep-En)/DELE)).GE.1)
!      Ep = En
!      En = C - zpe + Beta*zpe*(LOG(5.0)-LOG(Ep/zpe))/19
!      END DO
!      gxsolver = En
  
      gxsolver = x-1
      
      RETURN
      END FUNCTION gxsolver
 
