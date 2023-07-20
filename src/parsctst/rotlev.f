**==rotlev.spg  processed by SPAG 5.11R  at 20:12 on 26 Mar 2001
 
 
      !SUBROUTINE ROTLEV(Ttot,AT,DELE,JMAX,B,NDIM,NSYMM)
      SUBROUTINE ROTLEV(DELE,JMAX,B,NDIM,NSYMM)
 
      USE decl_alloc
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!    
!    LICENSE NOTICE
!
!    Copyright (C) 2019 Michele Ceotto, Chiara Aieta, Fabio Gabas, 
!                       Thanh Lam Nguyen, and John R. Barker
!
!    Contact:
!
!    Michele Ceotto  (email: michele.ceotto@unimi.it)
!    Dipartimento di Chimica
!    Università degli Studi di Milano
!    via Golgi 19, Milano 20133, ITALY
!
!    or:
!
!    John R. Barker   (email: jrbarker@umich.edu)
!    Department of Atmospheric, Oceanic, and Space Sciences
!    College of Engineering
!    University of Michigan
!    Ann Arbor, MI 48109
!
!    http://clasp-research.engin.umich.edu/multiwell
!    This program is free software; you can redistribute it and/or
!    modify it under the terms of the GNU General Public License (version 2)
!    as published by the Free Software Foundation.
!   
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!   
!    See the 'ReadMe' file for a copy of the GNU General Public License,
!    or contact:
!   
!    Free Software Foundation, Inc.
!    59 Temple Place - Suite 330
!    Boston, MA 02111-1307, USA.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                 
!                         PROGRAM parsctst    
!
!                               by      
!
!           by Michele Ceotto, Chiara Aieta, Fabio Gabas,
!               Thanh Lam Nguyen, and John R. Barker    
!                                                                  
!           ***PARalle Semi-Classical Transition State Theory***
!
!                             based on
!
!          The theory of W. H. Miller and coworkers* and the
!          parallel implementation** of the Wang-Landau algorithms 
!          for densities of states***
!
!    Literature Citations:                                                    
!    *Semi-Classical Transition State Theory
!    W. H. Miller, J. Chem. Phys. 62, 1899-1906 (1975).
!    W. H. Miller, Faraday Discuss. Chem. Soc. 62, 40-46 (1977).
!    W. H. Miller, R. Hernandez, N. C. Handy, D. Jayatilaka, and A. Willets,
!      Chem. Phys. Letters 172, 62-68 (1990).
!    R. Hernandez and W. H. Miller, Chem. Phys. Lett. 214, 129-136 (1993).
!    J. F. Stanton, J. Phys. Chem. Lett. 7, 2708-2713 (2016)
!
!    **Semi-Classical Transition State Theory parallel implementation
!    C. Aieta, F. Gabas and M. Ceotto, J. Chem. Theory Comput.,
!      15, 2142−2153 (2019).
!                                                                  
!    ***Density of states algorithms
!    F. Wang and D. P. Landau, Phys. Rev. Letters 86, 2050-2053 (2001).           
!    M. Basire, P. Parneix, and F. Calvo, J. Chem.Phys. 129, 081101 (2008).    
!    T. L. Nguyen and J.R. Barker, J. Phys. Chem. A., 114, 3718–3730 (2010).
!    C. Aieta, F. Gabas and M. Ceotto, J. Phys. Chem. A., 120(27), 4853-4862 (2016).
!                                                                  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C      CALCULATES ENERGY LEVELS FOR FREE ROTOR TO BE USED
C          WITH STERAB
C
c      DELE      = energy grain size
c      JMAX      = index of Emax
c      B         = rotational constant (cm-1)
c      NDIM      = rotor dimension
c      IR        = vector of energy level indices
c      NMAX      = ceiling of array
c
c
      IMPLICIT REAL(8)(A-H,O-Z)
      DIMENSION IR(20001)
      !DIMENSION T(JMAX), AT(JMAX)
!      SAVE 

      Emax=(JMAX-1)*DELE 
      RMAX = SQRT(Emax/B)  ! Maximum bound level
      NMAX = INT(RMAX)  
      IF ( NMAX.GT.20001 ) NMAX = 20001
 
      DO 100 J = 1 , NMAX
         IF ( NDIM.EQ.1 ) THEN
             R = B*J*J/DELE
           ELSE
             R = B*J*(J+1)/DELE
         ENDIF
         IR(J) = NINT(R)              ! Nearest INTEGER(4): number of grains
 100  CONTINUE

CCCCCCCCCCCCCCCCCCCCCCCCC

      DO J = 1 , NMAX ! J = rot quantum number;  NMAX is the number of energy levels
        IF ( J.EQ.0 ) THEN
           F = 1d+00             ! F is the multiplicity of the rotational level (see J.L. McHale, Molecular Spectroscopy (Prentice-Hall, 1999), 216f.
        ELSEIF ( NDIM.EQ.1 ) THEN   ! 1-D rotor (e.g. single-axis internal rotor)
           F = 2.0d+00
        ELSEIF ( NDIM.EQ.2 ) THEN   ! 2-D rotor
           F = 2.0d+00*J + 1.0d+00
        ELSEIF ( NDIM.EQ.3 ) THEN   ! 3-D rotor (shperical top)
           F = (2.0d+00*J+1.0d+00)*(2.0d+00*J+1.0d+00)
        ENDIF
        DO K = 1 , JMAX
           KARG = K - IR(J)
           IF ( KARG.GT.0 .AND. KARG.LE.JMAX) THEN
               AT(K) = AT(K) + F*Ttot(KARG)
           ENDIF
        END DO
      END DO

      DO K = 1 , JMAX
          Ttot(K) = AT(K) + Ttot(K)
          Ttot(K) = Ttot(K) / NSYMM
          AT(K) = 0.0D+00
      END DO

CCCCCCCCCCCCCCCCCCCCCCCC
 
      RETURN
      END
