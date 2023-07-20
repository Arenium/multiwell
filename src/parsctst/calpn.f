      SUBROUTINE CalPN(PN, DE, omega, xFF)
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
c
c         PN = Semi-classical tunneling transmission probability
c         DE = deltaE (Eq. 7b in Miller et al. 1990) ***with new sign convention***
c         omega = imaginary frequency, including coupling with orthogonal DOF: omega(F) in Eq. 7c (Miller et al. 1990)
c         xFF = diagonal anharmonic coupling coefficient for rxn coordinate
c
c         Latest revision: 12/2015
c
      Double precision PN, THETA, DE, omega, xFF
      Double precision PI, tp , D
      PI=acos(-1.0d0)

      D = (omega**2) / (4.0d0*xFF)                ! when xFF<0, -D = barrier height from VPT2
      IF ( omega.GT. 0.d+00 ) THEN
        IF( D.GE.0.0 .AND. DE.GE.D) THEN          ! failure of VPT2 when DE too positive 
          PN = 1.0d+00
        ELSEIF( D.LE.0.0 .AND. DE.LE.D) THEN      ! failure of VPT2 when DE too negative 
          PN = 0.0d+00
        ELSE
          tp = sqrt ( 1.0d0 - DE/D )
          THETA = -(PI*DE*2.0d0/omega)/( 1.0d0 + tp )      ! Eq. 7a in Miller et al. 1990 with new sign of DE
          PN = 1.0d0 / (1.0d0 + exp(2.0d0*THETA) )         ! Eq. 5b in Miller et al. 1990
        ENDIF
      ELSE
        IF (DE .GT. 0) THEN             ! logic: if omega .LE. 0, then barrier is very thick and hence no tunneling
          PN = 1.0d+00
        ELSE
          PN = 0.0d+00
        ENDIF
      ENDIF

      RETURN
      END 
