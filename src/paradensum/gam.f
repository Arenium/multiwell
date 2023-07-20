**==gam.spg  processed by SPAG 5.11R  at 20:12 on 26 Mar 2001
      DOUBLE PRECISION FUNCTION GAM(X)
 
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
!                         PROGRAM paradensum    
!
!                               by      
!
!           by Michele Ceotto, Chiara Aieta, Fabio Gabas,
!               Thanh Lam Nguyen, and John R. Barker    
!                                                                  
!           ***PARallel Anharmonic DENsities and SUM of states***
!
!                             based on
!
!          The parallel implementation* of the Wang-Landau algorithms 
!          for densities of states**
!
!    Literature Citations:                                                    
!    *parallel implementation of anaharmonic density of states
!    C. Aieta, F. Gabas and M. Ceotto, J. Phys. Chem. A., 120(27), 
!    4853-4862 (2016).
!                                                                  
!    **Density of states algorithms
!    F. Wang and D. P. Landau, Phys. Rev. Letters 86, 2050-2053 (2001).           
!    M. Basire, P. Parneix, and F. Calvo, J. Chem.Phys. 129, 081101 (2008).    
!    T. L. Nguyen and J.R. Barker, J. Phys. Chem. A., 114, 3718–3730 (2010).
!    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
c
c    gamma function of X, from FUNCTION ALNGAM
c
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      EXTERNAL ALNGAM
!      SAVE 
      XX = X
      GAM = EXP(ALNGAM(XX))
 100  RETURN
      END
