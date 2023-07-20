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
      SUBROUTINE rotunits( VROTIN , X , VROTOUT )
      
      IMPLICIT NONE
      CHARACTER(len=4) VROTIN , VROTOUT     
      REAL(8) X
      SAVE
     
c      CONVERT UNITS FROM VROTIN TO VROTOUT
c
c       VROTIN or VROTOUT
c               = 'AMUA' for moment of inertia in units of amu*Ang^2
c       or      = 'GMCM' for moment of inertia units of gram*cm^2
c       or      = 'CM-1' for rotational constant in units of cm^-1
c       or      = 'MHZ' for rotational constant in units of MHz
c       or      = 'GHZ' for rotational constant in units of GHz
c         
c         Convert VROTIN to AMUA
c
      IF ( VROTIN .EQ. 'AMUA' .OR. VROTIN .EQ. 'amua' ) THEN
      ELSEIF ( VROTIN .EQ. 'GMCM' .OR. VROTIN .EQ. 'gmcm' ) THEN
        X = X / 1.660538782d-040
      ELSEIF ( VROTIN .EQ. 'CM-1' .OR. VROTIN .EQ. 'cm-1' ) THEN
        X = 16.85763D+00 / X
      ELSEIF ( VROTIN .EQ. 'MHZ' .OR. VROTIN .EQ. 'MHz' 
     &      .OR. VROTIN .EQ. 'mhz' .OR. VROTIN .EQ. 'Mhz' ) THEN
        X = 5.05379D+005 / X
      ELSEIF ( VROTIN .EQ. 'GHZ' .OR. VROTIN .EQ. 'GHz' 
     &      .OR. VROTIN .EQ. 'ghz' .OR. VROTIN .EQ. 'Ghz' ) THEN
        X = 5.05379D+002 / X
      ELSE
        write (*,*) 'FATAL: units (VROTIN) for rotations not recongized'
        STOP
      ENDIF
c
c         Convert AMUA to VROTOUT and RETURN
c
      IF ( VROTOUT .EQ. 'AMUA' .OR. VROTOUT .EQ. 'amua' ) THEN
        RETURN
      ELSEIF ( VROTOUT .EQ. 'GMCM' .OR. VROTOUT .EQ. 'gmcm' ) THEN
        X = X * 1.660538782d-040
        RETURN
      ELSEIF ( VROTOUT .EQ. 'CM-1' .OR. VROTOUT .EQ. 'cm-1' ) THEN
        X = 16.85763D+00 / X
        RETURN
      ELSEIF ( VROTOUT .EQ. 'MHZ' .OR. VROTOUT .EQ. 'MHz' 
     &      .OR. VROTOUT .EQ. 'mhz' .OR. VROTOUT .EQ. 'Mhz' ) THEN
        X = 5.05379D+005 / X
        RETURN
      ELSEIF ( VROTOUT .EQ. 'GHZ' .OR. VROTOUT .EQ. 'GHz' 
     &      .OR. VROTOUT .EQ. 'ghz' .OR. VROTOUT .EQ. 'Ghz' ) THEN
        X = 5.05379D+002 / X
        RETURN
      ELSE
        write (*,*) 'FATAL: rotation units (VROTOUT) not recongized'
        STOP
      ENDIF

      END
      
      
