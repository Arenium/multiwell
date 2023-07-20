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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c To compute the first derivative of Ev energy function 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE CKDERIVXYZ(ns, ntest)
    
      USE decl_alloc

      IMPLICIT NONE
!      INTEGER(4) ns , i , j , k , l, ntest , nv(100)
      INTEGER(4) ns , i , j , k , l, ntest 
c      REAL(8) wa(100) , xa(100,100), deriv , sum 
      REAL(8) deriv 
!      REAL(8) ya(100,100,100) , za(100,100,100,100) 

        REAL(8) XX, A1, A2, A3, A4, tp

!      SAVE


      ntest = 0
      k = 1
      DO WHILE ( (k .LE. ns) .AND. (ntest .EQ. 0) )    ! start checking

        A1=4.0d0*za(k,k,k,k)

        A2=ya(k,k,k)
        DO j=1, ns
                IF(j.ne.k) A2=A2+za(k,k,k,j)*(nv(j)+0.5d0)
        ENDDO
        A2=3.0d0*A2

        A3=xa(k,k)
        DO j=1, ns
                IF(j.ne.k) A3=A3+ya(k,k,j)*(nv(j)+0.5d0)
        ENDDO
        DO i=1, ns
                IF(i.ne.k) THEN
                        DO j=1, i
        IF(j.ne.k) A3=A3+za(k,k,i,j)*(nv(i)+0.5d0)*(nv(j)+0.5d0)
                        ENDDO
                ENDIF
        ENDDO
        A3=2.0d0*A3

        A4=wa(k)
        DO i=1, ns
                IF(i.ne.k) A4=A4+xa(k,i)*(nv(i)+0.5d0)
        ENDDO
        DO i=1, ns
                IF(i.ne.k) THEN
                        DO j=1, i
        IF(j.ne.k) A4=A4+ya(k,i,j)*(nv(i)+0.5d0)*(nv(j)+0.5d0)
                        ENDDO
                ENDIF
        ENDDO
        DO i=1, ns
                IF(i.ne.k) THEN
                DO j=1, i
                        IF(j.ne.k) THEN
                                DO l=1, j
                                IF(l.ne.k) THEN
        tp=(nv(i)+0.5d0)*(nv(j)+0.5d0)*(nv(l)+0.5d0)
        A4=A4+za(k,i,j,l)*tp
                                ENDIF
                                ENDDO
                        ENDIF
                ENDDO
                ENDIF
        ENDDO

        XX = 0.5d0 + nv(k) 
        deriv = A1*(XX**3) + A2*(XX**2) + A3*XX + A4

         IF ( deriv .LE. 0.0d+00 ) ntest = 1           ! failed test
         k = k + 1
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


