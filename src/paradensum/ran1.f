**==ran1.spg  processed by SPAG 5.11R  at 20:12 on 26 Mar 2001
 
 
      subroutine ran1(idum, iv, iy, randn)
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
C---------------------------------------------------------------------------------
C      "Minimal" random number generator of Park and Miller with Bays-Durham shuffle
C      and added safeguards.  Returns a uniform random deviate between 0.0 and 1.0
C      (exclusive of the endpoint values).  Call with idum a negative INTEGER(4)(KIND=4) to
C      initialize; thereafter do not alter idum between successive deviates in a
C      sequence.  RNMX should approximate the largest floating value that is less
C      than 1.
C
C      W. H. Press and S. A. Teukolsky, Computers in Physics, 6(5), 522-4 (1992).
C---------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(4) idum , IA , IM , IQ , IR , NTAB , NDIV
      REAL(KIND=8) AM , EPS , RNMX , T, randn
      PARAMETER (IA=16807,IM=2147483647,AM=1.0D+00/IM,IQ=127773,IR=2836,
     &        NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=2.3D-16,RNMX=1.0D+00-EPS)
      INTEGER(4) j , k , iv(NTAB) , iy 
!      SAVE iv , iy !si ricorda alla chiamata successiva dove era arrivato alla precedente nella sequenza
!      DATA iv/NTAB*0/ , iy/0/ !se è la prima chiamata iv e iy sono tutti zeri

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
      randn = min(T,RNMX)        ! Because users don't expect endpoint values
c11      format(d90.50)
c      write(6,11)RAN1
c      RETURN
      END SUBROUTINE ran1
