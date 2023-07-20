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
        REAL(8) FUNCTION energyxyz(N,freq,X,Y,Z,NN)

        IMPLICIT NONE
        INTEGER(4) Num
        PARAMETER (Num=100)        
        INTEGER(4) I, J, K, L, N, NN(Num)
        REAL(8) eng, zpe
        REAL(8) freq(Num), X(Num,Num)
        REAL(8) Y(Num,Num,Num), Z(Num,Num,Num,Num) 

        eng=0.0D0
        zpe=0.0d0
        DO I=1, N
                eng=eng + (NN(I)+0.5d0)*freq(I)
                zpe=zpe + 0.5d0*freq(I)
                DO J=1, I
                        eng=eng + (NN(I)+0.5d0)*(NN(J)+0.5d0)*X(I,J)
                        zpe=zpe + 0.25d0*X(I,J)
                DO K=1, J
                        eng=eng+(NN(I)+0.5d0)*
     &                  (NN(J)+0.5d0)*(NN(K)+0.5d0)*Y(I,J,K)
                        zpe=zpe + 0.125d0*Y(I,J,K)
                DO L=1, K 
                        eng=eng+(NN(I)+0.5d0)*(NN(J)+0.5d0)*
     &                          (NN(K)+0.5d0)*(NN(L)+0.5d0)*Z(I,J,K,L)
                        zpe=zpe + 0.0625d0*Z(I,J,K,L)
                ENDDO
                ENDDO
                ENDDO
        ENDDO   
        energyxyz = eng - zpe
                
        IF(energyxyz.LT.0.0d0) THEN
        IF(energyxyz.GE.(-1.0d-4)) energyxyz=0.0d0
        ENDIF
        
        RETURN
        END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


