
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Copyright (C) 2009 John R. Barker
c
c Authors: Thanh Lam Nguyen and John R. Barker
c          nguyenlt@umich.edu   
c          Aug. 9, 2009
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
c See the 'ReadMe' file for a copy of the GNU General Public License,
c or contact:
c
c Free Software Foundation, Inc.
c 59 Temple Place - Suite 330
c Boston, MA 02111-1307, USA.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        Double precision FUNCTION energyxyzdoloops(freq,N,NN,X,Y,Z)
        
        Integer I, J, K, L, N, NN(100)
        Real(8) eng, zpe
        Real(8) freq(100), X(100,100)
        Real(8) Y(100,100,100), Z(100,100,100,100) 

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
                        zpe=zpe + 0.5d0*0.5d0*0.5d0*Y(I,J,K)
                DO L=1, K 
                        eng=eng+(NN(I)+0.5d0)*(NN(J)+0.5d0)*
     &                          (NN(K)+0.5d0)*(NN(L)+0.5d0)*Z(I,J,K,L)
                        zpe=zpe + 0.5d0*0.5d0*0.5d0*0.5d0*Z(I,J,K,L)
                ENDDO
                ENDDO
                ENDDO
        ENDDO   
        energyxyzdoloops = eng - zpe
                
        IF(energyxyzdoloops.LT.0.0d0) THEN
        IF(energyxyzdoloops.GE.(-1.0d-4)) energyxyzdoloops=0.0d0
        ENDIF
        
        RETURN
        END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


