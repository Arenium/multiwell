c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Copyright (C) 2009 John R. Barker and Thanh Lam Nguyen
c
c Authors: Thanh Lam Nguyen and John R. Barker
c          nguyenlt@umich.edu
c          September 2009
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

        REAL(8) FUNCTION energy(N,freq,X,NN)

        INTEGER(4) No
        PARAMETER (No=100)        
c        INTEGER(4) I, J, K, L, N, NN(No)
        INTEGER(4) I, J, N, NN(No)
        REAL(8) eng, zpe
        REAL(8) freq(No), X(No,No)

        eng=0.0d0
        zpe=0.0d0
        DO I=1, N
                eng=eng + (NN(I)+0.5d0)*freq(I)
                zpe=zpe + 0.5d0*freq(I)
                DO J=1, I
                        eng=eng + (NN(I)+0.5d0)*(NN(J)+0.5d0)*X(I,J)
                        zpe=zpe + 0.25d0*X(I,J)
                ENDDO
        ENDDO   
        energy = eng - zpe
                
        IF(energy.LT.0.0d0) THEN
        IF(energy.GE.(-1.0d-5)) energy=0.0d0
        ENDIF
        
        RETURN
        END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


