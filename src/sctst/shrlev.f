c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c SHRLEV: a code to calculate eigenvalues of symmetrical, rigid 1D-hindered internal rotation
c Copyright (C) 2009 Lam T. Nguyen and John R. Barker
c
c      DATE:      Feb. 13, 2009
c
c John R. Barker
c jrbarker@umich.edu
c University of Michigan
c Ann Arbor, MI 48109-2143
c (734) 763 6239
c
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
c See the "ReadMe" file for a copy of the GNU General Public License,
c or contact:
c
c Free Software Foundation, Inc.
c 59 Temple Place - Suite 330
c Boston, MA 02111-1307, USA.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE SHRLEV(T,AT,DELE,JMAX,B,V,n1,n2,IMAX,zpe,
     &          Vhr,Bhr,Phav,Phab,Nsig)

        IMPLICIT REAL(8)(A-H,O-Z), INTEGER(4)(I-N)
        PARAMETER (NN=2000, No=1)
        DIMENSION EV(NN), CV(No), CF(No), T(JMAX), AT(JMAX)
         DIMENSION IR(NN)
        CHARACTER Vhr*5, Bhr*5

        Emax=(JMAX-1)*DELE

        NCV=1
        NCF=1
        CV(1)=V
        CF(1)=B

        CALL GHRLEV(EV,NN,IMAX,Emax,B,NCV,NCF,CV,CF,n1,n2,
     &          Vhr,Bhr,Phav,Phab)

      zpe=EV(1)

C......      write(*,*) "Zero-point energy = ", zpe

      DO I=1, IMAX
            tp=EV(I)-zpe
            IR(I)=NINT(tp/DELE)
C......            write(*,*) I, tp

      ENDDO 

      DO J=1, IMAX
            LL=IR(J)
            DO L=1, JMAX - LL
                  AT(L+LL) = AT(L+LL) + T(L)
            ENDDO
      ENDDO
      DO J=1, JMAX
            T(J)=AT(J)/Nsig
            AT(J)=0.0d0
      ENDDO 

      RETURN
      END


