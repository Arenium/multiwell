c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c SHRLEV: a code to calculate eigenvalues of symmetrical, rigid 1D-hindered internal rotation
c Copyright (C) 2009 Lam T. Nguyen and John R. Barker
c
c      DATE: Feb. 13, 2009
c
c John R. Barker
c jrbarker@umich.edu
c University of Michigan
c Ann Arbor, MI 48109-2143
c (734) 763 6239
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE SHRLEV(Emax,EV,B,V,n1,n2,NIMAX,zpe,
     &      Vhr,Bhr,Phav,Phab,Vo,Vm)

      IMPLICIT NONE
      INTEGER(4) n1, n2, NIMAX, NN, No, NCV, NCF
      PARAMETER (NN=2000, No=1)
      REAL(8) Emax, EV, B,V, zpe, Phav, Phab, Vo, Vm
      REAL(8) CV, CF
      DIMENSION EV(NN), CV(No), CF(No) 
      CHARACTER(5) Vhr , Bhr
      SAVE

      NCV=1
      NCF=1
      CV(1)=V
      CF(1)=B

      CALL GHRLEV(EV,NN,NIMAX,Emax,B,NCV,NCF,CV,CF,n1,n2,
     &      Vhr,Bhr,Phav,Phab,Vo,Vm)

      zpe=EV(1)     ! zero-point energy (cm-1)

      RETURN
      END


