c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c UHRLEV: a code to calculate eigenvalues of unsymmetrical, non-rigid 1D-hindered internal rotation
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
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE UHRLEV(Emax,EV,NCB,NCV,CB,CV,n1,n2,NIMAX,
     &            zpe,Vhr,Bhr,Phav,Phab,Vo,Vm)

      REAL(8) EV, CV, CB, PI,Phav,Phab,Vo,Vm,Emax,B,zpe
      INTEGER(4) I,n2,NCV,NCB,NN,NIMAX
      PARAMETER (NN=2000)
      DIMENSION EV(NN), CV(NCV), CB(NCB) 
      CHARACTER(len=5) Vhr , Bhr
      SAVE

      PI=acos(-1.0d0)

      IF( Bhr.EQ.'BHRD2' ) THEN      ! Average rotational constant
            B=0.0d0
            DO I=1, n2
                  B=B + (CB(I)/I)*((PI*2)**(I-1))
            ENDDO
      ELSEIF( Bhr.EQ.'IHRD2' ) THEN      ! Average moment of inertia
                B=0.0d0
                DO I=1, n2
                        B=B + (CB(I)/I)*((PI*2)**(I-1))
                ENDDO
            B=16.85763d0 / B      ! convert from I to B
        ELSEIF( Bhr.EQ.'IHRD1' ) THEN      ! Average moment of inertia
                B=16.85763d0/CB(1)      ! convert from I to B
        ELSEIF( Bhr.EQ.'BHRD1' ) THEN      ! Average rotational constant                  
            B=CB(1)
      ELSE
            write(*,*) "ERROR at INPUT for Bhr or Ihr"
      ENDIF

      CALL GHRLEV(EV,NN,NIMAX,Emax,B,NCV,NCB,CV,CB,n1,n2,
     &            Vhr,Bhr,Phav,Phab,Vo,Vm)

      zpe=EV(1)      ! zero-point energy (cm-1)

      RETURN
      END


