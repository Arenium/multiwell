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

        SUBROUTINE UHRLEV(T,AT,DELE,JMAX,NCB,NCV,CB,CV,n1,n2,IMAX,
     &          zpe,Vhr,Bhr,Phav,Phab,Nsig)

      IMPLICIT REAL(8)(A-H,O-Z), INTEGER(4)(I-N)
      PARAMETER (NN=2000)
      DIMENSION EV(NN), CV(NCV), CB(NCB), T(JMAX), AT(JMAX)
      DIMENSION IR(NN)
        CHARACTER Vhr*5, Bhr*5

        PI=acos(-1.0d0)
        Emax=(JMAX-1)*DELE

        IF((Bhr.EQ.'Bhrd2').OR.(Bhr.EQ.'BHRD2')) THEN   ! Average rotational constant
                B=0.0d0
                DO I=1, n2
                        B=B + (CB(I)/I)*((PI*2)**(I-1))
                ENDDO
        ELSEIF((Bhr.EQ.'Ihrd2').OR.(Bhr.EQ.'IHRD2')) THEN       ! Average moment of inertia
                B=0.0d0
                DO I=1, n2
                        B=B + (CB(I)/I)*((PI*2)**(I-1))
                ENDDO
                B=16.85763d0 / B        ! convert from I to B
        ELSEIF((Bhr.EQ.'Ihrd1').OR.(Bhr.EQ.'IHRD1')) THEN       ! Average moment of inertia
                B=16.85763d0/CB(1)      ! convert from I to B
        ELSEIF((Bhr.EQ.'Bhrd1').OR.(Bhr.EQ.'BHRD1')) THEN       ! Average rotational constant
                B=CB(1)
        ELSE
                write(*,*) "ERROR at INPUT for Bhr or Ihr"
        ENDIF

        CALL GHRLEV(EV,NN,IMAX,Emax,B,NCV,NCB,CV,CB,n1,n2,
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


