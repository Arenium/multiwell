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
!       SUBROUTINE UHRLEV(Ttot,AT,DELE,JMAX,NCB,NCV,CB,CV,n1,n2,IMAX,
!     &          zpe,Vhr,Bhr,Phav,Phab,Nsig)
       SUBROUTINE UHRLEV(DELE,JMAX,NCB,NCV,CB,CV,n1,n2,IMAX_1,
     &          zpe_1,Vhr,Bhr,Phav_1,Phab_1,Nsig_1)
   
      USE decl_alloc

      IMPLICIT REAL(8)(A-H,O-Z), INTEGER(4)(I-N)
      PARAMETER (NN=2000)
      !DIMENSION EV(NN), CV(NCV), CB(NCB), T(JMAX), AT(JMAX)
      DIMENSION EV_2(NN), CV(NCV), CB(NCB)
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

        CALL GHRLEV(EV_2,NN,IMAX_1,Emax,B,NCV,NCB,CV,CB,n1,n2,
     &          Vhr,Bhr,Phav_1,Phab_1)

      zpe_1=EV_2(1)

C......      write(*,*) "Zero-point energy = ", zpe

      DO I=1, IMAX_1
            tp=EV_2(I)-zpe_1
            IR(I)=NINT(tp/DELE)
C......            write(*,*) I, tp
      ENDDO 

      DO J=1, IMAX_1
            LL=IR(J)
            DO L=1, JMAX - LL
                  AT(L+LL) = AT(L+LL) + Ttot(L)
            ENDDO
      ENDDO
      DO J=1, JMAX
            Ttot(J)=AT(J)/Nsig_1
            AT(J)=0.0d0
      ENDDO 

      RETURN
      END


