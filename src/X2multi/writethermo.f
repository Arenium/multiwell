      subroutine writethermo(nameG03,ZeroEnergy,Eunit,nTemp,TList,Punit)
       character*72 nameG03
       character*72 TList
       character*50 Stoichiometry
       character*4  Eunit, Punit
       character*4 STDSTATE
       integer nTemp,DofTot
       integer  countAtom, Multiplicity,DoF,ESN
       double precision Ix,Iy,Iz,Krot,ADrot,Ezpe
       double precision Freq(900), conver, Energy
       double precision MolecularMass,ZeroEnergy

c     Allowed standard states:
       IF ( PUNIT .NE. 'ATM' .OR. PUNIT .NE. 'atm' .OR. 
     &      PUNIT .NE. 'BAR' .OR. PUNIT .NE. 'bar' .OR. 
     &      PUNIT .NE. 'MCC' .OR. PUNIT .NE. 'mcc') THEN
         STDSTATE = 'BAR'
       ENDIF
       
10     conver=627.51
       if(Eunit.eq."KJOU") conver=2625.500
       if(Eunit.eq."CM-1") conver=219474.60

       open(unit=2,file="thermo.dat",status='unknown')
       write(2,*) Eunit, "  ",STDSTATE
       write(2,*) nTemp
       write(2,*) TList
       write(2,*) 1
       call extractDataMomInert(nameG03,countAtom,Ix,Iy,Iz)
       call calcRotor(Ix,Iy,Iz,Krot,ADrot)
       call extractDataOUT(nameG03,countAtom,Multiplicity,Stoichiometry,
     &  Ezpe, MolecularMass,ESN,Freq,Dof)

       Energy=(Ezpe-ZeroEnergy)*conver

       DoFTot=Dof+2              !  Vibrations + Krot + I2Drot
        
       if(Krot.eq.0)  DoFTot=DoFTot-1
       if(ADrot.eq.0)  DoFTot=DoFTot-1

       write(2,'(A6,A11,A1,F13.6)') "none  ",
     &  nameG03(1:lnblnk(nameG03)-4), " ", Energy
       write (2,*) Stoichiometry
       write (2,*) "! Comment line..."
       write (2,*) "! Up to 20 comment lines are allowed"
       write (2,*) "! Comment line..."
       write (2,'(I3,A9)') ESN, "   1   1 "
       write (2,*) "0.0 ", Multiplicity
       write (2,'(I4,A15)') DoFTot, " 'HAR'     AMUA"
       call writefreq(Dof,Freq,Krot,ADrot)
       close(2)
       return
      

1000   write(*,*) "Thermo units for set-up:   KCAL, 298 K, ATM"
       write(*,*) "   (manually change units, as appropriate)"
       Eunit="KCAL"
       nTemp=1
       TList="298 "
       Punit="ATM"
       goto 10
       end subroutine
