       subroutine extractDataOUT(nameG03,countAtom,Multiplicity,
     & Stoichiometry,Ezpe,MolecularMass,ESN,Freq,DoF)

       integer maxAtom, correction
       PARAMETER (maxAtom=300)
       integer nul,countAtom,i,numFreq,DoF,ESN
       double precision Ezpe
       double precision Freq1,Freq2,Freq3,Freq(MaxAtom*3)
       double precision MolecularMass
       double precision Temperature, Pressure
       character*72 Find, nameG03
       character*15 FindFreq
       character*50 Stoichiometry
       integer numFile, Multiplicity
       save
           
      open(unit=3,file=nameG03(1:lnblnk(nameG03)),status='old')
      ESN=1
       
       ! Find multiplicity 
      do
       read(3,FMT='(A7)',end=10000) Find
!      write(*,*) Find
       if(Find.eq." Charge") then
       backspace(3)
       read(3,*,err=100) find,find,find,find,find, Multiplicity
       exit
       endif
100   CONTINUE
      enddo
       ! Find Stoichiometry
      do
       read(3,FMT='(A14)',end=10000) Find
       if(Find.eq." Stoichiometry") then
        backspace(3)
        read(3,*) find, Stoichiometry
        Stoichiometry=ADJUSTL(Stoichiometry)
        IF(Stoichiometry(len(TRIM(Stoichiometry)):
     &len(TRIM(Stoichiometry))).EQ.")") 
     &      Stoichiometry=Stoichiometry(1:len(TRIM(Stoichiometry))-3)
        exit
       endif
      enddo

       ! Find and count frequencies. 
      i=1
      Dof=0
      do
       read(3,'(A14)') FindFreq
       if (FindFreq.eq." Frequencies -") then
        Freq2=0
        Freq3=0
        backspace(3)
        read(3,*,err=200) FindFreq,FindFreq,Freq1,Freq2,Freq3
        goto 300
200     backspace(3)
        backspace(3)
        read(3,*) FindFreq,FindFreq,Freq1
      
300     Freq(i)=Freq1

        Dof=Dof+1
        if (Freq1.lt.0) then
         i=i-1
         Dof=Dof-1
        endif

        if (Freq2.ne.0) then
         Freq(i+1)=Freq2
         Dof=Dof+1
        endif
        if (Freq3.ne.0) then
         Freq(i+2)=Freq3
         Dof=Dof+1
         endif
         i=i+3
        endif
       If (FindFreq.eq." Temperature  ") exit
      enddo

      ! Find Temperature, Pressure, Mass,External Symmetry Number, E+zpe
      backspace(3)
      read(3,*) Find, Temperature, Find, Find, Pressure
      do
       read(3,'(A14)',end=30000) FindFreq
       If (FindFreq.eq." Molecular mas") then 
        backspace(3)
        read(3,*) Find, Find, MolecularMass
       endif
       If (FindFreq.eq." ROTATIONAL SY".OR.           
     &  FindFreq.eq." Rotational sy")then
        backspace(3)
        read(3,'(A27,A4)') Find,Find
        READ (Find(1:(lnblnk(Find)-1)),* )  ESN
       endif
       If (FindFreq.eq." Sum of electr") then
        backspace(3)
        read(3,FMT='(A51,F13.6)') Find, Ezpe
        exit
       endif
      enddo

      
      close(3)
      return

10000  write(*,*) "Impossible to find CHARGE and MULTIPLICITY or 
     & Stoichiometry"
       stop
!20000  write(*,*) "Impossible to find the frequencies"
!       stop
30000  write(*,*) "Impossible to find MASS, ESN and ZPE"
       stop
  
      end subroutine
