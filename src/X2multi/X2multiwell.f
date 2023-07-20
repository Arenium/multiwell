       integer i,zero, nTemp
       character*72 nameG03, namefile, Find, nameGausslog
       integer numFile, Multiplicity
       character*4 Eunit , Punit, TypeInter ,namepath
       character*5 OStext,OStext2
       double precision conver, ZeroEnergy ,Energy
       character*72 multiDist,MULTIWELL, TList, Elimit
       character*50 Stoichiometry

!       call getenv('MULTIWELL',multiDist)    
       call OpSystem(OStext,OStext2)   
       
       call system (OStext2 //" energies")
       
       call GETARG(1,nameGausslog)
       if (nameGausslog.ne.'') then
         call singlefile(nameGausslog,OStext)
         stop
       endif

       call writeCFG()

       open(unit=4,file='gauss2multi.cfg',status='old')
       call Find1stwell(ZeroEnergy,OStext)  
       
       read(4,*) Eunit
       read(4,*) nTemp
       read(4,'(A72)') TList
       read(4,*) Punit
       read(4,*) 
       read(4,*) 
       read(4,'(A72)') Elimit
       
       zero=0 
       i=0
       
       do
       i=i+1
       read(4,*,end=10000,err=9000) numFile, nameG03, TypeInter
       write(*,'(A19,T21,A20)') " Processing file ",
     &nameG03(1:lnblnk(nameG03)-4)
       namefile= nameG03(1:lnblnk(nameG03)-4)
!---------------------------------------------------------------------------
!      create Mominert data file and run it
      call writemominert(nameG03)
!     MULTIWELL=multiDist(1:lnblnk(multiDist))//"bin/mominert"
!     call system (MULTIWELL)
      call system ("mominert")

      find=OStext//" mominert.dat "//namefile(1:(lnblnk(namefile)))//
     &".coords"
      call system  (Find)
      find=OStext//" mominert.out "//namefile(1:(lnblnk(namefile)))//
     &".coords.out"
      call system (Find)
       
!---------------------------------------------------------------------------
!      create densum and run it
       call writedensum(nameG03,Elimit)
!      MULTIWELL=multiDist(1:lnblnk(multiDist))//"bin/densum > NULL"
!      call system (MULTIWELL)
       call system ("densum > NULL")
       find=OStext//" densum.dat "//namefile(1:(lnblnk(namefile)))//
     &".vibs"
       call system (Find)
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
!      create thermo and run it
       call writethermo(nameG03,ZeroEnergy,Eunit,nTemp,TList,Punit)
!      MULTIWELL=multiDist(1:lnblnk(multiDist))//"bin/mominert"
!      call system (MULTIWELL)
       call system ("thermo")
       find=OStext//" thermo.dat "//nameG03(1:lnblnk(nameG03)-4)//
     & ".therm"
        call system (find)
       find=OStext//" thermo.out "//nameG03(1:lnblnk(nameG03)-4)//
     & ".therm.out"
        call system (find)

!---------------------------------------------------------------------------
       enddo

9000   write(*,*) "Error reading gauss2multi.cfg file."
       close(4)
!      PAUSE ' STOP'
       stop  'STOP'

10000  close(4)
!---------------------------------------------------------------------------
!      create multiwell.dat
        call writemultiwelldat
        write(*,*) "  DONE"
        write(*,*) 
        write(*,*) "  multiwell.dat template created " 
        write(*,*)
!---------------------------------------------------------------------------
      write(*,*)
      write(*,*) "------------------------------------------------------
     &-----------"
      write(*,*)
      write(*,*) " If you need to add internal rotors, you have to MANUA
     &LLY modify"
      write(*,*) " the .coords and the .vibs files and run them again"
      write(*,*)
      write(*,*) "------------------------------------------------------
     &-----------"
      write(*,*)
      write(*,*) " Please check the number of optical isomers and"
      write(*,*) " electronic level degeneracy in .therm files are corre
     &ct."
      write(*,*) " They cannot be determined automatically!"
      write(*,*)

       end




      subroutine Find1stWell(ZeroEnergy,OStext)
      double precision ZeroEnergy
      integer maxAtom, correction
      PARAMETER (maxAtom=300)
      integer nul,countAtom,i,numFreq,DoF,ESN
      double precision Ezpe
      double precision Freq1,Freq2,Freq3,Freq(MaxAtom*3)
      double precision MolecularMass
      character*72 Find, nameG03
      character*15 FindFreq
      character*50  Stoichiometry
      integer numFile, Multiplicity
      character*4  TypeInter
      character*5  OStext

     
      do i=1,7
      read(4,*) 
      enddo
    
      do
       read(4,*,end=20000,err=19000) numFile, nameG03, TypeInter 
       IF(TypeInter.eq."WELL") then
        call extractDataOUT(nameG03,countAtom,Multiplicity,
     & Stoichiometry, Ezpe, MolecularMass,ESN,Freq,DoF)
        ZeroEnergy=Ezpe
        rewind(4)
        return
       endif
      enddo 

19000  write(*,*) "Error reading gauss2multi.cfg file."
      close(4)
!     pause 
      stop "STOP"

20000  write(*,*) "Impossible to find wells. Check your 
     &gauss2multi.cgf file"
      close(4)
!     pause
      stop "STOP"
      end subroutine 
      
       

      subroutine singlefile(nameGausslog,OStext)
      character*72 nameGausslog,nameG03,namefile,find,TList
      character*72 Elimit
      character*5 Ostext
      character*4 Eunit, Punit
      double Precision ZeroEnergy, Ezpe
      INTEGER countAtom , DoF , ESN                       ! added explicit DECLARATION 
      CHARACTER*50  Stoichiometry                         ! added explicit DECLARATION
      PARAMETER (maxAtom=300)                             ! added explicit DECLARATION
      DOUBLE PRECISION MolecularMass , Freq(MaxAtom*3)    ! added explicit DECLARATION

      nameG03=nameGausslog
      call extractDataOUT(nameG03,countAtom,Multiplicity,Stoichiometry,
     & Ezpe, MolecularMass,ESN,Freq,Dof)
      ZeroEnergy=Ezpe
      
!      create Mominert data file and run it
      nameG03=nameGausslog
      namefile= nameG03(1:lnblnk(nameG03)-4)
      call writemominert(nameG03)
      call system ("mominert")
      find=OStext//"mominert.dat "//namefile(1:(lnblnk(namefile)))//
     &".coords"
      call system  (Find)
      find=OStext//"mominert.out "//namefile(1:(lnblnk(namefile)))//
     &".coords.out"
      call system (Find)
       
!---------------------------------------------------------------------------
!      create densum and run it
       Elimit="10    400     500     50000"
       call writedensum(nameG03,Elimit)
!     write(*,*) " Point 2"
       call system ("densum > NULL")
       find=OStext//" densum.dat "//namefile(1:(lnblnk(namefile)))//
     &".vibs"
       call system (Find)
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
!      create thermo and run it
       Eunit="KCAL"
       nTemp=1
       TList="298 "
       Punit="ATM"

       call writethermo(nameG03,ZeroEnergy,Eunit,nTemp,TList,Punit)
       call system ("thermo")
       find=OStext//" thermo.dat "//nameG03(1:lnblnk(nameG03)-4)//
     & ".therm"
        call system (find)
       find=OStext//" thermo.out "//nameG03(1:lnblnk(nameG03)-4)//
     & ".therm.out"
        call system (find)

      
      end subroutine
