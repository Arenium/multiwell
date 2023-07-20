       subroutine writeCFG()
       character*72 TempList,PressList, EnergyGrain
       character*3  Def
       character*4  Eunit, Punit
       character*72 nameG03, TypeInter
       integer i, nTemp, nPress

       open(unit=4,err=500,file='gauss2multi.cfg',status='new')
       
       write(*,*) "File gauss2multi.cfg not found "
       write(*,*) "It will be created now..."
10     write(*,*)
       write(*,*) "Do you want to use default value ? (Y/n)"
       write(*,*) "(Kcal/mol, 1 BAR, 298 K, Egrain=10, Emax2=50000" 
       read(*,*) Def
       Def=Def(1:lnblnk(Def))
       write(*,*) 

               Eunit="KCAL"
               Punit="BAR"
               nTemp=1
               TempList="298 "
               nPress=1
               PressList="1  "
               EnergyGrain="10    400     900     50000  "

       if (Def.eq."y".or.Def.eq."Y".or.Def.eq."Yes".or.Def.eq."YES".or.
     &Def.eq."yes") goto 200
       
       write(*,*) " Energy unit : (KCAL / KJOU / CM-1)"
       read(*,*) Eunit
       write(*,*) " Number of temperatures: "
       read(*,*) nTemp
       write(*,*) " List of ",nTemp," temperatures (use space to separat
     &e) :"
       read(*,'(A72)') TempList
       write(*,*) " Pressure units : (TOR / ATM / BAR / MCC) "
       read(*,*) Punit
       IF ( Punit. EQ. 'TOR' ) THEN
         write(*,*)'   (Note: BAR will be used instead of TOR in the'
         write(*,*)'          thermo.dat data file.)'
       ENDIF
       write(*,*) " Number of pressures: "
       read(*,*) nPress
       write(*,*) " List of ",nPress, " pressures (use space to separate
     &) :"
       read(*,'(A72)') PressList
       write(*,*) " Energy Grain, Imax1, Isize, Emax2 (separated by spac
     &es):"
       write(*,*) " (Recommended:  10  400  900  50000)"
       read(*,'(A72)') EnergyGrain

200    write(4,*) Eunit
       write(4,*) nTemp
       write(4,*) TempList
       write(4,*) Punit
       write(4,*) nPress
       write(4,*) PressList
       write(4,*) EnergyGrain 
      
       i=1
300    write(*,*)"full name of Gaussian output file (type 0 to finish):"
       read(*,*) nameG03 
       if (nameG03(1:lnblnk(nameG03)).ne."0") then 
        write(*,*) "Type of structure : (PROD, TS, WELL)"
        read(*,*) TypeInter
        write(4,*) i ,"  ",nameG03(1:lnblnk(nameG03)), "   ", 
     &  TypeInter(1:lnblnk(TypeInter))
        i=i+1
        goto 300
       endif
       close(4)
       return

500   write(*,*) " gauss2multi.cfg file already exists."
      write(*,*) " Do you want to overwrite ? (Y,N)"
      read(*,*) Def
      if (Def.eq."n".or.Def.eq."N".or.Def.eq."No".or.Def.eq."NO".or.
     & Def.eq."no") then
       return
      else
       open(unit=4,file='gauss2multi.cfg',status='old')
       goto 10
      endif
 
      end subroutine

