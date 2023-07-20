      subroutine writemultiwelldat 
      integer i,ntot,numvib,x, ESN(300), Multiplicity(300)
      character*72 Find,nameFile
      character*72 PressList,Egrain, nameG03(300),nameThermo
      character*4  Eunit, Punit ,TypeInt(300)    
      integer nPress , nWell ,nProd, nTS ,num(300)
      double precision ADrot(300)
      double precision Energy(300)

      open(unit=1,file='gauss2multi.cfg',status='old')
      open(unit=3,file='multiwell.dat',status='unknown')
      
      read(1,*) Eunit
      read(1,*)
      read(1,*)  
      read(1,*) Punit
      read(1,*) nPress
      read(1,'(A72)') PressList
      read(1,'(A72)') Egrain
      write(3,*) " Multiwell.dat auto-generated.  "
      write(3,*) Egrain(1:lnblnk(Egrain)), "   ", irand(time())
      write(3,*) "'",Punit,"'  '",Eunit,"'  'AMUA'"
      write(3,*) "298   298      !   <-  translational and initial vibra
     &tional temperatures."
      write(3,'(I3)') nPress
      write(3,*) PressList
      
      nWell=0
      nProd=0
      nTS=0 
      ntot=0
      x=0

      do
      ntot=ntot+1
      read(1,*,end=10000,err=10000)num(ntot),nameG03(ntot),TypeInt(ntot)
      nameG03(ntot)=nameG03(ntot)(1:lnblnk(nameG03(ntot))-4)
      if(TypeInt(ntot).eq."WELL") nWell=nWell+1
      if(TypeInt(ntot).eq."PROD") nProd=nProd+1
      if(TypeInt(ntot).eq."TS") nTS=nTS+1
      if(TypeInt(ntot).ne."TS".AND.TypeInt(ntot).ne."PROD".AND.TypeInt(n
     &tot).ne."WELL") then
      write(*,*) "Error in ", nameG03(ntot)(1:lnblnk(nameG03(ntot)))
      write(*,*) "Type of structure not recognized."
      write(*,*) "Check your gauss2multi.cfg file."
      close(1)
      close(3)
!     pause
      stop  "STOP"
      endif
      enddo
10000 close(1)
      ntot=ntot-1
      write(3,'(I3,A2,I3)') nWell,"  ", nProd
      

      do i=1,ntot
      nameThermo=nameG03(i)(1:lnblnk(nameG03(i)))//".therm"
      open(unit=2,file=nameThermo,status='old')
      do n=1,4
       read(2,*) 
      enddo

      read(2,*) Find,Find,Energy(i)
      do n=1,4
       read(2,*)
      enddo

      read(2,*) ESN(i)
      read(2,*) Find,Multiplicity(i)
      read(2,*) numvib
      do n=1,numvib-1
      read(2,*)
      enddo
      read(2,*) Find,Find,ADrot(i)
      close(2)
      enddo

!     write wells
      do i=1,ntot
      if(TypeInt(i).eq."WELL") then
       x=x+1
       nameFile="  '"//nameG03(i)(1:lnblnk(nameG03(i)))//"'  "
       write(3,100) x, nameFile(1:lnblnk(nameFile)),Energy(i),
     & ADrot(i), ESN(i), Multiplicity(i), "  1  "
      endif
      enddo
100   FORMAT(T1,I2,T3,A10,T23,F12.4,F10.4,T47,I2,I3,A5)

      
!     write products and reactants.
      do i=1,ntot
      if(TypeInt(i).eq."PROD") then
      x=x+1
      nameFile="  '"//nameG03(i)(1:lnblnk(nameG03(i)))//"' "
      write(3,200) x, nameFile(1:lnblnk(nameFile)),Energy(i)
      endif
      enddo
200   FORMAT(T1,I2,T3,A10,T23,1x,F12.4)


      write(3,*) "3.74   82.    28.0     42.  
     &          !  <- Check this line    ! N2 Collider "

      do i=1,nWell
      write(3,'(I3,A87)') i, "4.8   270.   1   175.   0.0  0.0  0.0  0.0
     &  0.0  0.0  0.0   !  <- Check these numbers"  
      write(3,*) "LJ" 
      enddo

!     write TS
      write(3,'(I3,A87)') nTS,  
     & "                                              !  <- Check this 
     &number "

      do i=1,ntot
      if(TypeInt(i).eq."TS") then
      nameFile="  '"//nameG03(i)(1:lnblnk(nameG03(i)))//"' "
      write(3,300) "X  Y",nameFile(1:lnblnk(nameFile)),ADrot(i),
     & ESN(i),Multiplicity(i),
     &"  1  1.0e+16  Del-E 'rev' 'NOTUN' 'FAST' 'cent2' 'sum'  ! <- Inse
     &rt indexes and del-E"
      endif
      enddo

      write(3,*) "1000  'COLL' 1000     'CHEMACT'   X   Y   0.              
     &          !  <- Check this line "
      write(3,*)
      write(3,*)
      write(3,*) " This is only a template, you have to check and insert
     & the right numbers !!"
      write(3,*)
      close(3)
 
300   format(A,A," ",F8.3,I3,I3,A)
      end subroutine
    
