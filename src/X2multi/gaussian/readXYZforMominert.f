      subroutine readXYZforMominert(countAtom,NumAtom,X,Y,Z)

      PARAMETER (maxAtom=900)
      integer nul,i, NumAtom(maxAtom),countAtom
      double precision X(maxAtom),Y(maxAtom),Z(maxAtom)
      CHARACTER Find*46

      countAtom=0

200   FORMAT(A46)     

      do 
      read(1,200) Find
      if (Find.eq."                         Standard orientation:") exit   ! Gaussian03
      if (Find.eq."                         Z-Matrix orientation:") exit   ! Gaussian98
      if (Find.eq."                          Input orientation:  ") exit   ! Gaussian03 nosymm
      end do

! Read Atom number, X,Y,Z
      read(1,200) 
      read(1,200) 
      read(1,200) 
      read(1,200) 
      do
      countAtom=countAtom+1
      read(1,*) nul,NumAtom(countAtom),nul,X(countAtom),Y(countAtom),
     &  Z(countAtom)
      if (NumAtom(countAtom).lt.1) countAtom=countAtom-1
      read(1,200) Find
      if (Find.eq." ---------------------------------------------") then
      exit 
      else
      backspace(1)
      endif
      enddo


      end subroutine


!------------------------------------------------------
! 
      subroutine CheckFreqLog(nameG03)
      character*72 NameG03
      character*24 FindFreq
      do
      read(1,'(A24)',end=100) FindFreq
      if (FindFreq.eq." Optimization completed.") 
     &  call ReadOUT_WriteMominert(nameG03)
      if (FindFreq(1:15).eq." Zero-point cor") then
      return
      endif
      enddo 

100   Write(*,*)  "  WARNING "
      Write(*,*) "  It does not seem a frequencies calculation output fi
     &le. "
      stop
   
      end subroutine
!"----------------------------------------------------------------
