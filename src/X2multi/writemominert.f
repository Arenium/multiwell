      subroutine writemominert(nameG03)
      character*72 nameG03
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
! read aces output file (nameG03): find cartesian coordinates, read X,Y,Z and
! write  momintert.dat
! LIMIT: 900 atoms
!
! 
!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      open(unit=1,file=nameG03(1:lnblnk(nameG03)),status='old')
!  controllo se si tratta di un calcolo di freq.
      call ReadOUT_WriteMominert(nameG03)
      call CheckFreqlog(nameG03)

      close(1) 
      end subroutine
c----------------------------------------------------------------------------------      
c----------------------------------------------------------------------------------      
! Find cartesian coordinates
      subroutine ReadOUT_WriteMominert(nameG03)

      character*2 AtomType(100)
      data AtomType/ ' H','He', 'Li','Be',' B',' C',' N',
     2 ' O',' F','Ne','Na','Mg','Al','Si',' P',' S','Cl','Ar',
     3 ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga',
     4 'Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru',
     5 'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe','Cs','Ba','La',
     6 'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     7 'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi',
     8 'Po','At','Rn','Fr','Ra','Ac','Th','Pa',' U','Np','Pu','Am','Cm',
     9 'Bk','Cf', 'XX',' '/

      PARAMETER (maxAtom=900)
      integer nul,i, NumAtom(maxAtom),countAtom
      double precision X(maxAtom),Y(maxAtom),Z(maxAtom)
      CHARACTER Find*46
      CHARACTER nameG03*72
      call readXYZforMominert(countAtom,NumAtom,X,Y,Z)
!     countAtom=0
!
!200   FORMAT(A46)     
!
!      do 
!      read(1,200) Find
!      if (Find.eq."  Symbol    Number           X              Y ") exit
! !     end do
!
! Read Atom number, X,Y,Z
!!      read(1,*) 
!      do
!      countAtom=countAtom+1
! !     read(1,*) Find,NumAtom(countAtom),X(countAtom),Y(countAtom),
!     &  Z(countAtom)
!      if (NumAtom(countAtom).lt.1) countAtom=countAtom-1
!      read(1,200) Find
!      if (Find.eq." ---------------------------------------------") then
!      exit 
!      else
!      backspace(1)
!      endif
!      enddo
!
!write mominert.dat 

      open(unit=2,file='mominert.dat',status='unknown')
      write(2,*) nameG03(1:lnblnk(nameG03)-4)
      write(2,*) 'ANGS'
      write(2,'(I3)') countAtom
      do i=1,countAtom
        write(2,'(A1,A2,A2,I3,3F12.6)') " ",AtomType(NumAtom(i)),"  ",i,
     &  X(i),Y(i),Z(i)
      enddo
      write(2,*) "0 , 0"
      write(2,*) 
      write(2,*) 
      write(2,*)  "Please check internal rotors"
      write(2,*) 

      close(2)
      
      end subroutine


!------------------------------------------------------
! 
!     subroutine CheckFreqLog(nameG03)
!     character*72 NameG03
!     character*24 FindFreq
!     do
!      read(1,'(A24)',end=100) FindFreq
!      if (FindFreq.eq."  Zero-point vibrational") return
!     enddo 
!
!100   Write(*,*)  "  WARNING "
!     Write(*,*) "  It does not seem a frequencies calculation output fi
!    &le. "
!     stop
!      Write(*,*) " Do you know what you are doing ?"
   
!     end subroutine
!"----------------------------------------------------------------
