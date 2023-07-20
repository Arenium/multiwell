c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Copyright (C) 2009 Andrea Maranzana and John R. Barker
c
c Andrea Maranzana
c andrea.maranzana@unito.it
c Department of General and Organic Chemistry
c University of Torino
c Corso Massimo D'Azeglio, 48
c Torino  10125
c ITALY
c ++39-011-670-7637
c
c John R. Barker
c jrbarker@umich.edu
c University of Michigan
c Ann Arbor, MI 48109-2143
c (734) 763 6239
c
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License (version 2)
c as published by the Free Software Foundation.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
c GNU General Public License for more details.
c
c See the 'ReadMe' file for a copy of the GNU General Public License,
c or contact:
c
c Free Software Foundation, Inc.
c 59 Temple Place - Suite 330
c Boston, MA 02111-1307, USA.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      Subroutine chemformula(i,formula)

      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      character(len=*) formula
      character(len=90) char2
      character(len=4)  char1, char3
      integer(4) y, z, LL, maxlen, MaxL
      parameter (MaxL=101)  
      integer(4) y1, y2, digit, ntype , niso
      integer(4) numtot, times
      integer(4) Natomtemp(MaxL)

      parameter(ntype=101)
      parameter(niso=18)           ! number of isotopes (e.g. 'Si30')
      character(len=4) ATYPETEMP(maxL)
      character(len=4) AtomType(ntype)

      data AtomType/' C12',' C13',' C14', ' C16', ' O16', ' O17',' O18',
     2 ' H1','N14','N15','Cl35','Cl37','Br79','Br81','S32','S34','Si29',
     3 'Si30', 
     3 ' H', ' D', ' T', ' He',  ' Li', ' Be',  ' B', ' C',  ' N', ' O',
     4 ' F','Ne','Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc',
     5 'Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga',
     6 'Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru',
     7 'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe','Cs','Ba','La',
     8 'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     9 'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg','Pb'/

      LL=0


!------------(R)n  ->  RRRRRR  ---------------------------------------
      char1="("
      char3=")"
      do
       y1=index(TRIM(formula),TRIM(char1))
       y2=index(TRIM(formula),TRIM(char3))
       IF(y1.eq.0.and.y2.eq.0) exit
       IF(y1.eq.0.and.y2.ne.0) STOP 'Error in chemical formula'
       IF(y1.ne.0.and.y2.eq.0) STOP 'Error in chemical formula'
       
       if(y2.eq.len(TRIM(formula))) then
        call delete(formula,y2,y2)
        call delete(formula,y1,y1)
        exit
       endif

       call checknumber(formula,y2+1,times)
       if(times.gt.1) then 
        call repeatText(formula,y1,y2,times)
       else
        call delete(formula,y2,y2)
        call delete(formula,y1,y1)
       endif

      enddo
!--------------------------------------------------------------------
      

100      k=0
      do 
         k=k+1
         if(k.gt.ntype) THEN
          OPEN (UNIT=10,STATUS='UNKNOWN',FILE='thermo.ERROR')          ! Error output
          write (10,*) 'Element not recognized in chemformula.f'
          write (*,*) 'Element not recognized in chemformula.f'
          CLOSE(10)
          STOP
         ENDIF

         char1=ADJUSTL(AtomType(k))
         char2=TRIM(ADJUSTL(formula))
         IF(char2.EQ.'') exit
         z=index(TRIM(formula),TRIM(char1))
         maxlen=len(TRIM(char1))
         if(z.eq.0) GOTO 1000
         IF(k.le.niso.AND.z.eq.1) GOTO 1000
         
         IF(z.gt.1) then
          IF(char2(z-1:z-1).ne."[".AND.k.le.niso) goto 1000
         endif

!   ---------to avoid C for Cr, Cu etc..-------------------------
         y1=ICHAR(char2(z+1:z+1))
         IF(y1.ge.97.and.maxlen.eq.1) goto 1000
!   -------------------------------------------------------------

         IF(k.le.niso) then    ! isotopes [XX]
                   y=1
         else
                   y=0
         endif

          LL=LL+1

          IF(LL.gt.MaxL) then 
           OPEN (UNIT=10,STATUS='UNKNOWN',FILE='thermo.ERROR')   !  Error output
           write (10,*) 'Too many elements in chemformula.f'
           write (*,*)  'Too many elements in in chemformula.f'
           CLOSE(10)
           STOP
          ENDIF


          ATYPETEMP(LL)=TRIM(char1)
  
          call checknumber(formula,z+maxlen+y,times)
          IF(times.le.9) then
           digit=1
          else
           digit=2
          endif
          if(times.le.0) then
             digit=0
             times=1
          endif
          Natomtemp(LL)=times
          call delete(formula,z-y,z+maxlen+digit+y-1)
          GOTO 100
         
1000  CONTINUE
      enddo

!---------count ATYPE and Nelement-----------------------------------

      !ATYPE(i,ne) Natom(i,ne) Nelement(i)
      
      Nelement(i)=0
      do y1= 1, ntype
       numtot=0
       do y2= 1, LL
        IF(TRIM(ADJUSTL(ATYPETEMP(y2))).eq.TRIM(ADJUSTL(AtomType(y1))))
     &    numtot=numtot+Natomtemp(y2)
       enddo
       IF(numtot.ne.0) then 
        Nelement(i)=Nelement(i)+1
        ATYPE(i,Nelement(i))=TRIM(ADJUSTL(AtomType(y1)))
        Natom(i,Nelement(i))=numtot
       endif
      enddo


      END subroutine


!--------------Delete characters from start to finish ---------------
      subroutine delete(formula,start,finish)
      CHARACTER(len=*) formula
      integer start,finish

       formula=(ADJUSTL(formula(1:start-1))
     &                    //ADJUSTL(formula(finish+1:len(formula))))
       
      return

      end subroutine


!--------------Repeat characters from start to finish  "times" times --
      subroutine repeatText(formula,start,finish,times)
      CHARACTER(len=*) formula
      CHARACTER(len=90) copytxt, txt
      integer start,finish,times, xx, digit
      IF(times.le.9) then
       digit=1
      else
       digit=2
      endif

      txt=""
      copytxt=TRIM(ADJUSTL(formula(start+1:finish-1)))
      do xx=1, times
       txt=TRIM(ADJUSTL(txt))//TRIM(ADJUSTL(copytxt)) 
      enddo
       formula=(ADJUSTL(formula(1:start-1))
     &  //TRIM(ADJUSTL(txt))//formula(finish+1+digit:len(formula)))

      return

      end subroutine


!--------------Check if character in start is a number ----------------
      subroutine Checknumber(formula,start,number)
      CHARACTER(len=*) formula
      CHARACTER(len=90) txt
      integer start,number, xx, numtemp(2)

      if(start.gt.LEN(TRIM(ADJUSTL(formula)))) then
        number=-1
        return
      endif

      do xx=0, 1
      if(start+xx.gt.LEN(TRIM(ADJUSTL(formula)))) then
        numtemp(2)=-1
        exit
      endif
       txt=TRIM(ADJUSTL(formula(start+xx:start+xx)))
       IF(ICHAR(TRIM(txt)).ge.48.AND.ICHAR(TRIM(txt)).le.57) then
        Read( txt,*)  numtemp(xx+1)
       ELSE
        numtemp(xx+1)=-1
       ENDIF
      enddo

       IF(numtemp(1).eq.0) GOTO 300
       IF(numtemp(1).lt.0) then
        number=-1
        return
       ENDIF
       IF(numtemp(2).lt.0) then
        number=numtemp(1) 
        return
       ELSE
        txt=TRIM(ADJUSTL(formula(start:start+1)))
        Read( txt,*) number
        return
       ENDIF

300    OPEN (UNIT=10,STATUS='UNKNOWN',FILE='thermo.ERROR')    ! Error output
       write (10,*) 'Error in the chemical formula (chemformula.f)'
       write (*,*)  'Error in the chemical formula (chemformula.f)'
       CLOSE(10)
       STOP

      return

      end subroutine


