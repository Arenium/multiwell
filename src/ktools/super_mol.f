c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c copyright (c) 2017 john r. barker, jason a. sonk
c
c john r. barker
c jrbarker@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c (734) 763 6239
c
c this program is free software; you can redistribute it and/or
c modify it under the terms of the gnu general public license (version 2)
c as published by the free software foundation.
c
c this program is distributed in the hope that it will be useful,
c but without any warranty; without even the implied warranty of
c merchantability or fitness for a particular purpose. see the
c gnu general public license for more details.
c
c see the 'readme' file for a copy of the gnu general public license,
c or contact:
c
c free software foundation, inc.
c 59 temple place - suite 330
c boston, ma 02111-1307, usa.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine super_mol
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      real*8     mw(rcnt+pcnt),amass,hfxn,sfxn
      character  for*20
      dimension  for(rcnt+pcnt)
 100  format (i3,2x,A20,2x,2(F10.4,2x))
 200  format (A12,1x,'mass',2x,5(F10.4,2x))
 300  format (19x,2(A10,2x)) 

      call sestamp('super_mol',1)

c  calculate molecular weights of all r,ts & p

      write(lunit,*)'molecular weights of r, ts & p'
      write(lunit,*)

      do i=1,rcnt+ntts+pcnt

         call chemformula(i,trim(formulal(i)))
         amu(i)=0.0d0

         do j=1,nelement(i)
            call element(atype(i,j),amass,hfxn,sfxn)
            amu(i)=amu(i)+natom(i,j)*amass
         end do

         write(lunit,100)i,trim(formulal(i)),amu(i)

      end do
      write(lunit,*)

c  find total & reduced  masses for each r & p


      if(rcnt.ge.1)then
         write(lunit,*)'super reactant'
         do i=1,rcnt
            write(lunit,100)i,trim(formulal(i)),amu(i)
            tamu(1)=tamu(1)+amu(i)
            ramu(1)=ramu(1)+(1.0d0/amu(i))
         end do
         ramu(1)=(1.0d0/ramu(1))
         write(lunit,*)
      end if

      if(pcnt.ge.1)then
         write(lunit,*)'super product'
         do i=rcnt+ntts+1,rcnt+ntts+pcnt
            write(lunit,100)i-rcnt-ntts,trim(formulal(i)),amu(i)
            tamu(2)=tamu(2)+amu(i)
            ramu(2)=ramu(2)+(1.0d0/amu(i))
         end do
         ramu(2)=(1.0d0/ramu(2))
         write(lunit,*)
      end if

      write(lunit,300)'reactant','product'
      write(lunit,200)'total',(tamu(i),i=1,2)
      write(lunit,200)'reduced',(ramu(i),i=1,2)
      write(lunit,*)
      write(lunit,*)

      call sestamp('super_mol',2)

      return

      end subroutine

