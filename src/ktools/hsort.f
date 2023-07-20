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
      subroutine hsort
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      real*8 temp(rcnt+ntts+pcnt),minh(2)
      real*8 tomp(rcnt+ntts+pcnt)
      real time   
      character uout*4
      parameter (uout = 'cm-1')

 111  format(5x,22(es20.8))
 222  format(a4,1x,1(i3,2x),4(F12.2,2x))
 333  format(a4,1x,1(a3,2x),4(a12,2x))

      call sestamp('hsort',1)

c sort potential energy data into arrays 
c (sorting was incorporated into read_inp)
c  for array delh of length rcnt+ntts+pcnt
c     array elements 1 will be the reactants
c     array elements 2-rcnt+1 through ntts will be trial transisition states
c     array elements rcnt+2 will be the products


      if((eunits.eq.'kcal').or.(eunits.eq.'kj'))then
         write(lunit,*)'converting pe from ',eunits,'/mol to cm-1'
      else
         write(lunit,*)'converting pe from ',eunits,' to cm-1'
      end if
      write(lunit,*)

      do i=1,rcnt+ntts+pcnt
         temp(i)=0.0d0
         tomp(i)=0.0d0
      end do

      do i=1,rcnt+ntts+pcnt
         if (reprod(i).eq.'reac')then
            temp(1)=temp(1)+delh(i)
         elseif(reprod(i).eq.'ctst')then
            temp(i)=delh(i)
         elseif(reprod(i).eq.'prod')then
            temp(rcnt+ntts+1)=temp(rcnt+ntts+1)+delh(i)
c         elseif(rcnt.le.1)then
c            temp(i)=delh(i)
c         elseif(rcnt.ge.2)then
c            temp(i-1)=delh(i)
         end if
      end do

      do i=1,ntts+rcnt+pcnt
         tomp(i)=temp(i)
      end do

c
c call: energy conversion
c   array delh(rcnt+ntts) contains user submitted energies for species
c   array eunits(rcnt+ntts) contains the units of delh for each species

      call nrgconvert(tomp,rcnt+ntts+pcnt,eunits,uout)

c log print

      write(lunit,*)
      if((eunits.eq.'kcal').or.(eunits.eq.'kj'))then
         if(pcnt.ne.0)then
            write(lunit,333)'type','#',eunits//'/mol','cm-1',
     $'forward cm-1','reverse cm-1'
         else
            write(lunit,333)'type','#',eunits//'/mol','cm-1',
     $'forward cm-1'
         end if
      else
         if(pcnt.ne.0)then
            write(lunit,333)'type','#',eunits,'cm-1','forward cm-1',
     $'reverse cm-1'
         else
            write(lunit,333)'type','#',eunits,'cm-1','forward cm-1'
         end if
      end if
      write(lunit,*)('-',i=1,66)

c  setting forward and reverse zeros of enthalpy is taken as point 1

      minh(1)=tomp(1)
      if(pcnt.ne.0)then
         minh(2)=tomp(rcnt+ntts+1)
      end if

      if(pcnt.ne.0)then
         do i=1, rcnt+ntts+pcnt
            delh(i)=tomp(i)-minh(1)
           delhf(i)=tomp(i)-minh(1)
           delhr(i)=tomp(i)-minh(2)
         end do

         temp(rcnt+ntts+pcnt)=temp(rcnt+ntts+1)
         tomp(rcnt+ntts+pcnt)=tomp(rcnt+ntts+1)
         delh(rcnt+ntts+pcnt)=delh(rcnt+ntts+1)
         delhf(rcnt+ntts+pcnt)=delhf(rcnt+ntts+1)
         delhr(rcnt+ntts+pcnt)=delhr(rcnt+ntts+1)

         do i=1,1
           write(lunit,222)reprod(i),i,temp(i),tomp(i),delhf(i),delhr(i)
         end do
         do i=rcnt+1,rcnt+ntts+pcnt
           write(lunit,222)reprod(i),i,temp(i),tomp(i),delhf(i),delhr(i)
         end do
      else
         do i=1, rcnt+ntts
            delh(i)=tomp(i)-minh(1)
           delhf(i)=tomp(i)-minh(1)
           write(lunit,222)reprod(i),i,temp(i),tomp(i),delhf(i)
         end do
      end if

C if there are "products" and their reaction energy is greater than that of the 
c final vts then use their energy as the vts energy

      if(pcnt.gt.0)then
         if(delhf(rcnt+ntts+1).gt.delhf(rcnt+ntts))then
            delh(rcnt+ntts)=delh(rcnt+ntts+1)
            delhf(rcnt+ntts)=delhf(rcnt+ntts+1)
            delhr(rcnt+ntts)=delhr(rcnt+ntts+1)
         end if
      end if


      call sestamp('hsort',2)

      return

      end subroutine hsort

