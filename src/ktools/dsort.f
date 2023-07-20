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
      subroutine  dsort(sortcnt)
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      integer sortarr(rcnt+ntts+pcnt)
      integer rankarr(rcnt+ntts+pcnt)
      integer  sortcnt
      logical sorted
 100  format(f3.1)
 200  format(i2,2x,i2,2x,f3.1,2x,f3.1,2x,l2)
 300  format(i2,2x,f5.1,2x,i2)

c  initallize sortarr
      do i=1,rcnt+ntts+pcnt
         if(reprod(i).eq.'prod')distl(i)=distl(i)+1000.0d0
         if(reprod(i).eq.'reac')distl(i)=distl(i)-1000.0d0
         sortarr(i)=i
      end do
      sorted=.true.

c  find index array which sorts data

      call indexx(rcnt+ntts+pcnt,distl,sortarr)
      call rank(rcnt+ntts+pcnt,sortarr,rankarr)

      do i=1,rcnt+ntts  
            k=sortarr(i)
            if(k.ne.i)sorted=.false.
      end do

      if(.not.sorted)call sort_input(rankarr,sortcnt)

      if(sortcnt.gt.10)then
         write(*,*)'Error in sorting. Check distances in input file'
         stop
      end if

      do i=1,rcnt+ntts+pcnt
         if(reprod(i).eq.'prod')distl(i)=distl(i)-1000.0d0
         if(reprod(i).eq.'reac')distl(i)=distl(i)+1000.0d0
      end do

      return

      end subroutine
