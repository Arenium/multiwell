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
      subroutine chkdens(dens,n,egrain1,imax1,isize,emax2)
      
      integer bestimax(6), imax1, isize, x,n
      double precision emax1, emax2, egrain1, dens(n)

      emax1 = (imax1-1)*egrain1
      x=1
      do i=1,6
       bestimax(i)=0
      enddo

      do i=(imax1-1), 2, -1 
       if(dens(i-1).lt.1e-2) goto 100
       err=abs((dens(i)-dens(i-1))/dens(i-1))*100
c      write(*,*) i,err
       if(i.eq.(imax1-1)) then
         if(err.gt.1.and.err.lt.2) then
           write(*,*) " NOTE: density fluctuation at emax1 1-2 %."
           bestimax(1)=(imax1-1)
           x=2
         elseif(err.gt.2.and.err.lt.3) then
           write(*,*) " NOTE: density fluctuation at emax1 2-3 %."
           bestimax(1)=(imax1-1)
           bestimax(2)=(imax1-1)
           x=3
         elseif(err.gt.3.and.err.lt.4) then
           write(*,*) " WARNING: density fluctuation at emax1 3-4 %."
           bestimax(1)=(imax1-1)
           bestimax(2)=(imax1-1)
           bestimax(3)=(imax1-1)
           x=4
         elseif(err.gt.4.and.err.lt.5) then
           write(*,*) " WARNING: density fluctuation at emax1 4-5 %."
           bestimax(1)=(imax1-1)
           bestimax(2)=(imax1-1)
           bestimax(3)=(imax1-1)
           bestimax(4)=(imax1-1)
           x=5
         elseif(err.gt.5) then
           write(*,*) " TOO LARGE!! DENSITY FLUCTUATION AT emax1 > 5 %."
           bestimax(1)=(imax1-1)
           bestimax(2)=(imax1-1)
           bestimax(3)=(imax1-1)
           bestimax(4)=(imax1-1)
           bestimax(5)=(imax1-1)
           return
         endif

       endif

       if(err.gt.x) then              ! error > 1%-5%
        bestimax(x)=i
!     bestimax(1): error< 1%, bestimax(2): error< 2%, bestimax(3): error< 3%
!     bestimax(4): error< 4%, bestimax(5): error< 5%
        x=x+1
        if(x.eq.6) goto 100
       endif
      

      enddo
c      write(*,*) "wrong definition of energy parameters"
      write(*,*)"cannot evaluate fluctuations using these e parameters"
c      stop

100   if(imax1.lt.bestimax(5)) then
       write(*,'(a45,i6)')" warning: miminum suggested value for imax1:"
     & , bestimax(5),imax1
      endif
      return
 

      end subroutine



