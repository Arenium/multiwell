c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c copyright (c) 2014 john r. barker, jason a. sonk
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

      subroutine get_viblo(ls,lf,viblo)
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      integer ls,lf
      real*8  viblo

      viblo=100000.0
      do i=ls,lf
         do j=1,ndof(i)
            if(idofl(i,j).eq.'vib')then
               viblo=min(viblo,wel(i,j))                                
            elseif(idofl(i,j).eq.'hra')then
               viblo=min(viblo,wel(i,j))                                
            elseif(idofl(i,j).eq.'hrb')then
               viblo=min(viblo,wel(i,j))
            elseif(idofl(i,j).eq.'hrc')then
               viblo=min(viblo,wel(i,j))
            else
               continue
            end if
         end do  
      end do        

      return
      end subroutine     
