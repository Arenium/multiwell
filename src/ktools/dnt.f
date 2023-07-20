c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c copyright (c) 2014 jason a. sonk
c
c jason a. sonk
c jsonk@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c
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

c  !==========================================================!
c  ! retrieve date and time                                   !
c  !==========================================================!
      subroutine dnt(iout)
      character date*8, time*10, zone*5, month*3, day*3, ampm*2,st*2
      integer values(8),iout

 100  format(a,x,i2,x,i4)
 200  format(i2,':',i2,':',i2,':',i3,x,a)
 300  format(a,'.',x,i2,a,x,i4,3x,i2,':',i2.2,':',i2.2,':',i3.3,x,a)
 400  format(a,x,i2,x,i4,3x,i2,':',i2.2,':',i2.2,':',i3.3,x,a)

c write to iout formatted date and time

      call date_and_time(date,time,zone,values)
      select case (values(2))
         case (1)
                 month = 'jan'
         case (2)
                 month = 'feb'
         case (3)
                 month = 'mar'
         case (4)
                 month = 'apr'
         case (5)
                 month = 'may'
         case (6)
                 month = 'jun'
         case (7)
                 month = 'jul'
         case (8)
                 month = 'aug'
         case (9)
                 month = 'sep'
         case (10)
                 month = 'oct'
         case (11)
                 month = 'nov'
         case (12)
                 month = 'dec'
      end select

      if (values(5) >= 12) then
              ampm = "pm"
              values(5) = mod(values(5),12)
              if(values(5).eq.0)values(5)=12
      else if (values(5) == 0) then
              ampm = "am"
              values(5) = 12
      else
              ampm = "am"
      end if

      select case (values(3))
         case (1)
                  st = 'st'
         case (21) 
                  st = 'st'
         case (31) 
                  st = 'st'
         case (2) 
                  st = 'nd'
         case (22) 
                  st = 'nd'
         case (3) 
                  st = 'rd'
         case (23) 
                  st = 'rd'
         case default
                  st = 'th'
      end select


      write(iout,300)month,values(3),st,values(1),values(5),values(6),
     $ values(7),values(8),ampm

      end subroutine dnt

