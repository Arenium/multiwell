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
      subroutine bsort
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      real*8 temparr1(datpts),temparr2(datpts)
      logical kfound
      character uout*4
      parameter (uout = 'cm-1')

 111  format(5x,22(es20.8))
 222  format(1x, 2(i4, 2x), 3(f10.3, 2x) )
 333  format(1x,i3,2x,4(f9.4,2x))
 444  format(1x,2(a4,2x),3(a9,2x))
 555  format(1x,a3,2x,4(a9,2x))

      call sestamp('bsort',1)

c  look for 2d external j-rotor dof labeled 'jro'

      write(lunit,*)'converting j-rotational constant into cm-1'
      write(lunit,*)
      do i=1, rcnt+ntts+pcnt
         temparr1(i)=0.0d0
         barr(i)=0.0d0
         do j=1,ndof(i)
            if(idofl(i,j).eq.'jro')then
               temparr1(i)=wel(i,j)
               call nrgconvert(wel(i,j),1,keyword(i,2),uout)
               barr(i)=wel(i,j)
               if(i.eq.1)then
                  write(lunit,444)'#','dof#',keyword(i,2),'cm-1','cm-1'
                  write(lunit,*)('-',k=1,47)
               end if
               write(lunit,222)i,j,temparr1(i),wel(i,j),barr(i)
            end if
          end do
      end do
      write(lunit,*)

c  look for 1d external k-rotor dof labeled 'kro' linear molecules will not have one

      write(lunit,*)'converting k-rotational constant into cm-1'
      write(lunit,*)
      do i=1, rcnt+ntts+pcnt
         kfound=.false.
         do j=1,ndof(i)
            if(idofl(i,j).eq.'kro')then
               kfound=.true.
               temparr2(i)=wel(i,j)
               call nrgconvert(wel(i,j),1,keyword(i,2),uout)
               aarr(i)=wel(i,j)
               if(i.eq.1)then
                  write(lunit,444)'#','dof#',keyword(i,2),'cm-1','cm-1'
                  write(lunit,*)('-',k=1,47)
               end if
               write(lunit,222)i,j,temparr2(i),wel(i,j),aarr(i)
            end if
         end do

         if(.not.kfound)then
            write(lunit,222)i,j-1,temparr2(i),wel(i,j),aarr(i)
         end if

      end do
      write(lunit,*)

c  calculate effective (a-b) rotational constant

      write(lunit,*)'(a-b) rotational constant'
      write(lunit,*)
      do i=1, rcnt+ntts+pcnt
         do j=1, ndof(i)
            if(idofl(i,j).eq.'kro')then
c              wel(i,j)=aarr(i)-barr(i)
             if(i.eq.1)then
                write(lunit,555)'#','a','b','a-b','check'
                write(lunit,*)('-',k=1,47)
             end if
cc             write(lunit,333)i,aarr(i),barr(i),aarr(i)-barr(i),aarr(i)
             write(lunit,333)i,aarr(i),barr(i),aarr(i)-barr(i),wel(i,j)
c             aarr(i)=wel(i,j)
            end if
         end do
      end do

      call sestamp('bsort',2)

      return

      end subroutine

