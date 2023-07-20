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
      subroutine multiwell_write
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      integer      tunit
      character*80 dir,command,sdir
      character*80 fname
      dimension    sdir(rcnt+ntts+pcnt)
      logical      there
      parameter    (tunit=24)      
 100  format(i3,1x,a3,2(2x,f9.4),2x,i5)
 200  format("     ",a,2x,i4,2x,f6.2,10(2x,f10.2))
 300  format(f4.0,2x,i4,2x,i5,2x,i6,2x)
 400  format(i3,2x,i2,2x,a4,2x,a5,2x)

      dir="densum"
      sdir(1)="react"
      do i=rcnt+2,ntts
         sdir(i)="ts"
      end do
      sdir(3)="product"

      inquire(file=dir,exist=there)
      if(.not.there)then
            command='mkdir '//trim(dir)
            write(*,*)'making ',trim(dir)
            call system(command)
      end if


      do i=1,rcnt+ntts+pcnt
         dir="densum"//'/'//trim(sdir(i))
         inquire(file=dir,exist=there)
         if(.not.there)then
            command='mkdir '//trim(dir)
            call system(command)
         end if
         fname=trim(dir)//'/'//"densum.dat"
         open(unit=tunit,file=fname)
         write(tunit,*)trim(title)
         write(tunit,*)trim(title),'-',trim(formulal(i)),'_'
         write(tunit,400)ndof(i),0,'HAR','cm-1'
         write(tunit,300)10.,1000,5000,etopuser
         do j=1,ndof(i)
            if((idofl(i,j).eq.'kro').or.(idofl(i,j).eq.'jro'))then
               write(tunit,100)j,'qro',wel(i,j),anhl(i,j),
     $                                                  ngl(i,j)

            else
               write(tunit,100)j,idofl(i,j),wel(i,j),anhl(i,j),
     $                                                  ngl(i,j)
               if(idofl(i,j).eq.'hrd')then
                  write(tunit,200)vhrl(i,j),nsvl(i,j),phavl(i,j),
     $                                       (cvl(i,j,k),k=1,ncvl(i,j))
                  write(tunit,200)bhrl(i,j),nsbl(i,j),phabl(i,j),
     $                                       (cbl(i,j,k),k=1,ncbl(i,j))
               end if
            end if
         end do
         write(tunit,*)
      end do

      end subroutine
            

