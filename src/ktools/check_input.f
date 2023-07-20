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
      subroutine check_input
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      integer errnum
      character frmt*60
 100  format('variational ts:  ',a,'  run started')
 101  format(1x,6(i2,10x))
 102  format("structure #",i4,1x,"dof #",1x,i4,a4," is not recognized.")
 200  format("     ",a,2x,i4,2x,f6.2,10(2x,f10.2))
 99   format(1x,6(i2,10x))
 98   format(a,2(2x,i2),2(2x,f10.3))
 97   format(i3,1x,a3,2(2x,f10.3),2x,i5)
 96   format(a3,1x,a3,2(2x,a9),2x,a5)
 95   format(3(a7,4x))
 94   format(3(i7,4x))

      call sestamp('check_input',1)

c  check if user supplied degrees-of-freedom are recognized

      errnum = 0 
      do i=1,rcnt+ntts+pcnt
         do j=1, ndof(i)
            select case(idofl(i,j))
               case('vib')
                  continue
               case('hra')
                  continue
               case('hrb')
                  continue
               case('hrc')
                  continue
               case('hrd')
                  continue
               case('kro')
                  continue
               case('jro')
                  continue
               case('qro')
                  continue
               case default
                  write(*,102)i,j,idofl(i,j)
                  errnum = errnum + 1
            end select
         end do
      end do

c  write copy of structural data to log file

      write(lunit,95)'# of r','# of ts','# of p'
      write(lunit,94)rcnt,ntts,pcnt
      write(lunit,*)
      write(lunit,*)

      do i=1,rcnt+ntts+pcnt
        if(i.le.rcnt) then
         write(lunit,98)'reactant #',i,ndof(i),delh(i),distl(i)
        else if(i.le.rcnt+ntts) then
         write(lunit,98)'ts #',i-rcnt,ndof(i),delh(i),distl(i)
        else
         write(lunit,98)'product #',i-rcnt-ntts,ndof(i),delh(i),distl(i)
        end if
        write(lunit,96)'#','dof','freq','anh','degen'
        write(lunit,*)('-',j=1,36)
        do j=1,ndof(i)
           write(lunit,97)j,idofl(i,j),wel(i,j),anhl(i,j),ngl(i,j)
           if(idofl(i,j).eq.'hrd')then
              write(lunit,200)vhrl(i,j),nsvl(i,j),phavl(i,j),
     $                                       (cvl(i,j,k),k=1,ncvl(i,j))
              write(lunit,200)bhrl(i,j),nsbl(i,j),phabl(i,j),
     $                                       (cbl(i,j,k),k=1,ncbl(i,j))
           end if
        end do
        write(lunit,*)
        write(lunit,*)
      end do

c exit if undefined dof is present

      if(errnum.ne.0)then
         write(*,*)errnum,' undefined dof found. exiting.'
         stop
      end if

300   format(i2,' reactants and ', i2,' products found.')
      if(rcnt+pcnt.gt.3)then
         write(*,*)'*** FATAL, exiting ***'
         write(*,300)rcnt,pcnt
         write(*,*)'This is not a recombination/dissociation reaction.'
         write(*,*)'  '
         stop
      end if
      
      

      call sestamp('check_input',2)

      return

      end subroutine

