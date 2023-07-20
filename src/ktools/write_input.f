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

      subroutine write_input
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      character logfile*80
      integer line
      logical there
 200  format(i4,1x,i2)
 300  format(a)
 400  format(5(i2.2))
 500  format(a,2x,a)
 600  format(i6,1x,f8.2)
 700  format(2(a,2x),2(f10.2,2x))
 800  format(100(f8.2,1x))
 900  format(3(i3,2x))
1000  format(2(f5.2,2x),i3,2x)
1100  format(f5.2,2x,i3,2x)
1200  format(i3,2x,a,2x,f10.2,2x,f10.2,2x,i3,2x)

c  get user supplied input file name, otherwise use ktools.dat default

      logfile=trim(tstmp)//".bck"

c  open input and log files

      open(unit=iunit,file=logfile)

c  begin writing data

      write(iunit,300)title
      DO line = 1, nlines
        write(iunit,*) commentline(line)
      END DO

      write(iunit,500)eunits,sunits
      write(iunit,500)whatdo
      write(iunit,600)etopuser,de
      write(iunit,200)jtopuser,dj
      write(iunit,200)nt
      write(iunit,800)(temps(i),i=1,nt)
      write(iunit,900)rcnt,ntts,pcnt
      write(iunit,*)

c  write structure data

      do 2222 i=1, rcnt+ntts+pcnt
            write(iunit,700)reprod(i),trim(molname(i)),delh(i),distl(i)
            write(iunit,*)formulal(i)
            write(iunit,*)title1(i)
            write(iunit,*)title2(i)
            write(iunit,*)title3(i)
            write(iunit,1000)sym(i),sopt(i),nele(i)
            do j=1,nele(i)
               write(iunit,1100)elev(i,j),gele(i,j)
            end do

            write(iunit,*)ndof(i),keyword(i,1),keyword(i,2)

            do j=1,ndof(i)
               write(iunit,1200)j,idofl(i,j),wel(i,j),anhl(i,j),ngl(i,j)
               if(idofl(i,j).eq.'hrd')then
                  ncvl(i,j)=dint(wel(i,j))
                  ncbl(i,j)=dint(anhl(i,j))
                  write(iunit,*)vhrl(i,j),nsvl(i,j),phavl(i,j),
     $                                      (cvl(i,j,iv),iv=1,ncvl(i,j))
                  write(iunit,*)bhrl(i,j),nsbl(i,j),phabl(i,j),
     $                                      (cbl(i,j,ib),ib=1,ncbl(i,j))
               end if
            end do

      write(iunit,*)
 2222 continue
      close(iunit)

      return

      end subroutine
