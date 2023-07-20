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
      subroutine sort_input(sortarr,sortcnt)
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      integer sortarr(rcnt+ntts+pcnt)
      integer sortcnt
      real*8  tngl
      logical there
 100  format(a,': ',a,' run started.')
 200  format('*** ',a,' not found ***')
 300  format('file:',a,' found')

      sortcnt = sortcnt + 1

c  open input and log files

      open(unit=iunit,file=inputfile)

c  begin reading data

      read(iunit,*)title
      read(iunit,*)eunits,sunits
      read(iunit,*)whatdo,backup
      read(iunit,*)etopuser,de
      read(iunit,*)jtopuser,dj
      read(iunit,*)dimax1,disize,demax2
      read(iunit,*)nt
      read(iunit,*)(temps(i),i=1,nt)
      read(iunit,*)rcnt,ntts,pcnt

c  trim unit lengths and convert to lowercase

      eunits=trim(eunits)
      sunits=trim(sunits)
      whatdo=trim(whatdo)
      call lowerc(eunits)
      call lowerc(sunits)
      call lowerc(whatdo)
      call lowerc(backup)

c  read structure data

      do 2222 k=1, rcnt+ntts+pcnt
            i=sortarr(k)
            read(iunit,*)reprod(i),molname(i),delh(i),distl(i)
            read(iunit,*)formulal(i)
            read(iunit,*)title1(i)
            read(iunit,*)title2(i)
            read(iunit,*)title3(i)
            read(iunit,*)sym(i),sopt(i),nele(i)
            do j=1,nele(i)
               read(iunit,*)elev(i,j),gele(i,j)
            end do

            read(iunit,*)ndof(i),keyword(i,1),keyword(i,2)

            do j=1,ndof(i)
               read(iunit,*)mode(i,j),idofl(i,j),wel(i,j),anhl(i,j),tngl
               ngl(i,j)=int(tngl)
               call lowerc(idofl(i,j))
               if(idofl(i,j).eq.'hrd')then
                  ncvl(i,j)=dint(wel(i,j))
                  ncbl(i,j)=dint(anhl(i,j))
                  read(iunit,*)vhrl(i,j),nsvl(i,j),phavl(i,j),
     $                                      (cvl(i,j,iv),iv=1,ncvl(i,j))
                  read(iunit,*)bhrl(i,j),nsbl(i,j),phabl(i,j),
     $                                      (cbl(i,j,ib),ib=1,ncbl(i,j))
               end if
            end do

 2222 continue
      close(iunit)

      call dsort(sortcnt)

      return

      end subroutine
