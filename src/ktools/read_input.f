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
      subroutine read_input
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      character logfile*80
      logical there
      character date*8,time*10,zone*5
      character (len=1) ADUM
      integer line
      real*8   tngl
      integer values(8)
 100  format(a,': ',a,' run started.')
 200  format('*** ',a,' not found ***')
 300  format('file:',a,' found')
 400  format(5(i2.2))
99114 FORMAT (A180)

c  prepare unique time stamp for file identification

      call date_and_time(date,time,zone,values)
      write(tstmp,400)values(1)-2000,values(2),values(3),values(5),
     $                values(6)
      write(*,*)'run id: ',tstmp


c  get user supplied input file name, otherwise use ktools.dat default

      call get_command_argument(1,inputfile)
      inputfile=trim(inputfile)
      if(inputfile.eq.'')then
            inputfile='ktools.dat'
            logfile='logfile.log'
      else
            inputfile=inputfile
            fileroot=inputfile(1:len_trim(inputfile)-4)
            logfile=trim(fileroot)//".log"
      end if

c  check to see if inputfile exists in directory

      inquire(file=inputfile,exist=there)
      if(.not.there)then
            write(*,200)trim(inputfile)
            stop
      else
            write(*,300)trim(inputfile)
            write(*,*)'reading input'
      end if

c  open input and log files

      open(unit=iunit,file=inputfile)
      open(unit=lunit,file=logfile)

c  begin reading data

      read(iunit,99114)title
      ADUM = '!'
      line = 0
      DO WHILE ( ADUM .EQ. '!' )					! read up to 20 comment lines starting with "!"
         READ (iunit, *) ADUM
         IF ( ADUM .EQ. '!' ) THEN
            line = line + 1
            BACKSPACE (iunit)
            READ (iunit,99114) commentline( line )
         ENDIF
      END DO
      nlines = line								! number of comment lines (nlines ≤ 20)
      BACKSPACE (iunit)
      
      call stamp(lunit,1)
      write(lunit,*)
      write(lunit,*)
      write(lunit,100)prog,trim(title)
      call dnt(lunit)
      write(lunit,*)'run id: ',tstmp
      write(lunit,*)

      read(iunit,*)eunits,sunits
c      read(iunit,*)whatdo,backup				! read keywords
      read(iunit,*)whatdo				! read keywords
      read(iunit,*)etopuser,de
      read(iunit,*)jtopuser,dj
      read(iunit,*)dimax1,disize,demax2
      read(iunit,*)nt
      read(iunit,*)(temps(i),i=1,nt)
      read(iunit,*)rcnt,ntts,pcnt
      
c  find maximum temperature

      maxtemp=temps(1)
      do i=1,nt
         if(temps(i).gt.maxtemp)maxtemp=temps(i)
      end do

c  find maximum number of bins above TS energy

      binmax=nint(etopuser/de)+1

c  trim unit lengths and convert to lowercase

      eunits=trim(eunits)
      sunits=trim(sunits)
      whatdo=trim(whatdo)
      call lowerc(eunits)
      call lowerc(sunits)
      call lowerc(whatdo)
      call lowerc(backup)

c  read structure data

      do 2222 i=1, rcnt+ntts+pcnt
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

      return

      end subroutine
