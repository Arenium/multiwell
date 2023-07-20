c
c	Date & time subroutine
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Copyright (C) 2010 John R. Barker
c
c John R. Barker
c jrbarker@umich.edu
c University of Michigan
c Ann Arbor, MI 48109-2143
c (734) 763 6239
c
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License (version 2)
c as published by the Free Software Foundation.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
c GNU General Public License for more details.
c
c See the 'ReadMe' file for a copy of the GNU General Public License,
c or contact:
c
c Free Software Foundation, Inc.
c 59 Temple Place - Suite 330
c Boston, MA 02111-1307, USA.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      SUBROUTINE DateTime(KUNIT)
c
c      Date & time subroutine
c
      IMPLICIT NONE
      INTEGER(4) :: KUNIT
      REAL(8) :: duration
      CHARACTER(len=10) date, time, zone
      INTEGER values(8)

      call date_and_time(date,time,zone,values)
      WRITE (KUNIT,99001) values(1), values(2), values(3), values(5), 
     &                    values(6), values(7), values(8)

      CALL CPU_TIME(duration)
      WRITE (KUNIT,99002) duration, duration/3600.
      
      RETURN
      
99001 FORMAT (/'....................................................',/,
     &     10x, I4,'/',I2.2,'/',I2.2,3x,I2.2,':',I2.2,':',
     &          I2.2,'.',I3.3,'   (Local Time)',/
     &     10x, 'year/mm/dd   hr:mi:second'   )
 
99002 FORMAT (17x, 'CPU :   ',f10.3,' s (',f7.3,' hr)',/
     &        '....................................................')

      RETURN
      END
