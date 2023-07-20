c
c john r. barker
c jrbarker@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c (734) 763 6239
c
c Copyright 2019
c 
c Jason A. Sonk and john R. Barker
c University of Michigan
c Ann Arbor, MI 48109
c
c Contact: jsonk@umichedu and jrbarker@umich.edu
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
!
!    INCLUDE file for declaring variables and parameters in program KTOOLS.
!
      IMPLICIT NONE

      integer ::  i, j, k, l, m, n, iv, ib
c
c i/o units
c
      integer, parameter :: iunit = 40
      integer, parameter :: lunit = 41		 ! log file
      integer, parameter :: efunit = 43
      integer, parameter :: canunit = 42
      integer, parameter :: veffout = 44
      integer, parameter :: densunit = 45
      integer, parameter :: nminsunit = 4
      integer, parameter :: junit = 47
      integer, parameter :: checkunit = 48
c
c noise parameters
c
      real (8) , parameter :: canonh = 2.0d0
      real (8) , parameter :: microh = 1.5d0
c
c found in nrgconvert
c
      real (8), parameter :: converge = 1.0d-16
      real (8), parameter :: h = 6.62606957d-34
      real (8), parameter :: c = 29979245800.0d0
      real (8), parameter :: kb = 0.69503476d0
      real (8), parameter :: kbj = 1.3806488d-23
      real (8), parameter :: avo = 6.02214129d+23
      real (8), parameter :: kcaltoj = 4.184d+03
      real (8), parameter :: pi = 2.0d0*dacos(0.0d0)
      real (8), parameter :: amu2kg= (1.0d0/avo)*(1.0d0/1.0d+03)
      real (8), parameter :: 
     &                 kcal2nu=(1.0d0/avo)*kcaltoj*(1.0d0/h)*(1.0d0/c)
      real (8), parameter :: amua2nu=((h*avo*1.0d+23)/(8.0d0*c*(pi**2)))
      real (8), parameter :: gmcm2nu=((h*1.0d+7)/(8.0d0*c*(pi**2)))
      real (8), parameter :: ghz2nu=((1.0d+09)/c)
      real (8), parameter :: mhz2nu=((1.0d+06)/c)
      real (8), parameter :: rgasjol = 8.3144621d+00
      real (8), parameter :: rgascal = rgasjol/kcaltoj
      real (8), parameter :: rgascm1 = rgasjol/(h*c*avo)

c
c part of hrd
c
      real (8), dimension (datpts*10) :: cbb
      real (8), dimension (datpts*10) :: cvv
      real (8), dimension (datpts*20) :: ev

      
      
