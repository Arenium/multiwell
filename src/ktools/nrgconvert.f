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

      subroutine nrgconvert(array1,n,unitsin,unitsout)
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      double precision temp1, temp2, array1(n)
      character unitsin*4,unitsout*4
      save

c  convert unitsin to cm-1
c
c  conversion factors are stored as parameters in 'include.inc'

      call lowerc(unitsin)
      select case(unitsin)
         case ("amua")
            do i= 1, n
               array1(i) = amua2nu/array1(i)
            end do
         case ("gcmcm")
            do i= 1, n
               array1(i) = gmcm2nu/array1(i)
            end do
         case ("cm-1")
            do i= 1, n
               array1(i) = array1(i)
            end do
         case ("mhz")
            do i= 1, n
               array1(i) = array1(i)*(1.0d+06/c)
            end do
         case ("ghz")
            do i= 1, n
               array1(i) = array1(i)*(1.0d+09/c)
            end do
         case ("kcal")
            do i= 1, n
               array1(i) = array1(i)*kcal2nu
            end do
         case ("kj")
            do i= 1, n
               array1(i) = array1(i)*(kcal2nu/4.184d+00)
            end do
      end select

      
c  convert unitsin to unitsout
c
c
c      unitsout="cm-1"
      call lowerc(unitsout)
      select case(unitsout)
         case ("amua")
            do i= 1, n
               array1(i) = amua2nu/array1(i)
            end do
         case ("gcmcm")
            do i= 1, n
               array1(i) = gmcm2nu/array1(i)
            end do
         case ("cm-1")
            do i= 1, n
               array1(i) = array1(i)
            end do
         case ("mhz")
            do i= 1, n
               array1(i) = array1(i)/(1.0d+06/c)
            end do
         case ("ghz")
            do i= 1, n
               array1(i) = array1(i)/(1.0d+09/c)
            end do
         case ("kcal")
            do i= 1, n
               array1(i) = array1(i)/kcal2nu
            end do
         case ("kj")
            do i= 1, n
               array1(i) = array1(i)/(kcal2nu/4.184d+00)
            end do

      end select
      

      return


      end subroutine nrgconvert

