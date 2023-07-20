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

      subroutine stamp(iout,len)
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      integer iout,len


9911  FORMAT (
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'//,
     &'                    ',A6,'-',A6,' (',A8,')',//,
     &'  Copyright (C) 2016 Jason A. Sonk, Thanh Lam Nguyen, and'
     &' John R. Barker',//
     &'        CONTACT:         John R. Barker',/
     &'                      (jrbarker@umich.edu)',/
     &'                     University of Michigan',/
     &'                  Ann Arbor, Michigan 48109-2143',//
     &'         http://clasp-research.engin.umich.edu/multiwell/',//
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%',//'Suggested Literature Citations:'//4x,
     &'a) J.R. Barker, T.L. Nguyen, J.F. Stanton, C. Aieta, M. Ceotto,',
     &/7x,'F. Gabas, T.J.D. Kumar, C.G.L.Li, L.L.Lohr, A. Maranzana,',
     &/7x,'N.F. Ortiz, J.M. Preses, J.M. Simmie, J.A. Sonk, and ',
     &'P.J. Stimac; ',
     &/7x,'MultiWell-',A7,' Software Suite; University of Michigan, ',
     &/7x,'Ann Arbor, Michigan, USA, ',A8,'; ',
     &/7x,'http://clasp-research.engin.umich.edu/multiwell/.',//4x
     &'b) F. Wang and D.P. Landau, Phys. Rev. Letters 86, ',/7x,
     &'2050-2053 (2001).',//4x,
     &'c) M. Basire, P. Parneix, and F. Calvo, J. Chem. Phys. 129,',/7x,
     &'081101 (2008).',//4x,
     &'d) T.L. Nguyen and J.R. Barker, J. Phys. Chem. A., 114,',/7x,
     &'3718-3730 (2010).',
     &//,
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'/)


99016 FORMAT (
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'//8x,
     &'                              Jason A. Sonk & John R. Barker'/8x,
     &'        ',A6,'-',A6,1x'        University of Michigan'/8x,
     &'                              Ann Arbor, MI 48109-2143'/8x,
     &'          ',A8,'              Contact:',/8x,
     &'                                jrbarker@umich.edu'/8x,
     &'                                (734) 763 6239'//8x,
     &'      http://aoss.engin.umich.edu/multiwell/'//,
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'//8x,
     &'Suggested Literature Citation:'//8x,
     &'a) J.R. Barker, T.L. Nguyen, J.F. Stanton, C. Aieta, M. Ceotto,',
     &/7x,'F. Gabas, T.J.D. Kumar, C.G.L.Li, L.L.Lohr, A. Maranzana,',
     &/7x,'N.F. Ortiz, J.M. Preses, J.M. Simmie, J.A. Sonk, and ',
     &'P.J. Stimac; ',
     &/7x,'MultiWell-',A7,' Software Suite; University of Michigan, ',
     &/7x,'Ann Arbor, Michigan, USA, ',A8,'; ',
     &/7x,'http://clasp-research.engin.umich.edu/multiwell/.',//8x,
     &'b) J.R. Barker, Int. J. Chem. Kinetics, 33, 232-45 (2001)',//8x,
     &//,
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'/)


99002 FORMAT (
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'//,
     &'                               Jason A. Sonk & John R. Barker',/,
     &'       ',A6,'-',A6, '           University of Michigan',/,
     &'                               Ann Arbor, MI 48109-2143',/,
     &'         ',A8,  '              Contact:',/,
     &'                                 jrbarker@umich.edu',/,
     &'                                 (734) 763 6239',//,
     &'            http://aoss.engin.umich.edu/multiwell/'/
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'//)



      if(len.eq.1)then
         write(iout,9911)prog,aversion,adate,aversion,adate
      else
         write(iout,99002)prog,aversion,adate
      end if

      end subroutine

