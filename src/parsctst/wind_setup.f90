SUBROUTINE wind_setup(windbal_key)
    USE decl_alloc
    implicit none
    include "mpif.h"
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!    
!    LICENSE NOTICE
!
!    Copyright (C) 2019 Michele Ceotto, Chiara Aieta, Fabio Gabas, 
!                       Thanh Lam Nguyen, and John R. Barker
!
!    Contact:
!
!    Michele Ceotto  (email: michele.ceotto@unimi.it)
!    Dipartimento di Chimica
!    Università degli Studi di Milano
!    via Golgi 19, Milano 20133, ITALY
!
!    or:
!
!    John R. Barker   (email: jrbarker@umich.edu)
!    Department of Atmospheric, Oceanic, and Space Sciences
!    College of Engineering
!    University of Michigan
!    Ann Arbor, MI 48109
!
!    http://clasp-research.engin.umich.edu/multiwell
!    This program is free software; you can redistribute it and/or
!    modify it under the terms of the GNU General Public License (version 2)
!    as published by the Free Software Foundation.
!   
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!   
!    See the 'ReadMe' file for a copy of the GNU General Public License,
!    or contact:
!   
!    Free Software Foundation, Inc.
!    59 Temple Place - Suite 330
!    Boston, MA 02111-1307, USA.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                 
!                         PROGRAM parsctst    
!
!                               by      
!
!           by Michele Ceotto, Chiara Aieta, Fabio Gabas,
!               Thanh Lam Nguyen, and John R. Barker    
!                                                                  
!           ***PARalle Semi-Classical Transition State Theory***
!
!                             based on
!
!          The theory of W. H. Miller and coworkers* and the
!          parallel implementation** of the Wang-Landau algorithms 
!          for densities of states***
!
!    Literature Citations:                                                    
!    *Semi-Classical Transition State Theory
!    W. H. Miller, J. Chem. Phys. 62, 1899-1906 (1975).
!    W. H. Miller, Faraday Discuss. Chem. Soc. 62, 40-46 (1977).
!    W. H. Miller, R. Hernandez, N. C. Handy, D. Jayatilaka, and A. Willets,
!      Chem. Phys. Letters 172, 62-68 (1990).
!    R. Hernandez and W. H. Miller, Chem. Phys. Lett. 214, 129-136 (1993).
!    J. F. Stanton, J. Phys. Chem. Lett. 7, 2708-2713 (2016)
!
!    **Semi-Classical Transition State Theory parallel implementation
!    C. Aieta, F. Gabas and M. Ceotto, J. Chem. Theory Comput.,
!      15, 2142−2153 (2019).
!                                                                  
!    ***Density of states algorithms
!    F. Wang and D. P. Landau, Phys. Rev. Letters 86, 2050-2053 (2001).           
!    M. Basire, P. Parneix, and F. Calvo, J. Chem.Phys. 129, 081101 (2008).    
!    T. L. Nguyen and J.R. Barker, J. Phys. Chem. A., 114, 3718–3730 (2010).
!    C. Aieta, F. Gabas and M. Ceotto, J. Phys. Chem. A., 120(27), 4853-4862 (2016).
!                                                                  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
      INTEGER(4) :: i
      CHARACTER(4) :: windbal_key

      call MPI_COMM_RANK( MPI_COMM_WORLD, my_rank, err )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, num_procs, err )

           
      CALL ucase ( windbal_key )

      IF (windbal_key .EQ. 'COST') THEN
         !     Each window must have an integer number of bins
         ngrains_per_chunk = CEILING( DBLE(ngrains) / DBLE(nwind) )
 
         ngrains_corr = ngrains_per_chunk(1) * nwind
         IF ( ngrains_corr /= ngrains ) THEN
            Emax = (ngrains_corr-1) * Egrain1
            ngrains = ngrains_corr
            IF (my_rank == 0) THEN
               WRITE(61,*) '************************'
               WRITE(61,*) 'WARNING! New parameters:'
               WRITE(61,*) 'ngrains: ', ngrains
               WRITE(61,*) 'Emax: ', Emax
               WRITE(61,*) '************************'
            END IF
         END IF
!
         half_perc_wind_overlap = perc_wind_overlap / 2.d0
         chunk = ( (Emax + Egrain1) - Emin ) / nwind

!        If there's just one window do not add the overlap
         IF ( nwind /= 1) THEN
            ngrains_to_add = NINT( ngrains_per_chunk * half_perc_wind_overlap / 100.d0 )
         ELSE
            ngrains_to_add = 0
         END IF

         ngrains_per_wind = ngrains_per_chunk + ngrains_to_add

!        The last window exceedes the high energy limit 
!        The required energy range will be recovered at the end
         DO i=1, nwind
            lowbound(i) = Emin + (i-1)*chunk(i)
            upbound(i) = Emin + i*chunk(i) + ngrains_to_add(i) * Egrain1
         END DO

         Emin_global = Emin
         Emax_global = upbound(nwind)

         ngrains_per_wind_max = ngrains_per_wind(1)
         ngrains_to_add_max = ngrains_to_add(1)
   
      ELSE IF (windbal_key .EQ. 'LOW ') THEN

         IF (nwind == 1 ) THEN

            ngrains_per_chunk = CEILING( DBLE(ngrains) / DBLE(nwind) )
            ngrains_to_add = 0
            half_perc_wind_overlap = perc_wind_overlap / 2.d0
            chunk = ( (Emax + Egrain1) - Emin ) / nwind

         ELSE IF (nwind == 2) THEN

             ngrains_per_chunk(1) = CEILING ( DBLE(ngrains) / DBLE((nwind *2)))
             ngrains_per_chunk(2) = ngrains - CEILING ( DBLE(ngrains) /DBLE((nwind * 2)))
             half_perc_wind_overlap(1) = perc_wind_overlap / 4.d0
             half_perc_wind_overlap(2) = perc_wind_overlap / 2.d0
             DO i = 1, nwind
                chunk(i) = ngrains_per_chunk(i)*Egrain1
                ngrains_to_add(i) = NINT( ngrains_per_chunk(i) * half_perc_wind_overlap(i) / 100.d0 )
             END DO

         ELSE

             ngrains_per_chunk_small = CEILING ( DBLE(ngrains) / DBLE((nwind * 2)))
             ngrains_per_chunk_big   = CEILING ( DBLE((ngrains - ngrains_per_chunk_small*2)) / DBLE((nwind -2)))

             ngrains_per_chunk = ngrains_per_chunk_big
             ngrains_per_chunk(1) = ngrains_per_chunk_small
             ngrains_per_chunk(2) = ngrains_per_chunk_small

             ngrains_corr = ngrains_per_chunk_big * (nwind -2)

             IF ( ngrains_corr /= (ngrains - ngrains_per_chunk_small*2) ) THEN
                Emax = (ngrains_corr + ngrains_per_chunk_small*2 - 1) *Egrain1
                ngrains = ngrains_corr + ngrains_per_chunk_small*2
                IF (my_rank == 0) THEN
                   WRITE(61,*) '************************'
                   WRITE(61,*) 'WARNING! New parameters:'
                   WRITE(61,*) 'ngrains: ', ngrains
                   WRITE(61,*) 'Emax: ', Emax
                   WRITE(61,*) '************************'
                END IF
             END IF

             half_perc_wind_overlap = perc_wind_overlap / 2.d0
             half_perc_wind_overlap(1) = perc_wind_overlap / 4.d0
             half_perc_wind_overlap(2) = perc_wind_overlap / 4.d0
             DO i = 1, nwind
                chunk(i) = ngrains_per_chunk(i)*Egrain1
                ngrains_to_add(i) = NINT( ngrains_per_chunk(i) *half_perc_wind_overlap(i) / 100.d0 )
             END DO
         END IF

         DO i = 1, nwind
            ngrains_per_wind(i) = ngrains_per_chunk(i) + ngrains_to_add(i)
         END DO

         !     The last window exceedes the high energy limit
         !     The required energy range will be recovered at the end
         lowbound(1) = Emin
         DO i = 2, nwind
            lowbound(i) = lowbound(i-1) + ngrains_per_chunk(i-1)*Egrain1
         END DO

         DO i = 1, nwind
            upbound(i) = lowbound(i) + ngrains_per_wind(i) * Egrain1
         END DO

         Emin_global = Emin
         Emax_global = upbound(nwind)

         ngrains_per_wind_max = ngrains_per_wind(nwind)
         ngrains_to_add_max = ngrains_to_add(nwind)

      ELSE IF (windbal_key .EQ. 'HIGH') THEN      

         IF (nwind == 1 ) THEN

            ngrains_per_chunk = CEILING( DBLE(ngrains) / DBLE(nwind) )
            ngrains_to_add = 0
            half_perc_wind_overlap = perc_wind_overlap / 2.d0
            chunk = ( (Emax + Egrain1) - Emin ) / nwind

         ELSE IF (nwind == 2) THEN

             ngrains_per_chunk(2) = CEILING ( DBLE(ngrains) / DBLE((nwind *2)))
             ngrains_per_chunk(1) = ngrains - CEILING ( DBLE(ngrains) /DBLE((nwind * 2)))
             half_perc_wind_overlap(2) = perc_wind_overlap / 4.d0
             half_perc_wind_overlap(1) = perc_wind_overlap / 2.d0
             DO i = 1, nwind
                chunk(i) = ngrains_per_chunk(i)*Egrain1
                ngrains_to_add(i) = NINT( ngrains_per_chunk(i) * half_perc_wind_overlap(i) / 100.d0 )
             END DO

         ELSE

             ngrains_per_chunk_small = CEILING ( DBLE(ngrains) / DBLE((nwind * 2)))
             ngrains_per_chunk_big   = CEILING ( DBLE((ngrains - ngrains_per_chunk_small*2)) / DBLE((nwind -2)))

             ngrains_per_chunk = ngrains_per_chunk_big
             ngrains_per_chunk(nwind-1) = ngrains_per_chunk_small
             ngrains_per_chunk(nwind) = ngrains_per_chunk_small

             ngrains_corr = ngrains_per_chunk_big * (nwind -2)

             IF ( ngrains_corr /= (ngrains - ngrains_per_chunk_small*2) ) THEN
                Emax = (ngrains_corr + ngrains_per_chunk_small*2 - 1) *Egrain1
                ngrains = ngrains_corr + ngrains_per_chunk_small*2
                IF (my_rank == 0) THEN
                   WRITE(61,*) '************************'
                   WRITE(61,*) 'WARNING! New parameters:'
                   WRITE(61,*) 'ngrains: ', ngrains
                   WRITE(61,*) 'Emax: ', Emax
                   WRITE(61,*) '************************'
                END IF
             END IF

             half_perc_wind_overlap = perc_wind_overlap / 2.d0
             half_perc_wind_overlap(nwind-1) = perc_wind_overlap / 4.d0
             half_perc_wind_overlap(nwind) = perc_wind_overlap / 4.d0
             DO i = 1, nwind
                chunk(i) = ngrains_per_chunk(i)*Egrain1
                ngrains_to_add(i) = NINT( ngrains_per_chunk(i) *half_perc_wind_overlap(i) / 100.d0 )
             END DO
         END IF
         
         DO i = 1, nwind
            ngrains_per_wind(i) = ngrains_per_chunk(i) + ngrains_to_add(i)
         END DO

         !     The last window exceedes the high energy limit
         !     The required energy range will be recovered at the end
         lowbound(1) = Emin
         DO i = 2, nwind
            lowbound(i) = lowbound(i-1) + ngrains_per_chunk(i-1)*Egrain1
         END DO

         DO i = 1, nwind
            upbound(i) = lowbound(i) + ngrains_per_wind(i) * Egrain1
         END DO

         Emin_global = Emin
         Emax_global = upbound(nwind)

         ngrains_per_wind_max = ngrains_per_wind(1)
         ngrains_to_add_max = ngrains_to_add(1)

      END IF


! print windows grains and boundaries

      DO i = 1, nwind
         ngrains_per_wind(i) = ngrains_per_chunk(i) + ngrains_to_add(i)
         IF(my_rank .eq. 0) WRITE(61,*) 'wind', i, 'ngrains_per_wind', ngrains_per_wind(i)
         IF(my_rank .eq. 0) WRITE(61,*) 'wind', i, 'ngrains_per_chunk', ngrains_per_chunk(i)
      END DO

      IF(my_rank.eq.0) then
           WRITE(61,*) 'Lowbound array:'
           DO i=1, nwind
              WRITE(61,*) i, lowbound(i)
           END DO

           WRITE(61,*) 'Upbound array:'
           DO i=1, nwind
              WRITE(61,*) i, upbound(i)
           END DO
           CLOSE(61)
      END IF


 END SUBROUTINE 
