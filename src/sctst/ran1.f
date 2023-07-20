**==ran1.spg  processed by SPAG 5.11R  at 20:12 on 26 Mar 2001
 
 
      REAL(8) FUNCTION RAN1(idum)
C---------------------------------------------------------------------------------
C      "Minimal" random number generator of Park and Miller with Bays-Durham shuffle
C      and added safeguards.  Returns a uniform random deviate between 0.0 and 1.0
C      (exclusive of the endpoint values).  Call with idum a negative INTEGER(4)(KIND=4) to
C      initialize; thereafter do not alter idum between successive deviates in a
C      sequence.  RNMX should approximate the largest floating value that is less
C      than 1.
C
C      W. H. Press and S. A. Teukolsky, Computers in Physics, 6(5), 522-4 (1992).
C---------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(4) idum , IA , IM , IQ , IR , NTAB , NDIV
      REAL(KIND=8) AM , EPS , RNMX , T
      PARAMETER (IA=16807,IM=2147483647,AM=1.0D+00/IM,IQ=127773,IR=2836,
     &        NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=2.3D-16,RNMX=1.0D+00-EPS)
      INTEGER(4) j , k , iv(NTAB) , iy 
      SAVE iv , iy
      DATA iv/NTAB*0/ , iy/0/
      IF ( idum.LE.0 .OR. iy.EQ.0 ) THEN
                                        ! Initialize
         idum = max(-idum,1)            ! Be sure to prevent idum=0
         DO 50 j = NTAB + 8 , 1 , -1    ! Load shuffle table after 8 warm-ups
            k = idum/IQ
            idum = IA*(idum-k*IQ) - IR*k
            IF ( idum.LT.0 ) idum = idum + IM
            IF ( j.LE.NTAB ) iv(j) = idum
 50      CONTINUE
         iy = iv(1)
      ENDIF
      k = idum/IQ                       ! Start here when not initializing
      idum = IA*(idum-k*IQ) - IR*k      ! Compute idum=mod(IA*idum,IM) without
      IF ( idum.LT.0 ) idum = idum + IM !    overflows by Schrage's Method
      j = 1 + iy/NDIV                   ! Will be in range 1:NTAB
      iy = iv(j)                        ! Output previously stored value
      iv(j) = idum                           ! Refill shuffle table
      T = AM*iy
      RAN1 = min(T,RNMX)        ! Because users don't expect endpoint values
c11      format(d90.50)
c      write(6,11)RAN1
      RETURN
      END
