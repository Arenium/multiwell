      SUBROUTINE CalPN(PN, DE, omega, xFF)
c
c         PN = Semi-classical tunneling transmission probability
c         DE = deltaE (Eq. 7b in Miller et al. 1990) ***with new sign convention***
c         omega = imaginary frequency, including coupling with orthogonal DOF: omega(F) in Eq. 7c (Miller et al. 1990)
c         xFF = diagonal anharmonic coupling coefficient for rxn coordinate
c
c         Latest revision: 12/2015
c
      Double precision PN, THETA, DE, omega, xFF
      Double precision PI, tp , D
      PI=acos(-1.0d0)

      D = (omega**2) / (4.0d0*xFF)                ! when xFF<0, -D = barrier height from VPT2
      IF ( omega.GT. 0.d+00 ) THEN
        IF( D.GE.0.0 .AND. DE.GE.D) THEN          ! failure of VPT2 when DE too positive 
          PN = 1.0d+00
        ELSEIF( D.LE.0.0 .AND. DE.LE.D) THEN      ! failure of VPT2 when DE too negative 
          PN = 0.0d+00
        ELSE
          tp = sqrt ( 1.0d0 - DE/D )
          THETA = -(PI*DE*2.0d0/omega)/( 1.0d0 + tp )      ! Eq. 7a in Miller et al. 1990 with new sign of DE
          PN = 1.0d0 / (1.0d0 + exp(2.0d0*THETA) )         ! Eq. 5b in Miller et al. 1990
        ENDIF
      ELSE
        IF (DE .GT. 0.0d0) THEN             ! logic: if omega .LE. 0, then barrier is very thick and hence no tunneling
          PN = 1.0d+00
        ELSE
          PN = 0.0d+00
        ENDIF
      ENDIF

      RETURN
      END 
