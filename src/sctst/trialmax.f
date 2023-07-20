      FUNCTION trialmax( KEYWORD , ngrains )
c
c   Trials per energy grain
c
      CHARACTER(len=6) KEYWORD
      INTEGER ngrains
      REAL(8) trialmax
      SAVE
      
      IF    ( KEYWORD .EQ. 'extra' ) THEN
        trialmax = 1.e+6
      ELSEIF ( KEYWORD .EQ. 'Extra' ) THEN
        trialmax = 1.e+6
      ELSEIF ( KEYWORD .EQ. 'EXTRA' ) THEN
        trialmax = 1.e+6
      ELSEIF ( KEYWORD .EQ. 'best' ) THEN
        trialmax = 1.e+5
      ELSEIF ( KEYWORD .EQ. 'Best' ) THEN
        trialmax = 1.e+5
      ELSEIF ( KEYWORD .EQ. 'BEST' ) THEN
        trialmax = 1.e+5
      ELSEIF ( KEYWORD .EQ. 'better' ) THEN
        trialmax = 1.e+4
      ELSEIF ( KEYWORD .EQ. 'Better' ) THEN
        trialmax = 1.e+4
      ELSEIF ( KEYWORD .EQ. 'BETTER' ) THEN
        trialmax = 1.e+4
      ELSEIF ( KEYWORD .EQ. 'good' ) THEN
        trialmax = 1.e+3
      ELSEIF ( KEYWORD .EQ. 'Good' ) THEN
        trialmax = 1.e+3
      ELSEIF ( KEYWORD .EQ. 'GOOD' ) THEN
        trialmax = 1.e+3
      ELSEIF ( KEYWORD .EQ. 'fair' ) THEN
        trialmax = 1.e+2
      ELSEIF ( KEYWORD .EQ. 'Fair' ) THEN
        trialmax = 1.e+2
      ELSEIF ( KEYWORD .EQ. 'FAIR' ) THEN
        trialmax = 1.e+2
      ELSEIF ( KEYWORD .EQ. 'poor' ) THEN
        trialmax = 10.
      ELSEIF ( KEYWORD .EQ. 'Poor' ) THEN
        trialmax = 10.
      ELSEIF ( KEYWORD .EQ. 'POOR' ) THEN
        trialmax = 10.
      ELSE                             ! KEYWORD not recognized
        trialmax = 1.e+2               ! default = 'Fair'
c        KEYWORD = 'Key???'
      ENDIF

      trialmax = ngrains * trialmax   
      IF ( TrialMax .GT. 2.d+09 ) TrialMax = 2.d+09       ! 32 bit machine (INTEGER(4)*4) upper limit 
      
      RETURN
      END FUNCTION
      
