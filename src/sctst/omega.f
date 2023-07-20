        Double precision FUNCTION omega(n,NN,FI,XI)
        Integer n, NN(n)
        Double precision FI, XI(n)
        Integer I
        Double precision tp
c
c         n         = number of orthogonal coupled DOF
c         NN        = old vector of quantum numbers
c         FI        = magnitude of imaginary frequency (REAL)
c         XI        = off-diagonal anharmonicities (cm-1) (i.e. X_kF) coupling rxn coord to other vibs.
c         omega     = imaginary frequency, corrected for coupling with orthogonal DOF (Eq. 7c in Miller et al. 1990)

        tp=0.0d0
        DO I=1, n
          tp = tp + XI(I) * ( NN(I)+0.5d0 ) 
        ENDDO
        omega = FI - tp

        RETURN
        END function
