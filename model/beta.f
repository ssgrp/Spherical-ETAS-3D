
C
C       ==========================================
C       Purpose: Compute the density of beta distribution
C       Input :  p  --- Parameter  ( p > 0 )
C                q  --- Parameter  ( q > 0 )
C       Output:   --- x^{p-1}*(1-x)^{q-1}/B(p,q)
C       Routine called: GAMMA for computing â(x)
C       ==========================================
C

        real*8 FUNCTION dbeta(x, p, q)
        implicit real*8 (a-h, o-z)

        dbeta = x**(p-1d0)*(1d0-x)**(q-1d0)/beta(p,q)

        return 
        end

        real*8 FUNCTION BETA(P,Q)
C
C       ==========================================
C       Purpose: Compute the beta function B(p,q)
C       Input :  p  --- Parameter  ( p > 0 )
C                q  --- Parameter  ( q > 0 )
C       Output:  BT --- B(p,q)
C       Routine called: GAMMA for computing â(x)
C       ==========================================
C
        IMPLICIT real*8 (A-H,O-Z)
        beta=exp(QLGAMMA(P)+QLGAMMA(Q)-QLGAMMA(P+Q))

        RETURN
        END


C=======================================================================
C
 
