      REAL FUNCTION AQ(X)
      COMMON /PAWPAR/PAR(3)

      AQ=PAR(1)*X**PAR(2)*(1.-X)**PAR(3)

      RETURN
      END

      
