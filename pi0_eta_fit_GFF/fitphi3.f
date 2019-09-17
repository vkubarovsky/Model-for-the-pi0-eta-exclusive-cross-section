      REAL FUNCTION FITPHI3(X)
      COMMON /PAWPAR/ P(10)
      FITPHI3=P(1)+P(2)*COS(2*X)+p(3)*2*cos(X)
      RETURN
      END
