      SUBROUTINE PUTPAR
      COMMON /PAWPAR/F(10)
      VECTOR P(10)
      DO I=1,10
         F(I)=P(I)
      ENDDO
      RETURN
      END
