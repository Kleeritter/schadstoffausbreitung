module gas
!-- The function gasdev() returns a normality pseudorandom value
FUNCTION GASDEV()
      REAL gasdev
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1,zv1,zv2
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       CALL RANDOM_NUMBER( zv1 )
        CALL RANDOM_NUMBER( zv2 )
        v1 = 2.0 * zv1 - 1.0
        v2 = 2.0 * zv2 - 1.0
!        v1=2.*ran1(idum)-1.
!        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
END 
END module gas
