! Programmierpraktikum zur Schadstoffausbreitung in Stadtgebieten
! Übugung 2

module berechnen
contains
    subroutine position()
          REAL:: rl,dt,tl
    INTEGER:: ui,sigu,sigw,zi,wi,xi
    rl= EXP(- dt/tl)
    ui= rl*ui + SQRT(1 - rl**2)*sigu* CALL ragas()
    xi= xi + ui*dt
    wi= rl*wi + SQRT(1 - rl**2)*sigw* CALL ragas()
    zi= zi + wi*dt
    print*,xi
    return
    end subroutine position
end module berechnen


module gas
    implicit none
    contains
    subroutine ragas()
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
        end subroutine ragas
end module gas

PROGRAM ubung2

    use berechnen
    use gas
    IMPLICIT NONE

    INTEGER :: n,tl,dt,sigu,ubalken,wbalken,zq,xq, counter,xgrenz
    REAL    ::sigw,xi,zi


    REAL, DIMENSION(:,:,:), ALLOCATABLE  :: c
    
    !Modellparameter
  n= 1000 !Anzahl Partikel
  tl = 100  !s Zeit
  dt = 4 !s Zeitschritt
  sigu= 0 !m/s
  sigw= 0.39 !m/s
  ubalken = 5 !m/s
  wbalken = 0 !m/s
  zq = 45 !m
  xq = 2000 !m
  counter=0
  xgrenz= 3000 !m

DO WHILE (counter <= n)
    xi=xq
    zi=zq
    DO WHILE(xi<= xgrenz)
    if (zi<0) THEN
    call position()
    zi=-zi
    wi= -wi
    ELSE
    call position()
    end if
    END DO
    counter=counter +1
END DO



END PROGRAM



