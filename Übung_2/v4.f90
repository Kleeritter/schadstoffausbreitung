! Programmierpraktikum zur Schadstoffausbreitung in Stadtgebieten
! Übugung 2
PROGRAM ubung2

    IMPLICIT NONE

    INTEGER :: n,tl,dt,sigu,ubalken,wbalken,zq,xq, counter,xgrenz
    REAL    ::sigw,xi,zi,wi,gasdev,ui


    REAL, DIMENSION(:), ALLOCATABLE  :: ges
    REAL, DIMENSION(:), ALLOCATABLE  :: posi
    !Modellparameter
  n= 1000 !Anzahl Partikel
    ubalken = 5 !m/s
    wbalken = 0 !m/s
  zq = 45 !m
  xq = 2000 !m
  counter=0
  xgrenz= 3000 !m

DO WHILE (counter <= n)
    xi=xq
    zi=zq
    ui=ubalken
    wi=wbalken
    print*,xi

    DO WHILE(xi<= xgrenz)
    if (zi<0) THEN
    call position(xi,zi,posi,ui,wi)
    zi=-zi
    wi= -wi
    ELSE
    call position(xi,zi,posi,ui,wi)
    end if
    END DO
    ges=[ges,posi]
    DEALLOCATE( posi)
    counter=counter +1
END DO
open(1, file = 'data1.dat', status = 'new')
write(1,*)ges
close(1)
 contains
    subroutine ragas()
         !REAL gasdev
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
        subroutine position(xi,zi,posi,ui,wi)

    REAL, DIMENSION(:), ALLOCATABLE  :: posi
    REAL, DIMENSION(:), ALLOCATABLE :: transporter
    REAL:: rl,dt,tl,xi,zi,sigw,wi,ui
    INTEGER:: sigu,count
      tl = 100  !s Zeit
     dt = 4 !s Zeitschritt
    sigu= 0 !m/s
    sigw= 0.39 !m/s

    rl= EXP(- dt/tl)
    CALL ragas()
    ui= rl*ui + SQRT((1 - rl**2))*sigu* gasdev
    xi= xi + ui*dt
    print*,xi
    wi= rl*wi + SQRT((1 - rl**2))*sigw* gasdev
    zi= zi + wi*dt
    print*,zi
    transporter=[xi,zi]
        print*,transporter
    posi =[posi,transporter]
    print*,posi
    DEALLOCATE(transporter)
    return

    end subroutine position
END PROGRAM



