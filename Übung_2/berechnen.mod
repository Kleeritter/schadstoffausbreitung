module berechnen
FUNCTION position()
    REAL: rl
    INTEGER:
    rl= EXP(- dt/tl)
    ui= rl*ui + SQRT(1- rl**2)*sigu* CALL GASDEV()
    xi= xi + ui*dt
    wi= rl*wi + SQRT(1- rl**2)*sigw* CALL GASDEV()
    zi= zi + wi*dt
    print*,xi
    return
END function position
END module berechnen
