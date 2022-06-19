PROGRAM ubung1

    USE netcdf

    IMPLICIT NONE

    INTEGER :: ca,co,indexcount,xgrenz,counter, nx, ny, nz, i, j,k, z, schritt, ubalken,zq,xq,yq,h,dx,dy,dz,wbalken,n,iz,ux,uz
    REAL    :: xi,zi,fg,fk,gg,gk,dey,dez,Pi,Q,gasdev,winew,zinew,xinew,ui,wi,wiold,ziold,conca,nj,xold,zold


    REAL, DIMENSION(:,:,:), ALLOCATABLE  :: c
    REAL, DIMENSION(:), ALLOCATABLE  :: xlist,zlist, difx,difz
    REAL, DIMENSION(:,:), ALLOCATABLE  :: gitter
    
    !Modellparameter
    nx = 2000
    ny = 2000
    nz=400
    dx = 2
    dy = 2
    dz = 2
    Pi = 3.1415927
    !Randbedingungen
    !Q= 1.5e+5 !540 kg/h also 5.4e+8mg/h und so 150000 
    ubalken= 5
    wbalken=0
    counter=0
    n=1500
    h= 100
    zq=45
    xq=51
    !yq=250
    xgrenz=2000
    indexcount=0
    !Schichtung 
    
    ALLOCATE( c(0:nx,0:ny,0:nz) ) 
    ALLOCATE( difx(0:1000) ) 
    ALLOCATE( difz(0:1000) ) 
    ALLOCATE(gitter(0:int(xgrenz/2 +1),0:int(nz/2 +1)))
    DO i =1,1999,2
        difx(indexcount)=i
        difz(indexcount)=i
        indexcount=indexcount+1
    ENDDO
!print*,SIZE(difx)
!print*,difx
ca=0
allocate( xlist(10**6) )
allocate( zlist(10**6) )
DO WHILE (counter < n)
!print *, counter
    xi=xq
    zi=zq
    ui=ubalken
    wi=wbalken
    !print*, "Gusto 1"
    
    DO WHILE(xi<= xgrenz)
    if (zi<0) THEN
    ziold= zi
    zi=-ziold
    wiold=wi
    wi= -wiold
    CALL positionen(xi,wi,zi)
 
    xi=xinew
    xold=xi
    CALL gitterkarotte(xold)
    !ux=iz
    zi=zinew
    zold=zi
    CALL gitterkarotte(zold)
    !uz=iz
    wi=winew
    !xold=xlist
    !zold=zlist
      CALL counterhe( int(xold),int(zold) )
    !xlist =[xold,xi]
    !zlist =[zold,zi]
    xlist(ca) =xi
    zlist(ca) =zi
    co=ca
    ca=  co +1
    ELSE
        CALL positionen(xi,wi,zi)
        
    xi=xinew
    xold=xi
    CALL gitterkarotte(xold)
    !ux=iz
    zi=zinew
    zold=zi
    CALL gitterkarotte(zold)
    !uz=iz
    wi=winew
    !xold=xlist
    !zold=zlist
      CALL counterhe( int(xold),int(zold) )
    
        !xold=xlist
        !zold=zlist
        
       ! xlist =[xlist,xi]
       ! zlist =[zlist,zi]
    xlist(ca) =xi
    zlist(ca) =zi
    co=ca
    ca=  co +1
       
    end if
    
    END DO
    !ges=[ges,posi]
    !DEALLOCATE( posi)
    counter=counter +1
    
END DO
open(unit=1,file='test.csv',status='unknown')
 write(1,*) xlist
 write(1,*) zlist
!PRINT*, xlist

PRINT*, "Montecarlo fertig"

!CALL gittergurke()
CALL conzentration()

PRINT*, "Konzentrationen fertig"
!CALL positionen(1,1,1)
!print*, zinew

!-- Initialize the NetCDF file
    CALL netcdfout( 'open' )
!
!-- Write concentration array to the NetCDF file
    CALL netcdfout( 'write' )
!
!-- Close the NetCDF file
    CALL netcdfout( 'close' )

    DEALLOCATE ( c )

 CONTAINS


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

 SUBROUTINE positionen(xi,wi,zi)
    IMPLICIT NONE
    REAL:: dt, sigu,sigw,rl,x,zi,wi,xi
    INTEGER::ui, tl
tl = 100  !#s Zeit
dt = 0.4 !# Zeitschritt
sigu= 0 !#m/s
sigw= 0.39! #m/s
rl= EXP(- dt/tl)
ui=5

    CALL ragas()
    xinew= xi + ui*dt
    winew= rl*wi + SQRT((1 - rl**2))*sigw* gasdev
    !print*,gasdev
    zinew= zi + winew*dt
    !print*, zinew
 END SUBROUTINE positionen


 SUBROUTINE koks()
    IMPLICIT NONE
  
 END SUBROUTINE koks
 
SUBROUTINE conzentration()
    IMPLICIT NONE
    !INTEGER:: dx,dz,
    REAL:: q,dx,dz,dt, gitterold
    print*, maxval	(gitter)
    
    dx=2
    dz=2
    q=150
    dt=0.4
    DO i = 0, int(xgrenz/2 +1)
       DO j = 0, int(nz/2 +1)
       	gitterold=gitter(i,j)
   	 gitter(i,j)= q*(gitterold*dt)/(n*dx*dz)
    
  	ENDDO
  ENDDO
  print*, maxval	(gitter)
  !print*,gitter
END SUBROUTINE conzentration
 
 SUBROUTINE counterhe(im,jm)
    IMPLICIT NONE
    INTEGER:: ix,jz,gitterold,im,jm
    !REAL, DIMENSION(:,:), ALLOCATABLE  :: gitter
    !im=int((ix-1)/dx)
    !jm=int((jz-1)/dy)
   ! print*, im, jm
    gitterold= gitter(im,jm)
    gitter(im,jm)=gitterold +1
    !print*,gitter(ix,jz)
 END SUBROUTINE counterhe




 SUBROUTINE gitterkarotte(ix)
    IMPLICIT NONE
    INTEGER ::runde,iz
    REAL::ix
         
         runde= NINT(ix)
 	if( runde >ix ) then
 		if( MODULO(runde, 2)==0 ) then
 		ix= runde-1
 		else
 		ix= runde
 		end if
   
   else if( runde <ix ) then
 		if( MODULO(runde, 2)==0 ) then
 		ix= runde+1
 		else
 		ix= runde
 		end if
  
   else
 		if( MODULO(runde, 2)==0 ) then
 		ix= runde-1
 		else
 		ix= runde
 		end if
      
   end if
        
    !print*, iz
 END SUBROUTINE gitterkarotte

SUBROUTINE gaus()
    IMPLICIT NONE
    !REAL:
    !INTEGER: Q,ubalken,nx,nz

    !Q= 150 #!540 kg/h also 5.4e+8mg/h und so 150000 
    !ubalken= 5

    !for i,indexei in zip(range(1 , nx),(range(0,nx))):
        !for j,indexej in zip(range(0,nz),(range(0,nz))):
            !#dez = gg*((dx*i)**gk)
       	    !#dez = 2*sigw**2*tl*((i/ubalken)- ubalken +ubalken*math.exp(- (i/(ubalken*tl))))
            !dez= 2*sigw**2*tl*((i/ubalken)-tl +tl*math.exp(-i/(ubalken*tl)))
            !cdkack.append(Q/(math.sqrt(2* math.pi)*math.sqrt(dez)*ubalken) *(math.exp((-((j*dx)-zq)**2)/(2*dez)) +math.exp((-((j*dz)+zq)**2)/(2*dez))))
            !#cdc.append(co)
    print*, "Alla"
END SUBROUTINE gaus
 SUBROUTINE netcdfout( action )

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT (IN) ::  action

    
    INTEGER, SAVE ::  id_dim_x, id_dim_y,id_dim_z, id_set, &
                      id_var_c, id_var_x, id_var_y,id_var_z,nc_stat

    REAL, DIMENSION(:), ALLOCATABLE ::  netcdf_data
	INTEGER:: indexcountnet	
!
!-- Initialize the NetCDF file before the first output is made
    IF ( action == 'open' )  THEN

!--    Deleting output file
       OPEN(1, FILE='scheibe.nc')
       CLOSE(1, STATUS='DELETE')
    
!
!--    Open the netcdf file for writing
       nc_stat = NF90_CREATE( 'scheibe.nc', NF90_NOCLOBBER, id_set )
       IF ( nc_stat /= NF90_NOERR )  PRINT*, '+++ netcdf error'

!
!--    Define some global attributes of the dataset
       nc_stat = NF90_PUT_ATT( id_set, NF90_GLOBAL, 'Conventions', 'COARDS' )
       nc_stat = NF90_PUT_ATT( id_set, NF90_GLOBAL, 'title', &
                               'Bruder' )

!
!--    Define the spatial dimensions and coordinates
!--    Define x + y coordinate
       nc_stat = NF90_DEF_DIM( id_set, 'x', nx+1, id_dim_x )
       nc_stat = NF90_DEF_VAR( id_set, 'x', NF90_DOUBLE, id_dim_x, id_var_x )
       nc_stat = NF90_PUT_ATT( id_set, id_var_x, 'units', 'meters' )

       
       nc_stat = NF90_DEF_DIM( id_set, 'y', ny+1, id_dim_y )
       nc_stat = NF90_DEF_VAR( id_set, 'y', NF90_DOUBLE, id_dim_y, id_var_y )
       nc_stat = NF90_PUT_ATT( id_set, id_var_y, 'units', 'meters' )
       

       
!
!--    Define the output array(s)
       nc_stat = NF90_DEF_VAR( id_set, 'c', NF90_DOUBLE,            &
                               (/ id_dim_x, id_dim_y /), id_var_c )
       nc_stat = NF90_PUT_ATT( id_set, id_var_c, 'units', 'gm^-3' )

!
!--    Leave NetCDF define mode
       nc_stat = NF90_ENDDEF( id_set )

!
!--    Write data for x and z axis
       ALLOCATE( netcdf_data(0:nx) )
       indexcountnet=0
       DO  i = 1, nx
          netcdf_data(i) = i
          indexcountnet=indexcountnet+1
       ENDDO
	
       nc_stat = NF90_PUT_VAR( id_set, id_var_x, netcdf_data, &
                               start = (/ 1 /), count = (/ nx+1 /) )
       DEALLOCATE( netcdf_data )

       ALLOCATE( netcdf_data(0:ny) )
       indexcountnet=0
       DO  i = 1, ny
          netcdf_data(i) = i
          indexcountnet=indexcountnet+1
       ENDDO

       nc_stat = NF90_PUT_VAR( id_set, id_var_y, netcdf_data, &
                               start = (/ 1 /), count = (/ ny+1 /) )
       DEALLOCATE( netcdf_data)
       !Z
     
       RETURN

    ENDIF

    IF ( action == 'close' )  THEN
!
!--    Close the NetCDF file
       nc_stat = NF90_CLOSE( id_set )

       RETURN

    ENDIF
    
    
    IF ( action == 'write' )  THEN

  print*, "goose"
!--    Write the array data
       nc_stat = NF90_PUT_VAR( id_set, id_var_c, gitter(0:nx,0:ny), &
                               start = (/ 1, 1/), &              
                               count = (/ nx+1, ny+1 /) )
    ENDIF

 END SUBROUTINE netcdfout

 END PROGRAM
