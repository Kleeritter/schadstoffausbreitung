! Programmierpraktikum zur Schadstoffausbreitung in Stadtgebieten
! Test zur Prüfung der Kompilierumgebung
!
! Compilation + running + checking:
! ifort -I /muksoft/packages/netcdf4_hdf5parallel/4411c_443f/hdf5-1.10.0-patch1/mvapich2-2.3rc1/intel/2018.1.163/include/ -c test.f90
! ifort -o test test.o -L/muksoft/packages/netcdf4_hdf5parallel/4411c_443f/hdf5-1.10.0-patch1/mvapich2-2.3rc1/intel/2018.1.163/lib64/ -lnetcdf -lnetcdff
! ./test
! ncview test_result.nc
 
PROGRAM test

    USE netcdf

    IMPLICIT NONE

    INTEGER :: nx, ny, i, j 
    REAL    :: dx, dy


    REAL, DIMENSION(:,:), ALLOCATABLE  :: c
    
    
    nx = 10
    ny = 30
    dx = 2.0
    dy = 1.0
    
    ALLOCATE( c(0:nx,0:ny) ) 

    DO i = 0, nx
       DO j = 0, ny
       
          c(i,j) = i*dx + j*dy
       
       ENDDO
    ENDDO
    


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

 SUBROUTINE netcdfout( action )

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT (IN) ::  action

    
    INTEGER, SAVE ::  id_dim_x, id_dim_y, id_set, &
                      id_var_c, id_var_x, id_var_y, nc_stat

    REAL, DIMENSION(:), ALLOCATABLE ::  netcdf_data

!
!-- Initialize the NetCDF file before the first output is made
    IF ( action == 'open' )  THEN

!--    Deleting output file
       OPEN(1, FILE='test_result.nc')
       CLOSE(1, STATUS='DELETE')
    
!
!--    Open the netcdf file for writing
       nc_stat = NF90_CREATE( 'test_result.nc', NF90_NOCLOBBER, id_set )
       IF ( nc_stat /= NF90_NOERR )  PRINT*, '+++ netcdf error'

!
!--    Define some global attributes of the dataset
       nc_stat = NF90_PUT_ATT( id_set, NF90_GLOBAL, 'Conventions', 'COARDS' )
       nc_stat = NF90_PUT_ATT( id_set, NF90_GLOBAL, 'title', &
                               'UDM test' )

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
       nc_stat = NF90_PUT_ATT( id_set, id_var_c, 'units', 'kg' )

!
!--    Leave NetCDF define mode
       nc_stat = NF90_ENDDEF( id_set )

!
!--    Write data for x and z axis
       ALLOCATE( netcdf_data(0:nx) )
       DO  i = 0, nx
          netcdf_data(i) = i * dx
       ENDDO

       nc_stat = NF90_PUT_VAR( id_set, id_var_x, netcdf_data, &
                               start = (/ 1 /), count = (/ nx+1 /) )
       DEALLOCATE( netcdf_data )

       ALLOCATE( netcdf_data(0:ny) )
       DO  i = 0, ny
          netcdf_data(i) = i * dy
       ENDDO

       nc_stat = NF90_PUT_VAR( id_set, id_var_y, netcdf_data, &
                               start = (/ 1 /), count = (/ ny+1 /) )
       DEALLOCATE( netcdf_data )


       RETURN

    ENDIF

    IF ( action == 'close' )  THEN
!
!--    Close the NetCDF file
       nc_stat = NF90_CLOSE( id_set )

       RETURN

    ENDIF
    
    
    IF ( action == 'write' )  THEN

!
!--    Write the array data
       nc_stat = NF90_PUT_VAR( id_set, id_var_c, c(0:nx,0:ny), &
                               start = (/ 1, 1 /), &              
                               count = (/ nx+1, ny+1 /) )
    ENDIF

 END SUBROUTINE netcdfout

 END PROGRAM
