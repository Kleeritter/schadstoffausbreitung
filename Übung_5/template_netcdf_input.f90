PROGRAM netcdfin

   IMPLICIT NONE

   ...

   CALL netcdfin ( 'openread', 'u', u )

   ...


END PROGRAM 


 SUBROUTINE netcdfin( action, field, array )

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT (IN) ::  action, field

    
    INTEGER, SAVE ::  id_set, &
                      id_var, nc_stat

    REAL, INTENT (INOUT), DIMENSION(:,:) :: array

    
    IF ( action == 'openread' )  THEN
   
       nc_stat = NF90_OPEN( 'input.nc', NF90_NOWRITE, id_set )
       IF ( nc_stat /= NF90_NOERR )  PRINT*, '+++ netcdf error'
       
       nc_stat = NF90_INQ_VARID( id_set, field, id_var )
       IF ( nc_stat /= NF90_NOERR )  PRINT*, '+++ netcdf error'

       nc_stat = NF90_GET_VAR( id_set, id_var, array ) 
       IF ( nc_stat /= NF90_NOERR )  PRINT*, '+++ netcdf error'

       nc_stat = NF90_CLOSE( id_set )
       IF ( nc_stat /= NF90_NOERR )  PRINT*, '+++ netcdf error'

    ENDIF
    
    RETURN
END SUBROUTINE
