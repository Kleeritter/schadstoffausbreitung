 SUBROUTINE SORT(array, num)
 
    IMPLICIT NONE
    
    REAL, INTENT(INOUT), DIMENSION(:) :: array
    REAL  :: temp
    INTEGER, INTENT(IN) :: num
    
    LOGICAL :: swapped = .TRUE.
    
    INTEGER  :: s1, s2
    
    
    DO s1 = num-1, 1, -1
    
       swapped = .FALSE.
       
       DO s2 = 1, s1
       
          IF ( array(s2) > array(s2+1) )   THEN
             temp = array(s2)
             array(s2) = array(s2+1)
             array(s2+1) = temp
             swapped = .TRUE.
          ENDIF
      ENDDO
      
      IF (.NOT. swapped)  THEN
         EXIT
      ENDIF
    
   ENDDO    
 
 
 END SUBROUTINE
