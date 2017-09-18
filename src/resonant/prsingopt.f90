! ******************************************************************
!     read options and initial conditions to call rhsloop.f90             
! ******************************************************************
SUBROUTINE prsingopt(elem,omfreq,hmax) 
  USE orbit_elements
  USE fund_const
  IMPLICIT NONE 
  TYPE(orbit_elem),INTENT(IN) :: elem ! asteroid orbital element
  REAL(KIND=dkind),INTENT(OUT) :: omfreq,hmax 
! ----------- end interface ----------------
! rhsloop output                                                  
  REAL(KIND=dkind) :: ddd(4),eee(4)
  INTEGER nnn(4) 
  INCLUDE 'pldata.h90' ! planet data
  CALL rhsloop(elem,ddd,eee,nnn) 
  omfreq=-ddd(1) 
  hmax = abs(pig/(omfreq*100)) 
! ========= MAXIMUM INITIAL STEP ===================================
  IF(hmax.gt.10) THEN 
!  IF(hmax.gt.500) THEN 
     WRITE(*,*)'MAXIMUM STEP EXCEEDED!' 
     hmax = 10 
!     hmax=500
  ENDIF
  write(*,*)'STEP =',hmax 
  RETURN 
END SUBROUTINE prsingopt
