SUBROUTINE pla_dat(ap)
!  USE planet_masses
  IMPLICIT NONE
  DOUBLE PRECISION ap(9)
! ============ planet semimajor axis ==================                 
  ap(1)= .3870992058d0 
  ap(2)= .7233274811d0 
  ap(3)= 1.0000036214d0 
  ap(4)= 1.5235973464d0 
  ap(5)= 5.2024107723d0 
  ap(6)= 9.5575876779d0 
  ap(7)= 19.3008879212d0 
  ap(8)= 30.2722024706d0 
  ap(9)= 39.7533710065d0 
  RETURN 
END SUBROUTINE pla_dat
