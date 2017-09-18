MODULE planetorb_circ
  USE fund_const
  USE orbit_elements
  IMPLICIT NONE 
  PRIVATE
  TYPE(orbit_elem),DIMENSION(8) :: elpl
  PUBLIC :: elpl,placirc
CONTAINS
  
  SUBROUTINE placirc
    INTEGER :: i,j
    elpl(1:8)= undefined_orbit_elem
    DO i=1,8
       DO j=1,6
          elpl(i)%coord(j)=0.d0
       ENDDO
! **** this options causes problems in the sign of the moid ****
!       elpl(i)%coo='COM'
!       elpl(i)%t=55000.d0 !default
! ***************************************8
       elpl(i)%coo='KEP'
       elpl(i)%t=0.d0 !default
    ENDDO
! ============ planet semimajor axes ==================                 
    elpl(1)%coord(1)= 0.3870992058d0   ! Mercury
    elpl(2)%coord(1)= 0.7233274811d0   ! Venus
    elpl(3)%coord(1)= 1.0000036214d0   ! Earth
    elpl(4)%coord(1)= 1.5235973464d0   ! Mars
    elpl(5)%coord(1)= 5.2024107723d0   ! Jupiter
    elpl(6)%coord(1)= 9.5575876779d0   ! Saturn
    elpl(7)%coord(1)= 19.3008879212d0  ! Uranus
    elpl(8)%coord(1)= 30.2722024706d0  ! Neptune
!  elpl(9)%coord(1)= 39.7533710065d0 ! Pluto (dwarf-planet)
  
  END SUBROUTINE placirc

END MODULE planetorb_circ
