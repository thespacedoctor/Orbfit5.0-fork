PROGRAM vastdat2fla
  USE fund_const
!  USE planet_orbits
  USE orbit_elements
  IMPLICIT NONE

  INTEGER,PARAMETER :: nstepx=1000000
  TYPE(orbit_elem) :: evolast(nstepx)
  DOUBLE PRECISION :: tast(nstepx),teph0
  INTEGER :: iun,iunout
  INTEGER :: i,ifin
  REAL(KIND=dkind) :: t0 
  CHARACTER*60 :: commen

! **********************************************************************

  CALL filopn(iunout,'vastdat.fla','unknown')

  evolast = undefined_orbit_elem
  evolast%coo='EQU'
  CALL filopn(iun,'vast.dat','old')
  READ(iun,*)
! Reference time and comment 
  CALL reaflc(iun,'t0',t0,commen) 
  teph0=t0-2400000.5d0
  READ(iun,*)
  READ(iun,*)
! Read asteroid ephemerides
  DO i=1,50000
     READ(iun,5,END=333) tast(i) ! difference of time (yr) from teph0
     READ(iun,101,END=333) evolast(i)%coord(1:5), &
          & evolast(i)%coord(6)
!100     FORMAT(f12.9,4f11.7,f11.8,1x,i9,1p,e12.4)
101  FORMAT(f12.9,4f11.7,f11.8)
     evolast(i)%t = teph0 + tast(i)*365.25d0
  ENDDO
5 FORMAT(f15.5)

  write(*,*)'do loop too short!'
  STOP
333 CONTINUE
  CALL filclo(iun,' ')
  ifin=i-1

  DO i=1,ifin
     WRITE(iunout,100) tast(i),evolast(i)%coord(1:5), &
          & evolast(i)%coord(6)
  ENDDO
100 FORMAT(f15.5,2x,f12.9,1x,4(f11.7,1x),f11.8)  

  CALL filclo(iunout,' ')
  
END PROGRAM vastdat2fla


