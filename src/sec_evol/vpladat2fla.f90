PROGRAM vpladat2fla
  USE fund_const
!  USE planet_orbits
  USE orbit_elements
  IMPLICIT NONE

  INTEGER,PARAMETER :: nstepx=1000000
  TYPE(orbit_elem) :: evolpla(9,nstepx)
  DOUBLE PRECISION :: tpla(nstepx),teph0
  INTEGER :: iun,iun2,iun3,iun4,iun5,iun6
  INTEGER :: i,j,ifin
  CHARACTER*3 :: ext !extension file planets secular elements
                     !'vpla.dat' or 'vpla.fil'

  REAL(KIND=dkind) :: t0 
  CHARACTER*60 :: commen

! **********************************************************************

  CALL filopn(iun2,'venus.fla','unknown')
  CALL filopn(iun3,'earth.fla','unknown')
  CALL filopn(iun4,'mars.fla','unknown')
  CALL filopn(iun5,'jupiter.fla','unknown')
  CALL filopn(iun6,'saturn.fla','unknown')


  evolpla = undefined_orbit_elem
  evolpla%coo='EQU'
  CALL filopn(iun,'vpla.dat','old')
  CALL skip(iun,1) ! skipping planet names
! Reference time and comment 
  CALL reaflc(iun,'t0',t0,commen) 
  teph0=t0-2400000.5d0
  CALL skip(iun,5)  !skip def of reference system and column headers
! Read planetary ephemerides
  DO i=1,50000
     READ(iun,5,END=333) tpla(i) ! difference of time (yr) from teph0
     DO j=2,6
        READ(iun,101,END=333) evolpla(j,i)%coord(1:5), &
             & evolpla(j,i)%coord(6)
!100     FORMAT(f12.9,4f11.7,f11.8,1x,i9,1p,e12.4)
101     FORMAT(f12.9,4f11.7,f11.8)
        evolpla(j,i)%t = teph0 + tpla(i)*365.25d0
     ENDDO 
  ENDDO
5 FORMAT(f15.5)

  write(*,*)'do loop too short!'
  STOP
333 CONTINUE
  CALL filclo(iun,' ')
  ifin=i-1

  DO i=1,ifin
     WRITE(iun2,100) tpla(i),evolpla(2,i)%coord(1:5), &
          & evolpla(2,i)%coord(6)
     WRITE(iun3,100) tpla(i),evolpla(3,i)%coord(1:5), &
          & evolpla(3,i)%coord(6)
     WRITE(iun4,100) tpla(i),evolpla(4,i)%coord(1:5), &
          & evolpla(4,i)%coord(6)
     WRITE(iun5,100) tpla(i),evolpla(5,i)%coord(1:5), &
          & evolpla(5,i)%coord(6)
     WRITE(iun6,100) tpla(i),evolpla(6,i)%coord(1:5), &
          & evolpla(6,i)%coord(6)
100  FORMAT(f15.5,2x,f12.9,1x,4(f11.7,1x),f11.8)
  ENDDO
  
  CALL filclo(iun2,' ')
  CALL filclo(iun3,' ')
  CALL filclo(iun4,' ')
  CALL filclo(iun5,' ')
  CALL filclo(iun6,' ')
  
END PROGRAM vpladat2fla


