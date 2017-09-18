SUBROUTINE read_pla(ext)
  USE fund_const
  USE orbit_elements
  USE planet_orbits 
  IMPLICIT NONE
! ext='fil' to read filtered planet ephemerides
! ext='dat' to read unfiltered data  
  CHARACTER*3,INTENT(IN) :: ext 
! ------------------- end enterface ------------------------------------ 
  INTEGER :: iun,status,i,j
  INTEGER :: ngf(9,nstepx)
! To read the header
  REAL(KIND=dkind) :: t0  ! is in MODULE planet_orbits
  CHARACTER*60 :: commen
  REAL(KIND=dkind) :: gmv(6),gjyr,gjyr2,g
! File names
  CHARACTER*60 :: elefil
  INTEGER :: le
! **********************************************************************
! planet data                                                       
  INCLUDE 'pldata.h90'
! ***********************************************************************

  evolpla = undefined_orbit_elem
  evolpla%coo='EQU'

!! Open input files
  CALL filnam('.','vpla',ext,elefil,le) 
  CALL filopn(iun,elefil(1:le),'old')
  
  CALL skip(iun,1) ! skipping planet names
!  CALL skip(iun,2) !if present all planets

! Reference time and comment 
  CALL reaflc(iun,'t0',t0,commen) 
  teph0=t0-2400000.5d0
!  write(*,*)'teph0, comments',teph0,commen

  CALL skip(iun,5)  !skip def of reference system and column headers
!  CALL skip(iun,7) !if present all planets

!! Read masses
!  READ(iun,101,end=3,err=3)(gmv(j),j=1,6) 
!101 FORMAT(1p,3d24.16)
!!  write(*,101) gmv
!  DO j=1,6 
!     gmv(j)=1.d0/gmv(j)
!  ENDDO
!!  Gauss constant: from fund_const.mod gk, gms=gk*gk
!!  conversion to internal units: 1AU, 1JYR=365.25 d(Julian year)        
!  gjyr=365.25d0 
!  gjyr2=gjyr*gjyr 
!  g=gms*gjyr2 
!  DO i=1,6 
!     gmv(i)=gmv(i)*g
!  ENDDO
!3 write(9,*)' error in mass input, gm=',gmv 

! Read planetary ephemerides
  DO i=1,50000
     READ(iun,5,END=333) tpla(i) ! difference of time (yr) from teph0
     DO j=inpl,ioupl
        READ(iun,100,END=333) evolpla(j,i)%coord(1:5), &
             & evolpla(j,i)%coord(6),ngf(j,i)
100     FORMAT(f12.9,4f11.7,f11.8,1x,i9)
        evolpla(j,i)%t = teph0 + tpla(i)*365.25d0
     ENDDO 
  ENDDO
5 FORMAT(f15.5)

  write(*,*)'do loop too short!'
  STOP
333 CONTINUE

  CALL filclo(iun,' ')

END SUBROUTINE read_pla
