
PROGRAM test_diff
  USE orbit_elements
  USE fund_const
  IMPLICIT NONE

  TYPE(orbit_elem) :: el1,el2
  REAL(KIND=dkind) :: dist
  REAL(KIND=dkind) :: uakm
  
  uakm=1.4959787d8
!write(*,*)'uakm,radeg',uakm,radeg

! Iridium 33
  el1=undefined_orbit_elem
  el1%coo="KEP"
  el1%t=0.d0
!  el1%coord(1)=7174.6984d0/uakm
!  el1%coord(2)=0.0002288d0
!  el1%coord(3)=86.399d0*radeg
!  el1%coord(4)=121.7d0*radeg
!  el1%coord(5)=90.d0*radeg
  el1%coord(1)=1.1d0
  el1%coord(2)=0.2d0
  el1%coord(3)=85.d0*radeg
  el1%coord(4)=120.d0*radeg
  el1%coord(5)=90.d0*radeg
  el1%coord(6)=0.d0

! Cosmos 2251
  el2=undefined_orbit_elem
  el2%coo="KEP"
  el2%t=0.d0
!  el2%coord(1)=7169.649d0/uakm
!  el2%coord(2)=0.0016027d0
!  el2%coord(3)=74.0355d0*radeg
!  el2%coord(4)=19.5d0*radeg
!  el2%coord(5)=270.d0*radeg
  el2%coord(1)=1.1d0
  el2%coord(2)=0.02d0
  el2%coord(3)=85.d0*radeg
  el2%coord(4)=120.d0*radeg
  el2%coord(5)=90.d0*radeg
  el2%coord(6)=0.d0
  
  CALL orbit_diff(el1,el2,dist)
  WRITE(*,*)'MS_distance = ',dist
  
  CALL SH_orbdist(el1,el2,dist)
  WRITE(*,*)'SH_distance = ',dist


END PROGRAM test_diff

! ***************************************************************
! distance in the trajectory space for bounded Keplerian orbits,
! introduced in Maruskin (2010), CMDA, vol. 108, pp.265-274
! ---------------------------------------------------------------
! TO BE DONE: compute also covariance of this distance
SUBROUTINE orbit_diff(el1,el2,dist)
  USE orbit_elements
  USE fund_const
  IMPLICIT NONE
  TYPE(orbit_elem),INTENT(IN) :: el1,el2
  REAL(KIND=dkind),INTENT(OUT) :: dist
  ! end interface
  TYPE(orbit_elem) :: car1,car2
  INTEGER :: fail_flag
  REAL(KIND=dkind),DIMENSION(3) :: angmom1,lenz1,angmom2,lenz2
  REAL(KIND=dkind) :: vsize,r1,ecc1,beta1,J1
  REAL(KIND=dkind) :: r2,ecc2,beta2,J2
  REAL(KIND=dkind),DIMENSION(3) :: e1,h1,eta1,xi1,e2,h2,eta2,xi2 
  REAL(KIND=dkind) :: Deltapsi
  REAL(KIND=dkind) :: eta1_eta2,xi1_xi2
  REAL(KIND=dkind) :: vel1,energy1,a1,vel2,energy2,a2

  rhs=1
  write(*,*)'gms:',gms

  CALL coo_cha(el1,'CAR',car1,fail_flag)
  IF(fail_flag.ge.5)THEN
     WRITE(*,*)'ERROR in orbit_diff.f90! fail_flag=',fail_flag  
     STOP
  ENDIF
!  write(*,*)'pos1:',car1%coord(1:3)
!  write(*,*)'vel1:',car1%coord(4:6)

  CALL prvec(car1%coord(1:3),gms*car1%coord(4:6),angmom1(1:3))
  CALL prvec(car1%coord(4:6),angmom1(1:3)/gms,lenz1(1:3))
  r1=vsize(car1%coord(1:3))
  e1 = lenz1(1:3) -  car1%coord(1:3)/r1
  ecc1=vsize(e1)
  beta1 = sqrt(1.d0-ecc1**2)
  J1=vsize(angmom1)
  h1=angmom1/J1*beta1

  CALL coo_cha(el2,'CAR',car2,fail_flag)  
  IF(fail_flag.ge.5)THEN
     WRITE(*,*)'ERROR in orbit_diff.f90! fail_flag=',fail_flag  
     STOP
  ENDIF
!  write(*,*)'pos2:',car2%coord(1:3)
!  write(*,*)'vel2:',car2%coord(4:6)

  CALL prvec(car2%coord(1:3),gms*car2%coord(4:6),angmom2(1:3))
  CALL prvec(car2%coord(4:6),angmom2(1:3)/gms,lenz2(1:3))
  r2=vsize(car2%coord(1:3))
  e2 = lenz2(1:3) -  car2%coord(1:3)/r2
  ecc2=vsize(e2)
  beta2 = sqrt(1.d0-ecc2**2)
  J2=vsize(angmom2)
  h2=angmom2/J2*beta2

  eta1= e1 + h1
  xi1 = e1 - h1
  eta2= e2 + h2
  xi2 = e2 - h2 
!  write(*,*)'|eta1|=',vsize(eta1)
!  write(*,*)'|xi1|=',vsize(xi1)
!  write(*,*)'|eta2|=',vsize(eta2)
!  write(*,*)'|xi2|=',vsize(xi2)

  eta1_eta2 = DOT_PRODUCT(eta1,eta2)
  xi1_xi2 = DOT_PRODUCT(xi1,xi2)
  Deltapsi = sqrt((acos(eta1_eta2)**2 + acos(xi1_xi2)**2)/2.d0)

  vel1 = vsize(car1%coord(4:6))
  energy1= 0.5d0*vel1**2 - gms/r1
  a1 = gms/(2.d0*energy1)

  vel2 = vsize(car2%coord(4:6))
  energy2= 0.5d0*vel2**2 - gms/r2
  a2 = gms/(2.d0*energy2)

  dist = sqrt(2.d0*(a1**2+a2**2 - 2.d0*a1*a2*cos(Deltapsi)))
  
END SUBROUTINE orbit_diff

! ====================================================
! Southworth and Hawkings criterion
SUBROUTINE SH_orbdist(el1,el2,dist)
  USE fund_const
  USE orbit_elements
  IMPLICIT NONE
  TYPE(orbit_elem),INTENT(IN) :: el1,el2
  REAL(KIND=dkind),INTENT(OUT) :: dist
  ! end interface
  TYPE(orbit_elem) :: com1,com2
  TYPE(orbit_elem) :: car1,car2
  INTEGER :: fail_flag
  REAL(KIND=dkind) :: q1,e1,inc1,om1,Omnod1
  REAL(KIND=dkind) :: q2,e2,inc2,om2,Omnod2
  REAL(KIND=dkind),DIMENSION(3) :: am1,am2
  REAL(KIND=dkind) :: vsize,pridif
  REAL(KIND=dkind) :: mutI,PI
  REAL(KIND=dkind) :: d1,d2,d3,d4

  CALL coo_cha(el1,'COM',com1,fail_flag)
  IF(fail_flag.ge.5)THEN
     WRITE(*,*)'ERROR in orbit_diff.f90! fail_flag=',fail_flag  
     STOP
  ENDIF
!  write(*,*)'pos1:',com1%coord(1:3)
!  write(*,*)'vel1:',com1%coord(4:6)

  CALL coo_cha(el2,'COM',com2,fail_flag)  
  IF(fail_flag.ge.5)THEN
     WRITE(*,*)'ERROR in orbit_diff.f90! fail_flag=',fail_flag  
     STOP
  ENDIF
!  write(*,*)'pos2:',com2%coord(1:3)
!  write(*,*)'vel2:',com2%coord(4:6)

  q1 =     com1%coord(1)
  e1 =     com1%coord(2)
  inc1 =   com1%coord(3)
  Omnod1 = com1%coord(4)
  om1 =    com1%coord(5)

  q2 =     com2%coord(1)
  e2 =     com2%coord(2)
  inc2 =   com2%coord(3)
  Omnod2 = com2%coord(4)
  om2 =    com2%coord(5)

  d1 = q1-q2
  d2 = e1-e2

!  CALL coo_cha(el1,'CAR',car1,fail_flag)
!  IF(fail_flag.ge.5)THEN
!     WRITE(*,*)'ERROR in orbit_diff.f90! fail_flag=',fail_flag  
!     STOP
!  ENDIF
!  CALL coo_cha(el2,'CAR',car2,fail_flag)
!  IF(fail_flag.ge.5)THEN
!     WRITE(*,*)'ERROR in orbit_diff.f90! fail_flag=',fail_flag  
!     STOP
!  ENDIF
!  CALL prvec(car1%coord(1:3),car1%coord(4:6),am1)
!  am1=am1/vsize(am1)
!  CALL prvec(car2%coord(1:3),car2%coord(4:6),am2)
!  am2=am2/vsize(am2)
!  mutI=acos(DOT_PRODUCT(am1,am2))
!  write(*,*)'mutI=',mutI

! find mutual inclination
  mutI = acos(cos(inc1)*cos(inc2) + sin(inc1)*sin(inc2)*&
       &cos(pridif(Omnod1,Omnod2)))
!  write(*,*)'mutI=',mutI
  d3 = 2.d0*sin(mutI/2.d0)

! angle between perihelion directions
  IF(abs(Omnod1-Omnod2).le.pig) THEN
     PI = (om1-om2) + 2.d0*asin(cos((inc1+inc2)/2.d0)* &
          & sin(pridif(Omnod1,Omnod2)/2.d0)/cos(mutI/2.d0))
  ELSE
     PI = (om1-om2) - 2.d0*asin(cos((inc1+inc2)/2.d0)* &
          & sin(pridif(Omnod1,Omnod2)/2.d0)/cos(mutI/2.d0))
  ENDIF
  d4 = (e1+e2)*sin(PI/2.d0)

!  write(*,*)'d1,d2,d3,d4',d1,d2,d3,d4
  dist = sqrt(d1**2 + d2**2 + d3**2 + d4**2)

END SUBROUTINE SH_orbdist
