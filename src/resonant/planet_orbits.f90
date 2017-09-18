
MODULE planet_orbits

  USE orbit_elements
  IMPLICIT NONE
  PRIVATE
  
  TYPE(orbit_elem),DIMENSION(9) :: el_pla
  INTEGER,PARAMETER :: nstepx=1000000
  TYPE(orbit_elem) :: evolpla(9,nstepx)
  DOUBLE PRECISION :: tpla(nstepx),teph0
  LOGICAL,PARAMETER :: force_circ=.false. ! tests planet on circular orbits
!  LOGICAL,PARAMETER :: force_circ=.true. ! tests planet on circular orbits

  ! variables  
  PUBLIC :: el_pla,evolpla,tpla,nstepx,teph0,force_circ
  ! routines
  PUBLIC :: placar,planet_elems

CONTAINS
! ======================================================================
! PLACAR - get planet Cartesian coordinates (ecliptic J2000)             
! ======================================================================
SUBROUTINE placar(n,t0,xpl,ifla) 
  USE fund_const
  USE planet_masses
! input: planet num., epoch time, 
! flag for getting Planet (heliocentric; ifla=1) or Sun (barycentric; ifla=2)
  INTEGER,INTENT(IN) :: n   
  DOUBLE PRECISION, INTENT(IN) :: t0 
  INTEGER, INTENT(IN) :: ifla
! output: heliocentric state vector of Planet (equinoctal, ecliptic)     
!      or barycentric state vector of Sun (equinoctal, ecliptic)        
  DOUBLE PRECISION, INTENT(OUT) :: xpl(6) 
! =============JPL EPHEM===============                                 
! data for masses                                                      
  INCLUDE 'jplhdr.h90' 
! output of JPL routine, Julian date, rotation matrix                   
  double precision et(2),rrd(6) 
! integers for call to JPl routines                                     
  integer ntarg,ncent,istate
  DOUBLE PRECISION ap(nplax)
! ====================================
! circular planet case
     IF(force_circ)THEN
        write(*,*)'placar: ERROR! circular case not implemented'
        STOP
     ENDIF
!  write(*,*)'Planet on elliptical keplerian orbits'
! JPL Earth vector at observation time                                  
  et(1)=2400000.5d0 
  et(2)=t0 
  if (ifla.eq.1) then 
     ntarg=n 
     ncent=11 
  else 
     ntarg=11 
     ncent=12 
  endif
! duplicate computation of gmse, in case masjpl has not been called yet 
  gmse=gms*(1.d0+cval(8+ntarg)/cval(18)) 
! ****** added on Sat Jun 14 1997 ******                                
! first istate need to be=2  (dpleph calculates also vel.)              
  istate=2 
  call dpleph(et,ntarg,ncent,rrd,istate)
! Change of reference system EQUM00 ---> ECLM00                         
!      call rotpn(rot,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0) 

  xpl(1:3)=MATMUL(roteqec,rrd(1:3)) 
  xpl(4:6)=MATMUL(roteqec,rrd(4:6))

END SUBROUTINE placar

! ======================================================================
! PLANET_ELEMS - interpolate planet coordinates from output of orbit9
! ======================================================================
SUBROUTINE planet_elems(npla,time,elpla)
  USE fund_const,ONLY: dkind,dpig
  USE planet_masses, ONLY: nplax
  INTEGER,INTENT(IN) :: npla
  REAL(KIND=dkind),INTENT(IN) :: time !interpolation time
  REAL(KIND=dkind),DIMENSION(6),INTENT(OUT) :: elpla !EQU elems of planet
! ------------------------- END INTERFACE -------------------------
  REAL(KIND=dkind) :: t1 !first extremum of the interpolation interval
  REAL(KIND=dkind) :: lambda !coeff of convex combination
  REAL(KIND=dkind) :: intstep ! planet integration time step (in yr)
  INTEGER :: nstep ! no. of time steps to reach t1
  INTEGER :: j,i !loop indexes
  TYPE(orbit_elem),DIMENSION(2) :: el
  REAL(KIND=dkind) :: ap(nplax)
  DOUBLE PRECISION :: princ,gaussk
! -----------------------------------------------------------------

! circular planet case
   IF(force_circ)THEN
      CALL pla_dat(ap)
      elpla(1)=ap(npla)
      elpla(2:5)=0.d0
! lambda (= elle in this case, since omega=Omnod=0)
!  Gauss constant ( internal units: 1AU, 1JYR=365.25 d(Julian year) )                                                      
      gaussk=0.01720209895d0*365.25d0 
      elpla(6)=evolpla(npla,1)%coord(6) + (time-teph0)* elpla(1)**(-1.5d0)*gaussk
      elpla(6) = MOD(elpla(6),dpig)
!      write(*,*)'teph0,time,elpla(6)',teph0,time,elpla(6)
!      stop

      RETURN
   ENDIF
! intstep is not a parameter, allowing to change its value 
   intstep= tpla(2)-tpla(1)
! *** assuming that time >= teph0 ***
  nstep= FLOOR((time-(teph0+tpla(1)*365.25d0))/(365.25d0*intstep))

!  *** t0+tpla(nstep+1) <= time <= t0+tpla(nstep+2) ***
  t1 = teph0+tpla(nstep+1)*365.25d0
  lambda=(time-t1)/(365.25d0*intstep) ! interpolating parameter
  IF(lambda.lt.-1.d-10.OR.lambda.gt.1.d0+1.d-10)THEN
     WRITE(*,*) 'error! lambda = ',lambda
     WRITE(*,*) 'time = ',time,'nstep = ',nstep
     STOP
  ENDIF
  DO i=1,2
     el(i) = evolpla(npla,nstep+i) ! HINT: coo='EQU'
  ENDDO
! interpolate EQU coordinates (1 to 6)
  elpla(1:6)= (1.d0-lambda)*el(1)%coord(1:6) + &
       & lambda*el(2)%coord(1:6) 
  
END SUBROUTINE planet_elems

END MODULE planet_orbits

