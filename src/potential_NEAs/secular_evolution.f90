MODULE secular_evolution
  USE fund_const
  IMPLICIT NONE
  PRIVATE
! types
  PUBLIC :: sec_evol,undefined_sec_evol

! routines
  PUBLIC :: sec_interp,comp_interp,evalspline,yr2mjd,mjd2yr,in_0_360
! ---------------------------------------
  TYPE sec_evol
     REAL(KIND=dkind) :: t
     REAL(KIND=dkind) :: a,om,nod,ecc,inc
     INTEGER :: cf   ! crossing flag, read from pro files  
                     ! cf=0 when crossings omega=k*pi/2
  END TYPE sec_evol
! ---------------------------------------
! undefined secular evolution
  TYPE(sec_evol), PARAMETER :: undefined_sec_evol = SEC_EVOL( &
&  0.d0,           & ! time (MJD)
&  0.d0,           & ! a    (au)
&  0.d0,           & ! omega (deg)
&  0.d0,           & ! Omnod (deg)
&  0.d0,           & ! eccentricity
&  0.d0,           & ! inclination (deg)
&  55              & ! default crossing flag
&     ) 

CONTAINS

  SUBROUTINE sec_interp(name,deltat,t0,a,ntot,secev,omega1,omega2,deltanod, &
       & Omnod2,b2,c2,d2,b3,c3,d3,b4,c4,d4,b5,c5,d5)
    IMPLICIT NONE
    INTEGER, PARAMETER :: hmax=10000 
    CHARACTER(LEN=9),INTENT(IN) :: name ! asteroid name 
    REAL(KIND=dkind),INTENT(OUT) :: deltat ! length of fund. interval J (yr)
    REAL(KIND=dkind),INTENT(OUT) :: t0,a ! init time (MJD), semimajor axis
    REAL(KIND=dkind),INTENT(OUT) :: omega1,omega2 ! extrema of omega in J
    REAL(KIND=dkind),INTENT(OUT) :: deltanod ! excursion of Omega in deltat yrs
    REAL(KIND=dkind),INTENT(OUT) :: Omnod2 ! last value of Omega in deltat yrs
! coefficients of cubic splines
    REAL(KIND=dkind),DIMENSION(hmax),INTENT(OUT) :: b2,c2,d2 !ecc
    REAL(KIND=dkind),DIMENSION(hmax),INTENT(OUT) :: b3,c3,d3 !inc
    REAL(KIND=dkind),DIMENSION(hmax),INTENT(OUT) :: b4,c4,d4 !nod
    REAL(KIND=dkind),DIMENSION(hmax),INTENT(OUT) :: b5,c5,d5 !om
    INTEGER,INTENT(OUT) :: ntot ! number of nodes in J
! end interface
    CHARACTER(LEN=9) :: astnam  ! asteroid name
    INTEGER :: lnam            
    REAL(KIND=dkind) :: Rbar,Z  ! auxiliary
! secular evolution data 
    TYPE(sec_evol),DIMENSION(hmax),INTENT(OUT) :: secev !(angles in deg)
    INTEGER,DIMENSION(hmax) :: hsrt ! index of sorted data
    INTEGER :: nfin             ! auxiliary
! time in fundamental interval J=[t1,t2] is rescaled:
! t_rescal(1) = 0; t_rescal(ntot) = t2-t1
    REAL(KIND=dkind),DIMENSION(hmax) :: t_rescal
    REAL(KIND=dkind) :: Omnod1 ! value at the extrema of J
    INTEGER :: h,nn             ! loop indexes

    secev=undefined_sec_evol

    astnam=name
     CALL rmsp(astnam,lnam) 

! *** mev files *** 
     OPEN (10,file='mev_files/'//astnam(1:lnam)//'.mev',status='old')
     READ(10,100) t0,a,Rbar,Z 
     secev(1:hmax)%a = a
     DO h=1,hmax
        READ(10,101,END=33)secev(h)%t,secev(h)%om,secev(h)%nod,secev(h)%ecc, &
             & secev(h)%inc  ! angles in degrees
        secev(h)%cf=99  ! dummy value, to remove duplicates 
        secev(h)%t=t0 + secev(h)%t*365.25d0
!        secev(h)%om=mod(secev(h)%om,360.d0)   !! NO MODULUS FOR omega !!  
!        secev(h)%nod=mod(secev(h)%nod,360.d0) !! NO MODULUS FOR Omnod !!  
     ENDDO
     WRITE(*,*)'hmax too small:',hmax
     STOP
33   CONTINUE
     CLOSE(10)
     nn=h
100  FORMAT(f13.6,1x,f18.15,1x,E22.15,1x,E22.15)   
101  FORMAT(f10.2,1x,f11.5,1x,f11.5,1x,f10.7,1x,f13.7,18x,i3)  
102  FORMAT(f10.2,1x,f11.5,1x,f11.5,1x,f10.7,1x,f13.7,i3)

! *** pro files ***
     OPEN (10,file='pro_files/'//astnam(1:lnam)//'.pro',status='old')
     deltat=0.d0
     DO h=nn,hmax
        READ(10,102,END=34)secev(h)%t,secev(h)%om,secev(h)%nod,secev(h)%ecc, &
             & secev(h)%inc,secev(h)%cf ! angles in degrees
        IF(secev(h)%cf.eq.0)THEN
           deltat=secev(h)%t-deltat
        END IF
        secev(h)%t=t0 + secev(h)%t*365.25d0  ! evolution times in MJD
!        secev(h)%om=mod(secev(h)%om,360.d0)    !! NO MODULUS FOR omega !!  
!        secev(h)%nod=mod(secev(h)%nod,360.d0) !! NO MODULUS FOR Omnod !!  
     ENDDO
     WRITE(*,*)'hmax too small:',hmax
     STOP
34   CONTINUE
     CLOSE(10)
     nn=h-1
     write(*,*)'length of fundamental interval (yr)',deltat

! sorting secev w.r.t time
     CALL heapsort(secev(1:nn)%t,nn,hsrt)

! delete duplicates
     CALL del_dupl(secev(1:nn),hsrt,nn,nfin)
!     DO h=1,nfin
!        write(*,*)'*',secev(h)%nod,secev(h)%om,secev(h)%t
!     ENDDO

! search for fundamental interval
     CALL fund_int(secev(1:nfin),nfin,ntot)
! rounding extremal values of omega
     secev(1)%om =    real(NINT(secev(1)%om))
     secev(ntot)%om = real(NINT(secev(ntot)%om))
!     write(*,*)'fundamental interval'
!     write(*,('(3(a12,10x))')) 'time (MJD)','Omnod','omega'  
!     DO h=1,ntot
!        write(*,*) secev(h)%t,secev(h)%nod,secev(h)%om
!     ENDDO
!     read(*,*)

! rescaling time of fundamental interval to [0,deltat]
     t_rescal(1:ntot)=secev(1:ntot)%t-secev(1)%t
     write(*,*)'t_trasl=',secev(1)%t,'t2-t1=',secev(ntot)%t-secev(1)%t

! rounding values of omega,Omnod at the extrema of the fundamental interval:
! these values must be in [0,360] 
!     omega1=mod(real(NINT(secev(1)%om)),360.d0) ! non viene usato!!
!     omega2=mod(real(NINT(secev(ntot)%om)),360.d0)
!
     omega1=secev(1)%om 
     omega2=secev(ntot)%om
     write(*,*)'omega1,omega2',omega1,omega2

     Omnod1=secev(1)%nod; Omnod2=secev(ntot)%nod
     write(*,*)'Omega1,Omega2',Omnod1,Omnod2
     deltanod = Omnod2-Omnod1
     write(*,*)'deltanod',deltanod

! HINT! *** actually I could be computed from a,e,Z ***

! spline interpolation in fundamental interval J
     CALL spline (ntot,t_rescal(1:ntot),secev(1:ntot)%ecc,b2(1:ntot), &
          & c2(1:ntot),d2(1:ntot)) 
     CALL spline (ntot,t_rescal(1:ntot),secev(1:ntot)%inc,b3(1:ntot), &
          & c3(1:ntot),d3(1:ntot)) 
     CALL spline (ntot,t_rescal(1:ntot),secev(1:ntot)%nod,b4(1:ntot), &
          & c4(1:ntot),d4(1:ntot)) 
     CALL spline (ntot,t_rescal(1:ntot),secev(1:ntot)%om,b5(1:ntot), &
          & c5(1:ntot),d5(1:ntot)) 

  END SUBROUTINE sec_interp

! ========================================================================
! compute orbital elements at a given time by spline interpolation, 
! using the data in the fundamental interval J
! written by G.F. Gronchi 2012
SUBROUTINE comp_interp(time,dtMJD,circ,hmax,ntot, &
     & secevin,omega2,deltanod,Omnod2,b2,c2,d2,b3,c3,d3,b4,c4,d4,b5,c5,d5,elem)
  USE fund_const
  IMPLICIT NONE
  REAL(KIND=dkind),INTENT(IN) :: time ! orbital elements time (yr)
  REAL(KIND=dkind),INTENT(IN) :: dtMJD ! length of fund. int. (d)
  LOGICAL,INTENT(IN) :: circ !.true. for circulation/.false. for libration
  INTEGER,INTENT(IN) :: hmax,ntot
  TYPE(sec_evol),DIMENSION(hmax),INTENT(IN) :: secevin
  REAL(KIND=dkind),INTENT(IN) :: omega2 ! 2nd extremum of omega in J
  REAL(KIND=dkind),INTENT(IN) :: deltanod ! excursion of Omega in deltat yrs
  REAL(KIND=dkind),INTENT(IN) :: Omnod2 ! second value of Omega in deltat yrs
  REAL(KIND=dkind),DIMENSION(hmax),INTENT(IN) :: b2,c2,d2
  REAL(KIND=dkind),DIMENSION(hmax),INTENT(IN) :: b3,c3,d3
  REAL(KIND=dkind),DIMENSION(hmax),INTENT(IN) :: b4,c4,d4
  REAL(KIND=dkind),DIMENSION(hmax),INTENT(IN) :: b5,c5,d5
  REAL(KIND=dkind),DIMENSION(4),INTENT(OUT) :: elem !ecc,inc,Omnod,omega
! end interface
  TYPE(sec_evol),DIMENSION(hmax) :: secev
  INTEGER :: h !loop index
  REAL(KIND=dkind) :: tMJD ! orbital elements time (d)
  REAL(KIND=dkind) :: tau,t_resid 
  INTEGER :: nquad ! number of quadrant: from 0 to 3 for circulators
                   !                     from 0 to ?? for librators
  LOGICAL :: fwdshift,bwdshift

  fwdshift=.false. !initialization
  bwdshift=.false. !initialization
  secev=secevin

  CALL yr2mjd(time,tMJD)
  tMJD = tMJD-secev(1)%t !rescaling orgin of time, for evalspline 
!  write(*,*)'tMJD,dtMJD',tMJD,dtMJD

  IF(circ)THEN
     ! tMJD must be a value between 0 and 4*dtMJD, otherwise we shift it
     IF(tMJD.lt.0.d0)THEN
        WRITE(*,*)'requested time in the past: shifting foreward!'
        WRITE(*,*)tMJD,tMJD+4.d0*dtMJD
        tMJD=tMJD+4.d0*dtMJD !use the known evolution in the future
        fwdshift=.true.
     ELSEIF(tMJD.gt.4.d0*dtMJD)THEN
        WRITE(*,*)'requested time in the future: shifting backward!'
        WRITE(*,*)tMJD,tMJD-4.d0*dtMJD
        tMJD=tMJD-4.d0*dtMJD !use the known evolution in the past
        bwdshift=.true.
     ENDIF
     IF(tMJD.lt.0.d0)THEN
        WRITE(*,*)'requested time again in the past (should not happen!)'
        STOP
     ELSEIF(tMJD.gt.4.d0*dtMJD)THEN
        WRITE(*,*)'requested time again in the future (should not happen!)'
        STOP
     ENDIF
  ELSE
     ! tMJD must be a value between 0 and 2*dtMJD, otherwise we shift it
     IF(tMJD.lt.0.d0)THEN
        WRITE(*,*)'requested time in the past: shifting foreward!'
        WRITE(*,*)tMJD,tMJD+2.d0*dtMJD
        tMJD=tMJD+2.d0*dtMJD !use the known evolution in the future
        fwdshift=.true.
     ELSEIF(tMJD.gt.2.d0*dtMJD)THEN
        WRITE(*,*)'requested time in the future: shifting backward!'
        WRITE(*,*)tMJD,tMJD-2.d0*dtMJD
        tMJD=tMJD-2.d0*dtMJD !use the known evolution in the past
        bwdshift=.true.
     ENDIF
     IF(tMJD.lt.0.d0)THEN
        WRITE(*,*)'requested time again in the past (should not happen!)'
        STOP
     ELSEIF(tMJD.gt.2.d0*dtMJD)THEN
        WRITE(*,*)'requested time again in the future (should not happen!)'
        STOP
     ENDIF
  ENDIF

  tau = mod(tMJD,dtMJD)
  CALL quadrante(circ,tMJD,dtMJD,t_resid,nquad)
!  write(*,('(a18,1x,3(es12.5,1x),i1)'))'tMJD,tau,t_resid,nquad',tMJD,tau,t_resid,nquad
!  read(*,*)


  IF (circ)THEN
     IF(mod(nquad,2).eq.0)THEN
        CALL evalspline(tau,ntot,secev(1:ntot),b2(1:ntot),c2(1:ntot),&
             & d2(1:ntot),b3(1:ntot),c3(1:ntot),d3(1:ntot), &
             & b4(1:ntot),c4(1:ntot),d4(1:ntot), &
             & b5(1:ntot),c5(1:ntot),d5(1:ntot),elem(1:4))
        IF(nquad.eq.0)THEN
           CALL in_0_360(elem(4)) 
              ! do nothing
        ELSEIF(nquad.eq.2)THEN
           elem(3)= elem(3) + 2.d0*deltanod
           elem(4)=mod(180.d0 + elem(4),360.d0) !assuming direct circulation
        ENDIF
     ELSEIF(mod(nquad,2).eq.1)THEN
        CALL evalspline(dtMJD-tau,ntot,secev(1:ntot), &
             & b2(1:ntot),c2(1:ntot),&
             & d2(1:ntot),b3(1:ntot),c3(1:ntot),d3(1:ntot), &
             & b4(1:ntot),c4(1:ntot),d4(1:ntot), &
             & b5(1:ntot),c5(1:ntot),d5(1:ntot),elem(1:4))
        IF(nquad.eq.1)THEN
!           elem(3)= 2.d0*(deltanod/90.d0)*secev(ntot)%om - elem(3)
           elem(3)= Omnod2 - (elem(3)-Omnod2)
           elem(4)=mod(2.d0*omega2 - elem(4),360.d0)
       ELSEIF(nquad.eq.3)THEN
!           elem(3)= 2.d0*(deltanod/90.d0)*secev(ntot)%om &
!                & - elem(3) + 2.d0*deltanod
           elem(3)= Omnod2 - (elem(3)-Omnod2) + 2.d0*deltanod
!assuming direct circulation
           elem(4)=mod(180.d0 + 2.d0*omega2 - elem(4),360.d0)
        ENDIF
     ENDIF
     IF(fwdshift)THEN
        elem(3) = elem(3)-4.d0*deltanod
     ENDIF
     IF(bwdshift)THEN
        elem(3) = elem(3)+4.d0*deltanod
     ENDIF
  ENDIF

  IF(.not.circ)THEN
     IF(mod(nquad,2).eq.0)THEN
        CALL evalspline(tau,ntot,secev(1:ntot),b2(1:ntot),c2(1:ntot),&
             & d2(1:ntot),b3(1:ntot),c3(1:ntot),d3(1:ntot), &
             & b4(1:ntot),c4(1:ntot),d4(1:ntot), &
             & b5(1:ntot),c5(1:ntot),d5(1:ntot),elem(1:4))
        CALL in_0_360(elem(4)) 
     ELSEIF(mod(nquad,2).eq.1)THEN
        CALL evalspline(dtMJD-tau,ntot,secev(1:ntot), &
             & b2(1:ntot),c2(1:ntot),&
             & d2(1:ntot),b3(1:ntot),c3(1:ntot),d3(1:ntot), &
             & b4(1:ntot),c4(1:ntot),d4(1:ntot), &
             & b5(1:ntot),c5(1:ntot),d5(1:ntot),elem(1:4))
        elem(3)= Omnod2 - (elem(3)-Omnod2)
        elem(4)=mod(2.d0*omega2 - elem(4),360.d0)
     ENDIF
    IF(fwdshift)THEN
        elem(3) = elem(3)-2.d0*deltanod
     ENDIF
     IF(bwdshift)THEN
        elem(3) = elem(3)+2.d0*deltanod
     ENDIF
  ENDIF
!     write(*,'(5(f15.5,2x))')time,elem(1:4)

END SUBROUTINE comp_interp

! ================================================================
! **** SUBROUTINES for PNEAs****
! ================================================================
! delete duplications
SUBROUTINE del_dupl(secev,hsrt,nn,nfin)
!  USE fund_const
!  USE secular_evolution
  IMPLICIT NONE
  TYPE(sec_evol),DIMENSION(nn),INTENT(INOUT) :: secev
  INTEGER,INTENT(IN),DIMENSION(nn) :: hsrt
  INTEGER,INTENT(IN) :: nn
  INTEGER,INTENT(OUT) :: nfin
! end interface
  TYPE(sec_evol),DIMENSION(nn) :: sectmp
  INTEGER :: htmp,h

  htmp=1
  sectmp(1)%t = secev(hsrt(1))%t
  sectmp(1)%om = secev(hsrt(1))%om
  sectmp(1)%nod = secev(hsrt(1))%nod
  sectmp(1)%ecc = secev(hsrt(1))%ecc
  sectmp(1)%inc = secev(hsrt(1))%inc
  sectmp(1)%cf = secev(hsrt(1))%cf

  DO h=2,nn
     IF(abs(secev(hsrt(h))%t-secev(hsrt(h-1))%t).le.1.d-3)THEN
        IF(secev(hsrt(h))%cf.eq.99)THEN
           ! skip line h
        ELSE
           sectmp(htmp)%t = secev(hsrt(h))%t
           sectmp(htmp)%om = secev(hsrt(h))%om
           sectmp(htmp)%nod = secev(hsrt(h))%nod
           sectmp(htmp)%ecc = secev(hsrt(h))%ecc
           sectmp(htmp)%inc = secev(hsrt(h))%inc
           sectmp(htmp)%cf = secev(hsrt(h))%cf
        ENDIF
     ELSE
        htmp=htmp+1
        sectmp(htmp)%t = secev(hsrt(h))%t
        sectmp(htmp)%om = secev(hsrt(h))%om
        sectmp(htmp)%nod = secev(hsrt(h))%nod
        sectmp(htmp)%ecc = secev(hsrt(h))%ecc
        sectmp(htmp)%inc = secev(hsrt(h))%inc
        sectmp(htmp)%cf = secev(hsrt(h))%cf
     ENDIF
  ENDDO
  nfin=htmp 
  DO h=1,nfin
!     secev(h)=undefined_sec_evol
     secev(h)%t = sectmp(h)%t 
     secev(h)%om = sectmp(h)%om 
     secev(h)%nod = sectmp(h)%nod 
     secev(h)%ecc = sectmp(h)%ecc 
     secev(h)%inc = sectmp(h)%inc 
     secev(h)%cf = sectmp(h)%cf 
  ENDDO

END SUBROUTINE del_dupl

! select fundamental interval
SUBROUTINE fund_int(secev,nfin,ntot)
!  USE fund_const
!  USE secular_evolution
  IMPLICIT NONE
  TYPE(sec_evol),DIMENSION(nfin),INTENT(INOUT) :: secev
  INTEGER,INTENT(IN) :: nfin
  INTEGER,INTENT(OUT) :: ntot
! end interface
  INTEGER :: htmp,h
  LOGICAL :: found1,found2

  htmp=0
  found1=.false.;found2=.false.;
  DO h=1,nfin
     IF(secev(h)%cf.eq.0)THEN
        IF(found1)THEN
           found2=.true.
        ENDIF
     ENDIF
     IF(secev(h)%cf.eq.0)THEN
        found1=.true.
     ENDIF
     IF(found1)THEN
        htmp=htmp+1
        secev(htmp)%t = secev(h)%t
        secev(htmp)%om = secev(h)%om
        secev(htmp)%nod = secev(h)%nod
        secev(htmp)%ecc = secev(h)%ecc
        secev(htmp)%inc = secev(h)%inc
        secev(htmp)%cf = secev(h)%cf
     ENDIF
     ntot=htmp 
     IF(found2) RETURN
  ENDDO

END SUBROUTINE fund_int

! compute averaged orbital elements by spline interpolation
SUBROUTINE evalspline(timein,ntot,secev,b2,c2,d2,b3,c3,d3,b4,c4,d4,b5,c5,d5,elem)
!  USE secular_evolution
!  USE fund_const
  IMPLICIT NONE
  REAL(KIND=dkind),INTENT(IN) :: timein
  INTEGER,INTENT(IN) :: ntot
  TYPE(sec_evol),DIMENSION(ntot),INTENT(IN) :: secev
  REAL(KIND=dkind),DIMENSION(ntot),INTENT(IN) :: b2,c2,d2
  REAL(KIND=dkind),DIMENSION(ntot),INTENT(IN) :: b3,c3,d3
  REAL(KIND=dkind),DIMENSION(ntot),INTENT(IN) :: b4,c4,d4
  REAL(KIND=dkind),DIMENSION(ntot),INTENT(IN) :: b5,c5,d5
  REAL(KIND=dkind),DIMENSION(4),INTENT(OUT) :: elem !ecc,inc,Omnod,omega
! end interface
  REAL(KIND=dkind),DIMENSION(ntot) :: t
  REAL(KIND=dkind) :: time,averdt
  INTEGER :: h

  time=timein

! rescaled fundamental interval
  t(1:ntot) = secev(1:ntot)%t - secev(1)%t 

  averdt=t(ntot)/ntot

! check if time is within the fundamental interval
  IF(time.lt.t(1).or.time.gt.t(ntot))THEN
     WRITE(*,*) 'evalspline WARNING: &
          & outside the fundamental interval!',t(1),time,t(ntot)
     WRITE(*,*) 'differences:',time-t(1),t(ntot)-time
     IF(abs(time-t(1)).lt.1.d-5*averdT)THEN
!        time = t(1)    !rounding
     ELSEIF(abs(t(ntot)-time).lt.1.d-5*averdT)THEN
!        time = t(ntot) !rounding
     ELSE
        WRITE(*,*) 'evalspline: FATAL ERROR!!'
        STOP
     ENDIF
  ENDIF

  DO h=1,ntot-1
     IF(time.ge.t(h).and.time.le.t(h+1))THEN
!        write(*,*)'time,th,th+1',time,t(h),t(h+1)
        elem(1) = secev(h)%ecc + b2(h)*(time-t(h)) + &
             & c2(h)*(time-t(h))**2 + d2(h)*(time-t(h))**3        
        elem(2) = secev(h)%inc + b3(h)*(time-t(h)) + &
             & c3(h)*(time-t(h))**2 + d3(h)*(time-t(h))**3        
        elem(3) = secev(h)%nod + b4(h)*(time-t(h)) + &
             & c4(h)*(time-t(h))**2 + d4(h)*(time-t(h))**3        
        elem(4) = secev(h)%om  + b5(h)*(time-t(h)) + &
             & c5(h)*(time-t(h))**2 + d5(h)*(time-t(h))**3        
     ENDIF
  ENDDO

END SUBROUTINE evalspline

! conversion from year to MJD
SUBROUTINE yr2mjd(tyr,tmjd)
!  USE fund_const
  IMPLICIT NONE
  REAL(KIND=dkind),INTENT(IN) :: tyr
  REAL(KIND=dkind),INTENT(OUT) :: tmjd
  tmjd = (tyr-1858.87953d0)*365.25d0
END SUBROUTINE yr2mjd

! conversion from MJD to year
SUBROUTINE mjd2yr(tmjd,tyr)
!  USE fund_const
  IMPLICIT NONE
  REAL(KIND=dkind),INTENT(IN) :: tmjd
  REAL(KIND=dkind),INTENT(OUT) :: tyr
  tyr = tmjd/365.25d0 + 1858.87953d0
END SUBROUTINE mjd2yr


! serach for mean motion resonances
SUBROUTINE mm_res(a,resflag)
!  USE fund_const
  IMPLICIT NONE
  REAL(KIND=dkind) :: a
  LOGICAL :: resflag
  ! end interface
  REAL(KIND=dkind),DIMENSION(8) :: ap
  INTEGER :: h,k

  resflag=.false.

  ap(1)= .3870992058d0;
  ap(2)= 0.7233274811d0; 
  ap(3)= 1.0000036214d0; 
  ap(4)= 1.5235973464d0; 
  ap(5)= 5.2024107723d0; 
  ap(6)= 9.5575876779d0; 
  ap(7)= 19.3008879212d0; 
  ap(8)= 30.2722024706d0; 

  DO h=1,5
     DO k=1,5
        IF(abs((a/ap(5))**1.5d0-dble(h)/dble(k)).lt.1.d-2)THEN
           WRITE(*,*) 'resonant case',h,k,a,ap(3)
           resflag=.true.
        ENDIF
     ENDDO
  ENDDO
  
END SUBROUTINE mm_res

SUBROUTINE in_0_360(ang)
!  USE fund_const
  IMPLICIT NONE
  REAL(KIND=dkind),INTENT(INOUT) :: ang
  ang=mod(ang,360.d0) 
  IF(ang.lt.0.d0)THEN
     ang=ang+360.d0
  ENDIF
END SUBROUTINE in_0_360

END MODULE secular_evolution

! =====================================
SUBROUTINE quadrante(circ,t,dt,resid,nquad)
  implicit none
  logical,intent(in) :: circ
  double precision,intent(in) :: t,dt
  integer,intent(out) :: nquad
  double precision,intent(out) :: resid
! END INTERFACE
  double precision :: ratio

  IF(dt.le.0.d0)THEN
     write(*,*)'non-positive value of dt:',dt
     STOP
  ENDIF

! computing residual of tmjd mod dt
  resid = mod(t,dt)
  IF(t.ge.0.d0)THEN
!     IF(circ)THEN
!        resid=resid
!     ELSE
!        resid=resid
!     ENDIF
  ELSEIF(t.lt.0.d0)THEN
     IF(circ)THEN
        resid = 4.d0*dt-abs(resid)
     ELSE
        resid = 2.d0*dt-abs(resid)
     ENDIF
  ENDIF

! computing quadrant
  nquad = NINT((t-resid)/dt)
  IF(t.ge.0.d0)THEN
     IF(circ)THEN
        nquad = mod(nquad,4)
     ELSE
        nquad = mod(nquad,2)
     ENDIF
  ELSEIF(t.lt.0.d0)THEN
     IF(circ)THEN
        nquad = 4 - mod(-nquad,4)
     ELSE
        nquad = 2 - mod(-nquad,2)
     ENDIF
  ENDIF
!  write(*,*)'nquad=', nquad

END SUBROUTINE quadrante
