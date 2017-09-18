! ================================================================
! **** SUBROUTINES for PNEAs****
! ================================================================
! delete duplications
SUBROUTINE del_dupl(secev,hsrt,nn,nfin)
  USE fund_const
  USE secular_evolution
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
  USE fund_const
  USE secular_evolution
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
SUBROUTINE evalspline(time,ntot,secev,b2,c2,d2,b3,c3,d3,b4,c4,d4,b5,c5,d5,elem)
  USE secular_evolution
  USE fund_const
  IMPLICIT NONE
  REAL(KIND=dkind),INTENT(IN) :: time
  INTEGER,INTENT(IN) :: ntot
  TYPE(sec_evol),DIMENSION(ntot),INTENT(IN) :: secev
  REAL(KIND=dkind),DIMENSION(ntot),INTENT(IN) :: b2,c2,d2
  REAL(KIND=dkind),DIMENSION(ntot) :: b3,c3,d3
  REAL(KIND=dkind),DIMENSION(ntot) :: b4,c4,d4
  REAL(KIND=dkind),DIMENSION(ntot) :: b5,c5,d5
  REAL(KIND=dkind),DIMENSION(4),INTENT(OUT) :: elem !ecc,inc,Omnod,omega
! end interface
  REAL(KIND=dkind),DIMENSION(ntot) :: t
  INTEGER :: h

! rescaled fundamental interval
  t(1:ntot) = secev(1:ntot)%t - secev(1)%t 

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
  USE fund_const
  IMPLICIT NONE
  REAL(KIND=dkind),INTENT(IN) :: tyr
  REAL(KIND=dkind),INTENT(OUT) :: tmjd
  tmjd = (tyr-1858.87953d0)*365.25d0
END SUBROUTINE yr2mjd

! conversion from MJD to year
SUBROUTINE mjd2yr(tmjd,tyr)
  USE fund_const
  IMPLICIT NONE
  REAL(KIND=dkind),INTENT(IN) :: tmjd
  REAL(KIND=dkind),INTENT(OUT) :: tyr
  tyr = tmjd/365.25d0 + 1858.87953d0
END SUBROUTINE mjd2yr


! serach for mean motion resonances
SUBROUTINE mm_res(a,resflag)
  USE fund_const
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
  USE fund_const
  IMPLICIT NONE
  REAL(KIND=dkind),INTENT(INOUT) :: ang
  ang=mod(ang,360.d0) 
  IF(ang.lt.0.d0)THEN
     ang=ang+360.d0
  ENDIF
END SUBROUTINE in_0_360
