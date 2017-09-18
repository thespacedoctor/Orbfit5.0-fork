
PROGRAM propag_PNEAs
  USE fund_const
  USE output_control
  USE orbit_elements
  USE secular_evolution
  USE critical_points
  IMPLICIT NONE
!  CHARACTER*9 :: name
  REAL(KIND=dkind) :: time,tmjd 
  CHARACTER*9 :: astname
  TYPE(orbit_elem) :: elkep1,elkep2
  REAL(KIND=dkind) :: t0,a,Rbar,Z 
  INTEGER, PARAMETER :: hmax=10000
  REAL(KIND=dkind),DIMENSION(hmax) :: t,omega,Omnod,ecc,inc
!  REAL(KIND=dkind),DIMENSION(hmax) :: t_srt,omega_srt,Omnod_srt,ecc_srt,inc_srt
!  REAL(KIND=dkind),DIMENSION(hmax) :: t_tmp,omega_tmp,Omnod_tmp,ecc_tmp,inc_tmp
!  REAL(KIND=dkind) :: ttmp,omtmp,Omntmp,etmp,itmp
  INTEGER :: cfltmp
  REAL(KIND=dkind),DIMENSION(hmax) :: b2,c2,d2
  REAL(KIND=dkind),DIMENSION(hmax) :: b3,c3,d3
  REAL(KIND=dkind),DIMENSION(hmax) :: b4,c4,d4
  REAL(KIND=dkind),DIMENSION(hmax) :: b5,c5,d5
  REAL(KIND=dkind) :: deltat,dtMJD ! length of fundamental interval (yr,d)
  REAL(KIND=dkind),DIMENSION(hmax) :: t_rescal
  INTEGER,DIMENSION(hmax) :: cflag,cflag_srt,cflag_tmp
  CHARACTER*9 :: name
  INTEGER :: lnam,nn,htmp,nfin,ntot
  INTEGER,DIMENSION(hmax) :: hsrt
  REAL(KIND=dkind),DIMENSION(hmax) :: q,incl
  INTEGER :: h,j,n ! loop indexes
  TYPE(sec_evol),DIMENSION(hmax) :: secev
  REAL(KIND=dkind) :: tau !
  INTEGER :: nstep 
  INTEGER :: iunout,iunpneas 
!  REAL(KIND=dkind),DIMENSION(4) :: elem !omega,Omnod,ecc,inc
  REAL(KIND=dkind),DIMENSION(4,hmax) :: elem !omega,Omnod,ecc,inc
  REAL(KIND=dkind) :: omega1,omega2
  LOGICAL :: circ ! circ=.true.  for circulators
                  ! circ=.false. for librators

  REAL(KIND=dkind),DIMENSION(nminx) :: dmintil ! local minima with sign
  INTEGER :: nummin
  LOGICAL :: resflag

  REAL(KIND=8) :: qmin,qmax,h_mag
  CHARACTER*8 :: type !nf = numbered, found in neolist
                      !nn = numbered, not in neolist
                      !mf = multiopp, found in neolist
                      !mn = multiopp, not in neolist
  CHARACTER*1 :: lc
  INTEGER,PARAMETER :: nastx=100000


! initialization
  elem = 0.d0 
!  circ=.false.
  q = 0.d0
  rhs=1

  elkep1=undefined_orbit_elem
  elkep1%coo='KEP'
  elkep1%coord=0.d0
  elkep1%coord(1)=1.0000036214d0 
!  elkep1%coord(1)=1.5235973464  ! Mars

  elkep2=undefined_orbit_elem
  elkep2%coo='KEP'
  elkep2%coord=0.d0

! *** inserire il do loop sui PNEAs ***
! *** leggere per ogni asteroide il flag C/L ***

  CALL filopn(iunpneas,'PNEAs','old')
106 FORMAT(1x,a9,2(3x,f10.5),2x,f8.2,3x,a8,3x,a1)

  DO n=1,nastx
     READ(iunpneas,106,END=111) astname,qmin,qmax,h_mag,type,lc
     IF (lc.eq.'C')THEN
        circ=.true.
     ELSEIF (lc.eq.'L')THEN
        circ=.false.
     ENDIF

!  astname='1951'
!  circ=.true.
!  astname='242450'
!  astname='105140'
!  circ=.false.

     secev = undefined_sec_evol
!  t=0.d0;omega=0.d0;Omnod=0.d0;ecc=0.d0;inc=0.d0 
  ! read file mev
     name=astname
     CALL rmsp(name,lnam) 
!  write(*,*)'mev_files/'//name(1:lnam)//'.mev'
     OPEN (10,file='mev_files/'//name(1:lnam)//'.mev',status='old')
     READ(10,100) t0,a,Rbar,Z 
     secev(1:hmax)%a = a
     
     GOTO 343 ! *** skip resonance check ***
  ! check whether it is resonant
     CALL mm_res(a,resflag)
     IF(resflag)THEN
        STOP
     ENDIF
343  CONTINUE

     DO h=1,hmax
! angles in degrees
        READ(10,101,END=33)secev(h)%t,secev(h)%om,secev(h)%nod,secev(h)%ecc, &
             & secev(h)%inc !,cflag(h)
        secev(h)%cf=99 ! dummy value, to remove duplicates 
        secev(h)%t=t0 + secev(h)%t*365.25d0
        secev(h)%om=mod(secev(h)%om,360.d0)
        secev(h)%nod=mod(secev(h)%nod,360.d0)
     ENDDO
     WRITE(*,*)'hmax too small:',hmax
     STOP
33   CONTINUE
     CLOSE(10)
     nn=h-1
     
100  FORMAT(f13.6,1x,f18.15,1x,E22.15,1x,E22.15)   
101  FORMAT(f10.2,1x,f11.5,1x,f11.5,1x,f10.7,1x,f13.7,18x,i3)  
102  FORMAT(f10.2,1x,f11.5,1x,f11.5,1x,f10.7,1x,f13.7,i3)
     
     OPEN (10,file='pro_files/'//name(1:lnam)//'.pro',status='old')
     deltat=0.d0
     DO h=nn+1,hmax
! angles in degrees
        READ(10,102,END=34)secev(h)%t,secev(h)%om,secev(h)%nod,secev(h)%ecc, &
             & secev(h)%inc,secev(h)%cf
        IF(secev(h)%cf.eq.0)THEN
           deltat=secev(h)%t-deltat
           dtMJD=deltat*365.25d0
        END IF
        ! converting into MJD
        secev(h)%t=t0 + secev(h)%t*365.25d0
        secev(h)%om=mod(secev(h)%om,360.d0)
        secev(h)%nod=mod(secev(h)%nod,360.d0)
     ENDDO
     WRITE(*,*)'hmax too small:',hmax
     STOP
34   CONTINUE
     CLOSE(10)
     nn=h-1
     write(*,*)'length of fundamental interval (yr,MJD)'
     write(*,'(2(f20.5,2x))') deltat,dtMJD

! sorting raws w.r.t time
     CALL heapsort(secev(1:nn)%t,nn,hsrt)
     
     ! delete duplicates
     CALL del_dupl(secev(1:nn),hsrt,nn,nfin)
!  DO h=1,nfin
!     write(*,101)secev(h)%t/365.25+1858.87953d0,secev(h)%om, &
!          & secev(h)%nod,secev(h)%ecc,secev(h)%inc,secev(h)%cf
!  ENDDO

! search for fundamental interval
     CALL fund_int(secev(1:nfin),nfin,ntot)
     write(*,*)'fundamental interval'

! check extremum values of omega
     IF(secev(1)%om.gt.359.9d0)THEN
        secev(1)%om=secev(1)%om-360.d0 
     ENDIF
     IF(secev(ntot)%om.lt.1d0)THEN
        secev(ntot)%om=secev(ntot)%om+360.d0 
     ENDIF

!  DO h=1,ntot
!     write(*,101)secev(h)%t/365.25+1858.87953d0,secev(h)%om, &
!          & secev(h)%nod,secev(h)%ecc,secev(h)%inc,secev(h)%cf
!  ENDDO
     t_rescal(1:ntot)=secev(1:ntot)%t-secev(1)%t
     write(*,*)'t_trasl=',secev(1)%t
!  write(*,*) 'rescaled times (MJD) in fundamental interval:', &
!       & t_rescal(1:ntot)

! values of omega at the extrema of the fundamental interval
     omega1=secev(1)%om
     omega2=secev(ntot)%om
     write(*,*)'omega1,omega2',omega1,omega2


! HINT! *** actually I could be computed from a,e,Z ***

! spline interpolation in fundamental interval
     CALL spline (ntot,t_rescal(1:ntot),secev(1:ntot)%ecc,b2(1:ntot), &
          & c2(1:ntot),d2(1:ntot)) 
     CALL spline (ntot,t_rescal(1:ntot),secev(1:ntot)%inc,b3(1:ntot), &
          & c3(1:ntot),d3(1:ntot)) 
     CALL spline (ntot,t_rescal(1:ntot),secev(1:ntot)%nod,b4(1:ntot), &
          & c4(1:ntot),d4(1:ntot)) 
     CALL spline (ntot,t_rescal(1:ntot),secev(1:ntot)%om,b5(1:ntot), &
          & c5(1:ntot),d5(1:ntot)) 
     
     CALL filopn(iunout,'AST'//name(1:lnam)//'.evol','unknown')
     write(iunout,108) name(1:lnam),a,Rbar,Z,0,0,0,0,0
108  FORMAT(a9,2x,f10.6,2x,es15.5,2x,f10.6,5(2x,i1))          

     DO h=1,2000
        time = 2010+(h-1)*100
        CALL yr2mjd(time,tmjd)
!     write(*,*)'tmjd:',tmjd
        tmjd = tmjd-secev(1)%t !rescaling

        IF(tmjd.lt.0.d0)THEN
!        WRITE(*,*)'requested time in the past'
           tmjd=tmjd+4.d0*dtMJD
        ENDIF
        tau=mod(tmjd,dtMJD)
        nstep = (tmjd-tau)/dtMJD
!     write(*,*)'tmjd,dtMJD,tau,nstep',tmjd,dtMJD,tau,nstep
!     nstep = FLOOR(tmjd/(2d0*dtMJD))
!     tau = tmjd-nstep*2.d0*dtMJD 
!    write(*,*)'tau,nstep',tau,nstep  

        IF (circ)THEN
           IF(mod(nstep,2).eq.0)THEN
              CALL evalspline(tau,ntot,secev(1:ntot),b2(1:ntot),c2(1:ntot),&
                   & d2(1:ntot),b3(1:ntot),c3(1:ntot),d3(1:ntot), &
                   & b4(1:ntot),c4(1:ntot),d4(1:ntot), &
                   & b5(1:ntot),c5(1:ntot),d5(1:ntot),elem(1:4,h))
              CALL in_0_360(elem(4,h)) 
              IF(mod(nstep,4).eq.0)THEN
              ! do nothing
              ! write(*,*)'0',elem(4,h)
              ELSEIF(mod(nstep,4).eq.2)THEN
!              write(*,*)'2',elem(4,h)
                 elem(4,h)=mod(180.d0 + elem(4,h),360.d0)
              ENDIF
           ELSEIF(mod(nstep,2).eq.1)THEN
              CALL evalspline(dtMJD-tau,ntot,secev(1:ntot), &
                   & b2(1:ntot),c2(1:ntot),&
                   & d2(1:ntot),b3(1:ntot),c3(1:ntot),d3(1:ntot), &
                   & b4(1:ntot),c4(1:ntot),d4(1:ntot), &
                   & b5(1:ntot),c5(1:ntot),d5(1:ntot),elem(1:4,h))
              CALL in_0_360(elem(4,h)) 
               IF(mod(nstep,4).eq.1)THEN
!              write(*,*)'1',elem(4,h)
                 elem(4,h)=mod(2.d0*omega2 - elem(4,h),360.d0)
              ELSEIF(mod(nstep,4).eq.3)THEN
!              write(*,*)'3',elem(4,h)
                 elem(4,h)=mod(180.d0 + 2.d0*omega2 - elem(4,h),360.d0)
              ENDIF
           ENDIF
        ENDIF
        
        IF(.not.circ)THEN
           IF(mod(nstep,2).eq.0)THEN
              CALL evalspline(tau,ntot,secev(1:ntot),b2(1:ntot),c2(1:ntot),&
                   & d2(1:ntot),b3(1:ntot),c3(1:ntot),d3(1:ntot), &
                   & b4(1:ntot),c4(1:ntot),d4(1:ntot), &
                   & b5(1:ntot),c5(1:ntot),d5(1:ntot),elem(1:4,h))
           ELSEIF(mod(nstep,2).eq.1)THEN
              CALL evalspline(dtMJD-tau,ntot,secev(1:ntot), &
                   & b2(1:ntot),c2(1:ntot),&
                   & d2(1:ntot),b3(1:ntot),c3(1:ntot),d3(1:ntot), &
                   & b4(1:ntot),c4(1:ntot),d4(1:ntot), &
                   & b5(1:ntot),c5(1:ntot),d5(1:ntot),elem(1:4,h))
              elem(4,h)=mod(2.d0*omega2 - elem(4,h),360.d0)
           ENDIF
           CALL in_0_360(elem(4,h)) 
        ENDIF
!     write(*,'(5(f15.5,2x))')time,elem(1:4,h)

! check whether it is NEA at current epoch

        q(h) = a*(1.d0-elem(1,h))
!     write(*,*)'*', Z/sqrt(gms*(365.25d0**2)*a*(1.d0-elem(1,h)**2))
        incl(h) = acos(Z/sqrt(gms*(365.25d0**2)*a*(1.d0-elem(1,h)**2))) 
        elkep2%coord(1)=a
        elkep2%coord(2)=elem(1,h)
        elkep2%coord(3)=incl(h) 
        elkep2%coord(5)=elem(4,h)*radeg
        dmintil=0.d0 ! initialization
        CALL dmintil_rms(elkep1,elkep2,nummin,dmintil)

!     IF(q(h).lt.1.3d0)THEN
        IF(nummin.eq.1)THEN
           write(iunout,313) time,q(h),elem(1,h),incl(h)*degrad,elem(4,h), &
                & nummin,dmintil(1),0.d0,0.d0
        ELSEIF(nummin.eq.2)THEN
           write(iunout,313) time,q(h),elem(1,h),incl(h)*degrad,elem(4,h), &
                & nummin,dmintil(1:2),0.d0
        ELSEIF(nummin.eq.3)THEN
           write(iunout,313) time,q(h),elem(1,h),incl(h)*degrad,elem(4,h), &
                & nummin,dmintil(1:3)
        ELSE
           WRITE(*,*)'ERROR! nummin=',nummin
           STOP
        ENDIF
        ! ENDIF
     ENDDO

313  FORMAT(f10.2,1x,4(f13.5,1x),i4,1x,3(f13.5,1x))

!     IF(time.lt.t(1).or.time.ge.t(nn))THEN
!        write(*,*)'outside time boundary'
!        STOP
!     ENDIF
     
     CALL filclo(iunout,' ')
     
  ENDDO
  WRITE(*,*)'too small nastx:',nastx
  STOP
111 CONTINUE

  CALL filclo(iunpneas,' ')
  
END PROGRAM propag_PNEAs


! ================================================================
! **** SUBROUTINES ****
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
