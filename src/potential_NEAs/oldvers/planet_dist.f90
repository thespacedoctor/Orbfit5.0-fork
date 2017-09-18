! compute the orbit distance with the planets from Mercury to Saturn
! taking into account the secular evolution of the asteroid orbit
PROGRAM planet_dist
  USE fund_const
  USE output_control
  USE orbit_elements
  USE critical_points
  USE planetorb_circ
  USE secular_evolution

  IMPLICIT NONE
  INTEGER, PARAMETER :: norbx=1000000
  INTEGER, PARAMETER :: npla=8
  INTEGER, PARAMETER :: hmax=10000
  CHARACTER*80 :: catnam
  TYPE(orbit_elem) :: eq,elkep
  CHARACTER(LEN=9) :: name ! asteroid name 
  REAL(KIND=dkind) :: time,tmjd 
  LOGICAL ::  eof
  INTEGER :: norb,loop
  INTEGER :: iunin,iunout
  TYPE(orb_uncert) :: unc
! to read mev/pro files
  INTEGER :: lnam,nfin,ntot
  REAL(KIND=dkind) :: t0,a,Rbar,Z 



  REAL(KIND=dkind),DIMENSION(hmax) :: b2,c2,d2
  REAL(KIND=dkind),DIMENSION(hmax) :: b3,c3,d3
  REAL(KIND=dkind),DIMENSION(hmax) :: b4,c4,d4
  REAL(KIND=dkind),DIMENSION(hmax) :: b5,c5,d5
  REAL(KIND=dkind) :: deltat,dtMJD ! length of fundamental interval (yr,d)
  REAL(KIND=dkind),DIMENSION(hmax) :: t_rescal
  TYPE(sec_evol),DIMENSION(hmax) :: secev
  REAL(KIND=dkind) :: tau 
  INTEGER :: nstep 
  INTEGER,DIMENSION(hmax) :: hsrt
  REAL(KIND=dkind),DIMENSION(4,hmax) :: elem !ecc,inc,Omnod,omega
  REAL(KIND=dkind) :: omega1,omega2
  LOGICAL :: circ ! circ=.true.  for circulators
                  ! circ=.false. for librators

! for dmintil_rms
  INTEGER :: nummin
  REAL(KIND=dkind):: dmintil(nminx)
  REAL(KIND=dkind):: dmin(npla)

  INTEGER :: j,h,nn ! loop indexes

  rhs=1

  CALL filopn(iunout,'pladist.out','unknown')

! opening catalog of orbits
  catnam='./neodys.ctc'
  CALL filopn(iunin,catnam,'OLD')
  CALL oporbf(catnam,iunin)

! reading planet data (elpl%coo='COM')
  CALL placirc

! main loop
  norb=0
  DO loop=1,1!norbx
     CALL read_elems(eq,name,eof,catnam,iunin,unc)
     IF(eof) EXIT
     write(*,*)name     

     CALL rmsp(name,lnam) 
     OPEN (10,file='mev_files/'//name(1:lnam)//'.mev',status='old')
     READ(10,100) t0,a,Rbar,Z 
     secev(1:hmax)%a = a

! *********************
     eq%coo='KEP'
     eq%coord(1)=a
     eq%coord(6)=0.d0
! *********************

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

! *** pro files ***
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

! search for fundamental interval
     CALL fund_int(secev(1:nfin),nfin,ntot)
!     write(*,*)'fundamental interval'

! check extremum values of omega
     IF(secev(1)%om.gt.359.9d0)THEN
        secev(1)%om=secev(1)%om-360.d0 
     ENDIF
     IF(secev(ntot)%om.lt.1d0)THEN
        secev(ntot)%om=secev(ntot)%om+360.d0 
     ENDIF

     t_rescal(1:ntot)=secev(1:ntot)%t-secev(1)%t
     write(*,*)'t_trasl=',secev(1)%t

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


     DO h=1,FLOOR(deltat/10)+1 !2000
        time = 2010+(h-1)*10
        write(*,*)time
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

! elem(1:4) = (ecc,inc,Omnod,omega)
!        write(*,*)elem(1:4,h)
        eq%coord(2)=elem(1,h)
        eq%coord(3)=elem(2,h)*radeg
        eq%coord(4)=elem(3,h)*radeg
        eq%coord(5)=elem(4,h)*radeg

        DO j=1,npla
 !          write(*,*)eq%coord(1:5)
           CALL dmintil_rms(elpl(j),eq,nummin,dmintil)
           dmin(j)=dmintil(1)
        ENDDO
!        write(*,111) name,eq%coord(1),dmin(1:3)
        write(iunout,111) name,eq%coord(1),dmin(1:npla)
        
     ENDDO
     
  ENDDO
111 FORMAT(a9,2x,f10.5,2x,8(f10.5,1x))

  CALL filclo(iunout,' ')


END PROGRAM planet_dist
