! compute the orbit distance with the planets from Mercury to Saturn
! taking into account the secular evolution of the asteroid orbit
PROGRAM planet_dist
  USE fund_const
  USE output_control
  USE orbit_elements
  USE planetorb_circ
  USE secular_evolution
  USE critical_points, ONLY:nminx,dmintil_rms
  USE right_hand_side !, ONLY :inpl,ioupl
  
  IMPLICIT NONE
  INTEGER, PARAMETER :: norbx=1000000
  INTEGER, PARAMETER :: npla=8
  INTEGER, PARAMETER :: hmax=10000
!  CHARACTER*80 :: catnam

  CHARACTER(LEN=9) :: name,namsel! asteroid name 
  CHARACTER*1 :: lc
  CHARACTER*60 :: elefil !file names
  INTEGER :: le

  TYPE(orbit_elem) :: eq,elkep
  REAL(KIND=dkind) :: time0,time1,time,timed
  LOGICAL ::  eof
  INTEGER :: norb,loop
  INTEGER :: iunin,iunout,ioutel,iunt0,iuntdg,iuncro,iunmars,ierrsgn
  TYPE(orb_uncert) :: unc
! to read mev/pro files
  INTEGER :: lnam,nfin,ntot
  REAL(KIND=dkind) :: t0,a,Unorm(8)

  REAL(KIND=dkind),DIMENSION(hmax) :: b2,c2,d2
  REAL(KIND=dkind),DIMENSION(hmax) :: b3,c3,d3
  REAL(KIND=dkind),DIMENSION(hmax) :: b4,c4,d4
  REAL(KIND=dkind),DIMENSION(hmax) :: b5,c5,d5
  REAL(KIND=dkind) :: deltat,dtMJD ! length of fundamental interval (yr,d)
  REAL(KIND=dkind) :: omegaper ! period of circulation/libration of omega
  TYPE(sec_evol),DIMENSION(hmax) :: secev
  INTEGER :: nsamp 
  INTEGER,DIMENSION(hmax) :: hsrt
  REAL(KIND=dkind),DIMENSION(4) :: elem !ecc,inc,Omnod,omega
  REAL(KIND=dkind) :: omega1,omega2,deltanod,Omnod2
  LOGICAL :: circ ! circ=.true.  for circulators
                  ! circ=.false. for librators
! for dmintil_rms
  INTEGER :: nummin
  REAL(KIND=dkind):: dmintil(nminx),deriv(nminx),c1min(3,nminx),c2min(3,nminx)
  REAL(KIND=dkind):: tau_1(3,nminx),tau_2(3,nminx)
  REAL(KIND=dkind),DIMENSION(8):: dnodp,dnodm,dnodpprev,dnodmprev 
  REAL(KIND=dkind):: dmin(npla),dminder(npla) !dmin and its time derivative
  REAL(KIND=dkind):: dminprev(npla) !dmin at previous iteration
  INTEGER,DIMENSION(nminx) :: dminflag
  REAL(KIND=dkind):: tdanger(npla),t_sum,t_tot,time_len,time_spent
! for the second minimum value
  REAL(KIND=dkind):: dmin_2(npla),dminder_2(npla),tdanger_2(npla)
  INTEGER :: j,h,nn ! loop indexes
  REAL(KIND=dkind),DIMENSION(8) :: gm
  REAL(KIND=dkind),DIMENSION(8) :: delta,Tiss
  REAL(KIND=dkind) :: mult_fact
  LOGICAL :: first,cross!,crossmars
 
  CALL pladat
  rhs=1

! total time span
  time_len = 1.d6  

 WRITE(*,*)'choose asteroid'
 READ(*,*) namsel

! output files
!  CALL filopn(ioutel,'elements.out','unknown')
  CALL filopn(ierrou,'pladist.err','unknown')
  CALL filopn(iunt0,'tempoinit','unknown')
  CALL filopn(iuntdg,'dangertime','unknown')
  CALL filopn(iuncro,'NOcros','unknown')
!  CALL filopn(iunmars,'onlyMars','unknown')
  CALL filopn(ierrsgn,'sign.err','unknown')

! reading planet data 
  CALL placirc         ! planet trajectories (elpl%coo='KEP')
  INCLUDE 'plamas.h90' ! planet masses in Sun mass units

  CALL filopn(iunin,'neos.data','old')
  READ(iunin,*) ! skip header
105 FORMAT(a9,1x,a1)

! main loop
  norb=0
  DO loop=1,norbx

     dminflag=0 !initialization
!     CALL read_elems(eq,name,eof,catnam,iunin,unc)
!     IF(eof) EXIT
     READ(iunin,105,END=33) name,lc
     CALL rmsp(name,lnam)                                                      


! skipping these cases
     IF(name.eq.'2007VA85'.OR. &
          & name.eq.'2002AA29'.OR. &
          & name.eq.'2009AV')THEN        
        CYCLE
     ENDIF

     IF(name.ne.namsel)THEN        
        CYCLE
     ELSE
        WRITE(*,*)name
     ENDIF
     CALL filnam('.',name,'dmin',elefil,le)
     write(*,*)'planet_dist: opening output file ',elefil(1:le)
     CALL filopn(iunout,elefil(1:le),'unknown')

!     WRITE(iunout,151)name
     WRITE(ierrou,*)name
!  WRITE(ioutel,*)name
151  FORMAT (a9,2x,10('0.d0',5x))

     CALL sec_interp(name,deltat,t0,a,ntot,secev, &
          & omega1,omega2,deltanod,Omnod2,b2,c2,d2,b3,c3,d3,b4,c4,d4,b5,c5,d5)
     dtMJD=deltat*365.25d0 ! length of fund. interval (d) 

     IF(omega1.gt.omega2)THEN
        WRITE(ierrou,*)'backward circulation of omega'
        WRITE(ierrou,*)name,omega1,omega2
        CYCLE
     ENDIF

     IF(deltat.lt.0.d0)THEN
        WRITE(*,*)'ERROR! negative deltat',name,deltat
        WRITE(ierrou,*)'ERROR! negative deltat',name,deltat
        numerr=numerr+1
        CYCLE
!        STOP
     ENDIF


     IF (lc.eq.'C')THEN
        circ=.true.
        omegaper = 4.d0*deltat
     ELSEIF (lc.eq.'L'.OR.lc.eq.'l')THEN
        circ=.false.
        omegaper = 2.d0*deltat
     ENDIF

! compute critical distances delta_j 
     mult_fact=10.d0
     write(*,*)'critical distances delta_j'
     DO j=1,npla
        Tiss(j) = elpl(j)%coord(1)/a + &
             & 2.d0*sqrt(a/elpl(j)%coord(1)* &
             & (1.d0-secev(1)%ecc**2))*cos(secev(1)%inc*radeg)
        IF(Tiss(j).lt.3.d0)THEN
           Unorm(j) = sqrt(3.d0- Tiss(j))
           delta(j)=mult_fact*gm(j)/Unorm(j)**2
           delta(j) = delta(j)*elpl(j)%coord(1) ! conversion into AU
           write(*,*)j,delta(j)
        ELSE
           delta(j) = 0.d0
           CYCLE
        ENDIF
     ENDDO

! initialization of asteroid elements
     eq=undefined_orbit_elem
     eq%t=t0
     eq%coo='KEP'
     eq%coord(1)=a
     eq%coord(6)=0.d0

     t_tot=0.d0 !initialization     
     first=.true.
! time evolution loop
     CALL mjd2yr(t0,time0)         !starting time of mevfile (yr)
     CALL mjd2yr(secev(1)%t,time1) !starting time of fund interval (yr)
     write(iunt0,*) time0
     write(*,*)'starting time of fund interval (yr)=',time1
     nsamp=3000

!     crossmars=.false.
     cross=.false.
     DO h=0,nsamp+1 
        time = time1 + (h-1)*(deltat/nsamp)
!        write(*,*)'h=',h,'time=',time

! ********************************************************
! make computation only in the fundamental interval range
!        CALL yr2mjd(time,timed)
!        IF(timed.lt.secev(1)%t-(deltat/nsamp))THEN
!           WRITE(*,*)'skip time',secev(1)%t,timed
!           CYCLE
!        ENDIF
!        IF(time.lt.time1)THEN
!           WRITE(*,*)'time less than time1',time,time1
!           CYCLE
!        ELSEIF(time.gt.time1+deltat)THEN
!           WRITE(*,*)'time greater than time2',time,time1+deltat
!           CYCLE
!        ENDIF
! ********************************************************

! compute secular evoluion at t=time
        CALL comp_interp(time,dtMJD,circ,hmax,ntot,secev,omega2,deltanod,&
             & Omnod2, &
        & b2,c2,d2,b3,c3,d3,b4,c4,d4,b5,c5,d5,elem)
        eq%coord(2)=elem(1)         !ecc
        eq%coord(3)=elem(2)*radeg   !inc
        eq%coord(4)=elem(3)*radeg   !Omnod
        eq%coord(5)=elem(4)*radeg   !omega
!        write(ioutel,108)time,eq%coord(1:5)

! compute d_{min}
        t_sum=0.d0
        DO j=1,npla
! orbit distance
           CALL dmintil_rms(elpl(j),eq,nummin,dmintil,DMINFLAG=dminflag)
           !           write(*,*)'dminflag=',dminflag
           IF(dminflag(1).gt.0.OR.dminflag(2).gt.0)THEN
              write(ierrou,*)'critical points failure: planet=',j, &
                   & 'dminflag=',dminflag(1:2)
              write(ierrou,('(5(f12.5,1x))')) elpl(j)%coord(1:5)
              write(ierrou,('(5(f12.5,1x))')) eq%coord(1:5)
              IF(dminflag(1).gt.1.OR.dminflag(2).gt.1)THEN
                 write(ierrsgn,*)'sign_dmin failure:ast=',name,'planet=',j, &
                   & 'dminflag=',dminflag(1:2)
                 write(ierrsgn,('(5(f12.5,1x))')) elpl(j)%coord(1:5)
                 write(ierrsgn,('(5(f12.5,1x))')) eq%coord(1:5)
              ENDIF
              GOTO 333
           ENDIF
           dmin(j)=dmintil(1)
           IF(nummin.gt.1)THEN
              dmin_2(j)=dmintil(2)
           ENDIF
! nodal distances
           dnodp(j) = elpl(j)%coord(1)-a*(1.d0-eq%coord(2)**2)/ &
                & (1.d0+eq%coord(2)*cos(eq%coord(5)))
           dnodm(j) = elpl(j)%coord(1)-a*(1.d0-eq%coord(2)**2)/ &
                & (1.d0-eq%coord(2)*cos(eq%coord(5)))
 !          write(*,*) j,dnodp(j),dnodm(j)
           IF(.not.first)THEN
              IF(dnodp(j)*dnodpprev(j).lt.0.d0.OR. &
                   & dnodm(j)*dnodmprev(j).lt.0.d0)THEN
                 cross=.true.
                 
                 CALL dmintil_deriv(elpl(j),eq,nummin,dmintil,deriv) 
                 dminder(j) = deriv(1)
                 tdanger(j)=2.d0*delta(j)/abs(dminder(j))
                 write(*,*)'***j,time,dmin,dminprev,tdanger'
                 write(*,('(i2,2x,f12.5,1x,2(es12.5,1x),f15.5)'))j,time,dmin(j),dminprev(j),tdanger(j)
                 write(*,*)'j,dnodp,dnodpprev', &
                      & j,dnodp(j),dnodpprev(j)
                 write(*,*)'j,dnodm,dnodmprev', &
                      & j,dnodm(j),dnodmprev(j)
                 t_sum=t_sum + tdanger(j)
                 
                 IF(nummin.gt.1)THEN
                    IF(dnodp(j)*dnodpprev(j).lt.0.d0.AND. &
                         & dnodm(j)*dnodmprev(j).lt.0.d0)THEN
                       write(*,*)'using also the second minimum value'
                       dminder_2(j) = deriv(2)
                       tdanger_2(j)=2.d0*delta(j)/abs(dminder(j))
                       t_sum=t_sum + tdanger_2(j)
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
           dminprev(j)=dmin(j)
           dnodpprev(j)=dnodp(j)
           dnodmprev(j)=dnodm(j)
        ENDDO
        write(iunout,111) h,eq%coord(1),dmin(1:npla)
        first=.false.
        t_tot = t_tot + t_sum 

333     CONTINUE        
     ENDDO !end loop on time evolution
     IF(.not.cross)THEN
        write(iuncro,*)name
     ENDIF
     
     time_spent = time_len*t_tot/deltat 
     write(*,*)'omegaper,t_tot,deltat=',omegaper,t_tot,deltat
     write(*,*)'time_spent=',time_spent
     write(iuntdg,*)name,time_spent

  CALL filclo(iunout,' ')
  ENDDO !end loop on asteroids   

  WRITE(*,*)'error! norbx too small:',norbx
33 CONTINUE

108 FORMAT(f12.5,2x,5(f10.5,1x))
111 FORMAT(i4,2x,f10.5,2x,8(f10.5,1x))

  CALL filclo(iuntdg,' ')  
  CALL filclo(iunt0,' ')
  CALL filclo(ierrou,' ')
  CALL filclo(ierrsgn,' ')
!  CALL filclo(ioutel,' ')
  CALL filclo(iuncro,' ')
!  CALL filclo(iunmars,' ')


END PROGRAM planet_dist
