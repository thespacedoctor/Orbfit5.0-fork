MODULE right_hand_side
  USE fund_const
  USE orbit_elements
  IMPLICIT NONE
  PRIVATE
! asteroid/planet orbital elements
  REAL(KIND=dkind) :: aa,ea,beta,ci,si,apl,ep,betap,gmp
  REAL(KIND=dkind) :: aplr,epr,betapr,gmpr ! resonant
  TYPE(orbit_elem) :: elpl !resonant
  REAL(KIND=dkind) :: om,omnod,g,zl,a
! rotational matrix
  REAL(KIND=dkind),DIMENSION(3,3) :: R_omp,R_ip,Hp,R_om,R_i,R_on,Hast
  REAL(KIND=dkind),DIMENSION(3,3) :: R_ompr,R_ipr,Hpr
  REAL(KIND=dkind),DIMENSION(3,3) :: dR_om,dR_i,dR_on,dHdom,dHdI,dHdon
  REAL(KIND=dkind),DIMENSION(3,3) :: ddHdI
! coord vectors and derivatives
  REAL(KIND=dkind) :: x(3),y(3),cu,su,l,ckl,skl,usel
  REAL(KIND=dkind),DIMENSION(3) :: dax,dex,dIx,domx,donx
  REAL(KIND=dkind),DIMENSION(3) :: ddex,ddeIx,ddIx,ddomx,ddonx
! for selecting derivatives
  INTEGER :: ider,res_flag,sin_flag
  !resonance
  INTEGER :: kappaa, kappap, n_res, numpla
  REAL(KIND=dkind),PARAMETER :: epsabs = 1.d-9,epsrel=1.d-7

! planet data                                                       
  INCLUDE 'pldata.h90'
  
  PUBLIC :: kappaa,kappap,n_res,gmp,sin_flag
  
! routines
  PUBLIC :: rhs2,secpert,ffd_pert_res,rot_matrix !derrhs

CONTAINS
! =================================================================
! RIGHT SIDE COMPUTATION                                       
! INTEGRATION DOMAIN: T = [0,DPIG] x [0,DPIG]                       
!SUBROUTINE rhs2(n,elpl,om,omnod,g,zl,a,ddd,eee,nnn) 
 SUBROUTINE rhs2(n,elpl0,n_res,om0,omnod0,g0,zl0,bigelle,sigma,ddd,eee,nnn)    ! resonant cfr linae sopra
! --------------- interface ---------------------
  INTEGER,INTENT(IN) :: n,n_res !planets number,resonant planet number
  TYPE(orbit_elem),INTENT(IN) :: elpl0 !planet elems 
!  TYPE(orbit_elem),INTENT(IN) :: elplr0 !resonant planet elems 
  REAL(KIND=dkind),INTENT(IN) :: om0,omnod0,g0,zl0,bigelle,sigma   !resonant add sigma
  REAL(KIND=dkind),INTENT(OUT) :: ddd(6),eee(6)     !resonant 6->4
  INTEGER,INTENT(OUT) :: nnn(6)                      !resonant "
! -------------- end interface -------------------
! for dqags                                                         
  INTEGER :: limx,limx4 
  PARAMETER (limx=500,limx4=4*limx)
  INTEGER :: ier,limit,iwork(limx),lenw,last 
! function evaluations                                              
  INTEGER :: neval,neval_res_c,neval_res_s 
  REAL(KIND=dkind) :: abserr, abserr_res_c,abserr_res_s,work(limx4) 
! output of dqags                                                   
  REAL(KIND=dkind) :: rm,rm_res_c,rm_res_s
! loop indexes                                                      
  INTEGER :: k,i
! mean motion resonant planet
  REAL(KIND=dkind) :: np
  REAL(KIND=dkind) :: par_elim
 ! ======================================================

  par_elim = 1.d0
  elpl = elpl0
!  elplr = elplr0
  om = om0
  omnod = omnod0
  g = g0
  zl = zl0
  a = (bigelle/ky)**2
  numpla = n
  gmp = gm(numpla)   !gmp = (k^2)*m(pla)/m(Sun)

! compute rotational matrix and Cartesian coords for planet and asteroid
  CALL rot_matrix(elpl,om,omnod,g,zl,a)
! ---- initialization -----
!  DO k=1,4                       
  DO k=1,6   !resonant cfr linea sopra
     ddd(k)=0.d0 
     eee(k)=0.d0 
     nnn(k)=0 
  ENDDO
! ---- preparing dqags call ---
!  epsabs=1.d-8 
!  epsrel=1.d-5 
  limit=limx 
  lenw=limx4 
! ********************************************************************     
!                    COMPUTING DERIVATIVES 
! ********************************************************************     
! first computing integral with respect to u
!  l=u-ea*sin(u)                            !resonant added  
! derivative with respect to G
  ider=1 
  res_flag=0
  CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)                               
  rm_res_c=0
  rm_res_s=0
  abserr_res_c=0
  abserr_res_s=0
  neval_res_c=0
  neval_res_s=0
  IF (numpla.eq.n_res) THEN
     res_flag=1
     sin_flag=0 !multiply by cos(sigma)
     CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm_res_c,abserr_res_c,neval_res_c,ier,    &    
          &        limit,lenw,last,iwork,work)                                         
     sin_flag= 1 !multiply by sin(sigma)
     CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm_res_s,abserr_res_s,neval_res_s,ier,    &
          &        limit,lenw,last,iwork,work)                                         

  ENDIF

! CORRECT FORMULA
  ddd(1)=ddd(1)+(rm+rm_res_c*cos(sigma)+rm_res_s*sin(sigma))/(dpig**2) 

! DUMMY VALUE
!  ddd(1)=ddd(1)+rm/(dpig**2)
!  ddd(1) = rm_res_c
!  ddd(1) = rm_res_s
!  ddd(1)= (rm_res_c*cos(sigma)+rm_res_s*sin(sigma))/(dpig**2) 

  eee(1)=eee(1)+(abserr+abserr_res_c+abserr_res_s)/(dpig**2) 
  nnn(1)=nnn(1)+neval+neval_res_c+neval_res_s 

! derivative with respect to omega                                  
  ider=2 
  res_flag=0
  CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)                               
  rm_res_c=0
  rm_res_s=0
  abserr_res_c=0
  abserr_res_s=0
  neval_res_c=0
  neval_res_s=0
  IF (numpla.eq.n_res) THEN
     res_flag=1
     sin_flag=0
     CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm_res_c,abserr_res_c,neval_res_c,ier,    &      ! resonant added
          &        limit,lenw,last,iwork,work)                                         ! resonant added
     sin_flag= 1 
     CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm_res_s,abserr_res_s,neval_res_s,ier,    &      ! resonant added
          &        limit,lenw,last,iwork,work)                                         ! resonant added

  ENDIF
!      ddd(2)=ddd(2)+rm/(dpig**2)
! CORRECT FORMULA
  ddd(2)=ddd(2)+(rm+rm_res_c*cos(sigma)+rm_res_s*sin(sigma))/(dpig**2) 

! DUMMY VALUE
!  ddd(2) = rm_res_c
!  ddd(2)=ddd(2)+rm/(dpig**2)

  eee(2)=eee(2)+(abserr+abserr_res_c+abserr_res_s)/(dpig**2) 
  nnn(2)=nnn(2)+neval+neval_res_c+neval_res_s
  
  ! derivative with respect to Z                                      
  ider=3 
  res_flag=0
  CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)                               
  rm_res_c=0
  rm_res_s=0
  abserr_res_c=0
  abserr_res_s=0
  neval_res_c=0
  neval_res_s=0
  IF (numpla.eq.n_res) THEN
     res_flag=1
     sin_flag=0
     CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm_res_c,abserr_res_c,neval_res_c,ier,    &      ! resonant added
          &        limit,lenw,last,iwork,work)                                         ! resonant added
     sin_flag= 1 
     CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm_res_s,abserr_res_s,neval_res_s,ier,    &      ! resonant added
          &        limit,lenw,last,iwork,work)                                         ! resonant added
    
  ENDIF
! CORRECT FORMULA
  ddd(3)=ddd(3)+(rm+rm_res_c*cos(sigma)+rm_res_s*sin(sigma))/(dpig**2)
  eee(3)=eee(3)+(abserr+abserr_res_c+abserr_res_s)/(dpig**2) 
  nnn(3)=nnn(3)+neval+neval_res_c+neval_res_s

! DUMMY VALUE
!  ddd(3)=ddd(3)+rm/(dpig**2)

! derivative with respect to Omnod
  ider=4 
  res_flag=0
  CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)                               
  rm_res_c=0
  rm_res_s=0
  abserr_res_c=0
  abserr_res_s=0
  neval_res_c=0
  neval_res_s=0
  IF (numpla.eq.n_res) THEN
     res_flag=1
     sin_flag=0
     CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm_res_c,abserr_res_c,neval_res_c,ier,    &   
          &        limit,lenw,last,iwork,work)                                         
     sin_flag= 1 
     CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm_res_s,abserr_res_s,neval_res_s,ier,    &
          &        limit,lenw,last,iwork,work)                                         
     
  ENDIF
! CORRECT FORMULA
  ddd(4)=ddd(4)+(rm+rm_res_c*cos(sigma)+rm_res_s*sin(sigma))/(dpig**2) 

! DUMMY VALUE
!  ddd(4) = rm_res_c
!  ddd(4)=ddd(4)+rm/(dpig**2)

  eee(4)=eee(4)+(abserr+abserr_res_c+abserr_res_s)/(dpig**2) 
  nnn(4)=nnn(4)+neval+neval_res_c+neval_res_s

! derivative with respect to sigma                                                      ! resonant added
  
!  IF (numpla.eq.n_res) THEN
  ider=5
  res_flag=0
  CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)                               
  rm_res_c=0
  rm_res_s=0
  abserr_res_c=0
  abserr_res_s=0
  neval_res_c=0
  neval_res_s=0
  IF (numpla.eq.n_res) THEN
     res_flag=1
     sin_flag=0
     CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm_res_c,abserr_res_c,neval_res_c,ier,    &      ! resonant added
          &        limit,lenw,last,iwork,work)                                         ! resonant added
     sin_flag= 1 
     CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm_res_s,abserr_res_s,neval_res_s,ier,    &      
          &        limit,lenw,last,iwork,work)                            

! CORRECT FORMULA
     ddd(5)=ddd(5)+(rm_res_s*cos(sigma)-rm_res_c*sin(sigma))/(dpig**2)    
     eee(5)=eee(5)+(abserr_res_c+abserr_res_s)/(dpig**2) 
     nnn(5)=nnn(5)+neval_res_c+neval_res_s
  ENDIF

! derivative with respect to S (=bigelle)
  ider=6
  res_flag=0
  CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)                               
  rm_res_c=0
  rm_res_s=0
  abserr_res_c=0
  abserr_res_s=0
  neval_res_c=0
  neval_res_s=0
  IF (numpla.eq.n_res) THEN
     res_flag=1
     sin_flag=0
     CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm_res_c,abserr_res_c, &
          & neval_res_c,ier,limit,lenw,last,iwork,work)
     sin_flag= 1 
     CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm_res_s,abserr_res_s, &
          & neval_res_s,ier,limit,lenw,last,iwork,work)

     np = ky/elpl%coord(1)**(3.d0/2.d0)                                   

!     write(*,*)'elpl%coo, elpl%coord(1):',elpl%coo,elpl%coord(1)
!     stop

! CORRECT FORMULA
     ddd(6)=ddd(6) + ky**4*kappaa**(-2.d0)*bigelle**(-3.d0)+np*kappap

! DUMMY VALUE
!     ddd(6)= ky**4*kappaa**(-2.d0)*bigelle**(-3.d0)+np*kappap

  ENDIF

! CORRECT FORMULA
  ddd(6)=ddd(6)+(rm+rm_res_c*cos(sigma)+rm_res_s*sin(sigma))/(dpig**2)  
  eee(6)=eee(6)+(abserr+abserr_res_c+abserr_res_s)/(dpig**2) 
  nnn(6)=nnn(6)+neval+neval_res_c+neval_res_s

! DUMMY VALUE!
!  ddd(6)=ddd(6)+rm/(dpig**2)
!  ddd(6) =  rm_res_c  
!  ddd(6) = ((ky**4)*(kappaa**(-2.d0))*(bigelle**(-3.d0)) + np*kappap)/(ioupl-inpl)
!    write(*,134)'rhs2: numpla,ddd(6),ky,kappaa,kappap,bigelle=',numpla,ddd(6),ky,kappaa,kappap,bigelle
!134 format(a33,i2,2x,2(f10.5,2x),(2i2,2x),f10.5)

END SUBROUTINE rhs2
                                                                        
! ********************************************************************     
! subroutine computing unidim.average of the unidim.average         
DOUBLE PRECISION FUNCTION ffd(u) 
! --------------- interface ---------------------
  REAL(KIND=dkind) :: u !eccentric anomaly of the asteroid
! ------------ end interface ----------------------
  REAL(KIND=dkind) :: dy(3),day(3)                             ! resonant added day
! for dqags                                                         
  INTEGER :: limx,limx4 
  PARAMETER (limx=500,limx4=4*limx) 
  INTEGER :: neval,ier,limit,iwork(limx),lenw,last 
  REAL(KIND=dkind) :: abserr,work(limx4),risult 
! --------------------------------------------------

  cu=cos(u)
  su=sin(u)
  usel = u
! Asteroid coordinates
  y(1) = aa*(cu-ea)
  y(2) = aa*su*beta
  y(3) = 0.d0
  CALL prodmv(x,Hast,y)
! 1th-derivatives w.r.t. a,ecc,Inc,omega,Omnod for asteroid coord
! derivatives w.r.t. ecc
  dy(1) =-aa
  dy(2) =-aa*su*ea/beta
  dy(3) = 0.d0      
! derivatives w.r.t. a
  day(1) =cu-ea                  ! resonant added               
  day(2) =su*beta                 ! resonant added
  day(3) = 0.d0                      ! resonant added
  CALL prodmv(dax,Hast,day)         ! resonant added
  CALL prodmv(dex,Hast,dy)
  CALL prodmv(dIx,dHdI,y)
  CALL prodmv(domx,dHdom,y) 
  CALL prodmv(donx,dHdon,y)
! ------------ preparing dqags call ------------
!  epsabs=1.d-8 
!  epsrel=1.d-5 
  limit=limx 
  lenw=limx4
  IF (res_flag.eq.0) THEN ! parte non risonante
!  ------------ selecting derivative  ------------
     IF(ider.eq.1)THEN 
        call dqagsc(fg,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
             &        ier,limit,lenw,last,iwork,work)                           
        ! if(ier.gt.1)write(*,*)'d/G ',ier,u 
     ELSEIF(ider.eq.2)THEN 
        call dqagsc(fom,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
             &        ier,limit,lenw,last,iwork,work)                           
        ! if(ier.gt.1)write(*,*)'d/dom ',ier,u                              
     ELSEIF(ider.eq.3)THEN 
        call dqagsc(fz,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
             &        ier,limit,lenw,last,iwork,work)                           
        ! if(ier.gt.1)write(*,*)'d/dZ ',ier,u
     ELSEIF(ider.eq.4)THEN 
        call dqagsc(fon,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
             &        ier,limit,lenw,last,iwork,work)                           
        ! if(ier.gt.1)write(*,*)'d/donod ',ier,u
     ELSEIF(ider.eq.5)THEN                                               
! in this case the computation of the integral is done only for res_flag=1
        risult = 0.d0
        abserr = 0.d0
        neval = 0
     ELSEIF(ider.eq.6)THEN                                               
        call dqagsc(fs,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    & 
             &        ier,limit,lenw,last,iwork,work)                    
        ! if(ier.gt.1)write(*,*)'d/donod ',ier,u
     ELSE 
        WRITE(*,*)' ffd: ider=',ider 
        STOP 
     ENDIF
      
     ffd=risult         

  ELSE ! res_flag = 1
     !  l=u-ea*sin(u)  
     IF(ider.eq.1)THEN 
        call dqagsc(fg_res,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
             &        ier,limit,lenw,last,iwork,work)                           
        ! if(ier.gt.1)write(*,*)'d/G ',ier,u 
     ELSEIF(ider.eq.2)THEN 
        call dqagsc(fom_res,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
             &        ier,limit,lenw,last,iwork,work)                           
        ! if(ier.gt.1)write(*,*)'d/dom ',ier,u                              
     ELSEIF(ider.eq.3)THEN 
        call dqagsc(fz_res,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
             &        ier,limit,lenw,last,iwork,work)                           
        ! if(ier.gt.1)write(*,*)'d/dZ ',ier,u
     ELSEIF(ider.eq.4)THEN 
        call dqagsc(fon_res,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
             &        ier,limit,lenw,last,iwork,work)                           
        ! if(ier.gt.1)write(*,*)'d/donod ',ier,u
     ELSEIF(ider.eq.5)THEN
!        write(*,*)'x in ffd:',x

        call dqagsc(fsig_res,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
             &        ier,limit,lenw,last,iwork,work)
     ELSEIF(ider.eq.6)THEN
        call dqagsc(fs_res,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
             &        ier,limit,lenw,last,iwork,work)
        ! if(ier.gt.1)write(*,*)'d/donod ',ier,u
     ELSE 
        WRITE(*,*)' ffd: ider=',ider 
        STOP 
     ENDIF

     ffd=risult         
  ENDIF
  
  RETURN 
END FUNCTION ffd


! subroutine computing G derivative (d(1/d)/dG)
DOUBLE PRECISION FUNCTION fg(up)
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface -----------------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: d2,d3,drde,drdi
!  REAL(KIND=dkind) :: inc,f0,f0_inc,g,g_inc,maux(3,3) for check deriv
! --------------------------------------------------
  cup=cos(up)
  sup=sin(up)

! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  CALL prodmv(xp,Hp,yp)

  d2=DOT_PRODUCT(xp-x,xp-x)
  d3=d2**(1.5d0)

! ========= DERIVATA with respect to e ========================
! ((-1/2)*der(D^2)*(1-e2*cosu)+D^2*der(1-e2*cosu))*(1-e1*cos(up))
! (-1/2)*der(D^2) = <Xpl-Xast, R_on*R_i*R_om*dex>
  drde = DOT_PRODUCT(xp-x,dex)*(1.d0-ea*cu)-cu*d2


! ========= DERIVATA with respect to I =========================
! der(D^2)*(1-e2*cosu)*(1-e1*cos(up))
! der(D^2) = -2*<Xpl-Xast, R_on*dR_i*R_om*Y>
  drdi = DOT_PRODUCT(xp-x,dIx)*(1.d0-ea*cu)


! ========= DERIVATA with respect to G ==========================
     fg = -(gmp/d3)*(1.d0-ep*cup)*((-beta/(ky*dsqrt(aa)*ea))*drde + &
          &   (ci/(si*ky*dsqrt(aa)*beta))*drdi)


! **************************************************
!  CALL check_deriv(cup,xp,d2,fg)
! *************************************************

  RETURN 
END FUNCTION fg


! subroutine computing G derivative (d(1/d5-(r*r5)/|r5|^3)/dG)*(sin oppure cos)(critical angle)
DOUBLE PRECISION FUNCTION fg_res(up)
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface -----------------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: dist,d2,d3,drde,drdI,dindde,dinddI
  REAL(KIND=dkind) :: auxres,derauxres,derskl,derckl
  REAL(KIND=dkind) :: vsize 

!  REAL(KIND=dkind) :: inc,f0,f0_inc,g,g_inc,maux(3,3) for check deriv
! --------------------------------------------------
  cup=cos(up)
  sup=sin(up)

! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  CALL prodmv(xp,Hp,yp)

  d2=DOT_PRODUCT(xp-x,xp-x)
  d3=d2**(1.5d0)
  dist = sqrt(d2)

! ========= DERIVATA with respect to e ========================
!  der((1/D)*(1-ea*cu)) = -1/(2*D^3)*der(D^2)*(1-ea*cu)-cu*(1/D)= 
!  = 1/D^3*( (-1/2)*der(D^2)*(1-ea*cu) - cu*(D^2/D^3))
! (-1/2)*der(D^2) = <Xpl-Xast, R_on*R_i*R_om*dex>

! CORRECT FORMULA
  drde = DOT_PRODUCT(xp-x,dex)*(1.d0-ea*cu) -cu*d2

! DUMMY VALUE
!  drde = 0.d0

! d/de((r*r5)/|r5|^3*(1-e*cosu))
! CORRECT VALUE
  dindde = ( DOT_PRODUCT(dex,xp)*(1.d0-ea*cu) - &
       & DOT_PRODUCT(x,xp)*cu )/vsize(xp)**3

! DUMMY VALUE
!  dindde = 0.d0

! ========= DERIVATA with respect to I =========================
! der(1/D*(1-ea*cu)) = -1/(2*D^3)*der(D^2)*(1-ea*cu)
! (-1/2)*der(D^2) = <Xpl-Xast, R_on*dR_i*R_om*Y>
! CORRECT FORMULA
  drdI = DOT_PRODUCT(xp-x,dIx)*(1.d0-ea*cu)

! DUMMY VALUE
!  drdi = 0.d0


! d/dI((r*r5)/|r5|^3*(1-e*cosu))
! CORRECT VALUE
  dinddI = DOT_PRODUCT(dIx,xp)*(1.d0-ea*cu)/vsize(xp)**3

! DUMMY VALUE
!  dinddI = 0.d0


!(1/D - r*r5/|r5|^3)*(1-e*cu)*(1-ep*cup)
! CORRECT VALUE
 auxres = -2.d0*gmp*(1/dist - DOT_PRODUCT(x,xp)/vsize(xp)**3)* &
       & (1.d0-ea*cu)*(1.d0-ep*cup)

! DUMMY VALUE
! auxres = (1.d0-ea*cu)*(1.d0-ep*cup)
! auxres = -2.d0*gmp*(-DOT_PRODUCT(x,xp)/vsize(xp)**3)* &
!       & (1.d0-ea*cu)*(1.d0-ep*cup)
!  auxres = -2.d0*gmp*(1/dist)*(1.d0-ea*cu)*(1.d0-ep*cup)


! ========= DERIVATA with respect to G ==========================
! CORRECT VALUE
  derauxres = -2.d0*gmp*(1.d0-ep*cup)*( &
       & (-beta/(ky*sqrt(aa)*ea)) * (drde/d3 - dindde) + &
       & (ci/(si*ky*sqrt(aa)*beta)) * (drdI/d3 - dinddI) )   

! DUMMY VALUE
!  derauxres= -cu*(1.d0-ep*cup)*(-beta/(ky*dsqrt(aa)*ea))
!  derauxres = -2.d0*gmp*(1.d0-ep*cup)*( &
!       & (-beta/(ky*dsqrt(aa)*ea)) * ( - dindde) + &
!       & (ci/(si*ky*dsqrt(aa)*beta)) * ( - dinddI) )   
!  derauxres = -2.d0*gmp*(1.d0-ep*cup)*( &
!       & (-beta/(ky*dsqrt(aa)*ea)) * (drde/d3) + &
!       & (ci/(si*ky*dsqrt(aa)*beta)) * (drdI/d3) )   

     IF (sin_flag.eq.0) THEN   
        ckl=cos(kappaa*(usel-ea*su)+kappap*(up-ep*sup)) 
        derckl = sin(kappaa*(usel-ea*su)+kappap*(up-ep*sup))*kappaa*su
        fg_res= derauxres*ckl + auxres*derckl*(-beta/(ky*sqrt(aa)*ea))
     ELSE                     
        skl=sin(kappaa*(usel-ea*su)+kappap*(up-ep*sup))   
        derskl = -cos(kappaa*(usel-ea*su)+kappap*(up-ep*sup))*kappaa*su
        fg_res= derauxres*skl + auxres*derskl*(-beta/(ky*sqrt(aa)*ea))  
     ENDIF   


! **************************************************
!  CALL check_deriv(cup,xp,d2,fg)
! *************************************************

  RETURN 
END FUNCTION fg_res


! ***********************************************************************
! subroutine computing omega derivative (d(1/d)/dg)
DOUBLE PRECISION FUNCTION fom(up)
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface -------------------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: d2,d3,drdom
! ----------------------------------------------------
  cup=cos(up)
  sup=sin(up)

! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  CALL prodmv(xp,Hp,yp)

  d2 = DOT_PRODUCT(xp-x,xp-x)
  d3 = d2**(1.5d0)

! (-1/2)*der(D^2) = -2*<Xpl-Xast, R_on*R_i*dR_om*Y>
  drdom=DOT_PRODUCT(xp-x,domx)

!write(*,*)'drdom:',drdom
!stop

! =========== DERIVATIVE with respect to omega ============== 
  fom = (gmp/d3)*(1.d0-ea*cu)*(1.d0-ep*cup)*drdom

! **************************************************
!  CALL check_deriv(cup,xp,d2,fom)
! **************************************************
  RETURN 
END FUNCTION fom
                
! subroutine computing omega derivative (d(1/d5-(r*r5)/|r5|^3)/dg)*(sin oppure cos)(critical angle)
DOUBLE PRECISION FUNCTION fom_res(up)
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface -------------------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: d2,d3,drdom,dindom
  REAL(KIND=dkind) :: auxres
  REAL(KIND=dkind) :: vsize
! ----------------------------------------------------
  cup=cos(up)
  sup=sin(up)

! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  CALL prodmv(xp,Hp,yp)

  d2 = DOT_PRODUCT(xp-x,xp-x)
  d3 = d2**(1.5d0)

! der(1/D) = -1/(2*D^3)*der(D^2) = 1/D^3*(-1/2)*der(D^2) = 
! (-1/2)*der(D^2) = <Xpl-Xast, R_on*R_i*dR_om*Y>
  drdom=DOT_PRODUCT(xp-x,domx)

! d/dom((r*r5)/|r5|^3)
  dindom=DOT_PRODUCT(domx,xp)/(vsize(xp)**3)


! =========== DERIVATIVE with respect to omega ============== 
! CORRECT FORMULA
  auxres = 2.d0*gmp*(1.d0-ea*cu)*(1.d0-ep*cup)*(drdom/d3 - dindom)

! DUMMY VALUE (without indirect)
! auxres = 2.d0*gmp*(1.d0-ea*cu)*(1.d0-ep*cup)*(drdom/d3)

!write(*,*)'drdom res:',drdom
!stop

  IF (sin_flag.eq.0) THEN    ! resonant added
     ckl=cos(kappaa*(usel-ea*su)+kappap*(up-ep*sup))    !resonant_added
     fom_res=auxres*ckl       ! resonant added
  ELSE                     ! resonant added
     skl=sin(kappaa*(usel-ea*su)+kappap*(up-ep*sup))    !resonant_added
     fom_res=auxres*skl        ! resonant added
  ENDIF

! DUMMY VALUE
!fom_res=auxres

! **************************************************
!  CALL check_deriv(cup,xp,d2,fom)
! **************************************************
  RETURN 
END FUNCTION fom_res

                                                    
! ********************************************************************      
! subroutine computing Z derivative (d(1/d)/dZ)
DOUBLE PRECISION FUNCTION fz(up) 
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface ------------------------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: d2,d3,drdi
! ---------------------------------------------------------
  cup=cos(up)
  sup=sin(up)

! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  CALL prodmv(xp,Hp,yp)

  d2 = DOT_PRODUCT(xp-x,xp-x)
  d3 = d2**(1.5d0)

! ============ DERIVATIVE with respect to I ===============
! (-1/2)*der(D^2) = <Xpl-Xast, R_on*dR_i*R_om*Y>
  drdi=DOT_PRODUCT((xp-x),dIx) 

! ========== DERIVATIVE with respect to Z ==================
  fz = -(gmp/d3)*(-1.d0/(ky*beta*sqrt(aa)*si))*(1.d0-ea*cu)*(1.d0-ep*cup)*drdi

! **************************************************
!  CALL check_deriv(cup,xp,d2,fz)
! **************************************************
  RETURN 
END FUNCTION fz


! ********************************************************************      
! subroutine computing Z derivative (d(1/d5-(r*r5)/|r5|^3)/dZ)*(sin oppure cos)(critical angle)
DOUBLE PRECISION FUNCTION fz_res(up) 
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface ------------------------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: d2,d3,drdI,dinddI
  REAL(KIND=dkind) :: auxres
  REAL(KIND=dkind) :: vsize
! ---------------------------------------------------------
  cup=cos(up)
  sup=sin(up)

! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  CALL prodmv(xp,Hp,yp)

  d2 = DOT_PRODUCT(xp-x,xp-x)
  d3 = d2**(1.5d0)

! ============ DERIVATIVE with respect to I ===============
! (-1/2)*der(D^2) = <Xpl-Xast, R_on*dR_i*R_om*Y>
  drdI=DOT_PRODUCT((xp-x),dIx) 

! d/dI((r*r5)/|r5|^3*(1-e*cosu))
  dinddI = DOT_PRODUCT(dIx,xp)/vsize(xp)**3

! ========== DERIVATIVE with respect to Z ==================
  auxres = -2.d0*gmp*(-1.d0/(ky*beta*sqrt(aa)*si))* &
       & (1.d0-ea*cu)*(1.d0-ep*cup)* (drdi/d3 - dinddI)

  IF (sin_flag.eq.0) THEN    ! resonant added
     ckl=cos(kappaa*(usel-ea*su)+kappap*(up-ep*sup))    !resonant_added
     fz_res=auxres*ckl       ! resonant added
  ELSE                     ! resonant added
     skl=sin(kappaa*(usel-ea*su)+kappap*(up-ep*sup))    !resonant_added
     fz_res=auxres*skl        ! resonant added
  ENDIF

! **************************************************
!  CALL check_deriv(cup,xp,d2,fz)
! **************************************************
  RETURN 
END FUNCTION fz_res




! *********************************************************************
! subroutine computing Omnod.nod derivative (d(1/d)/donod)
! This derivative is zero when the orbit of the planet is circular
DOUBLE PRECISION FUNCTION fon(up) 
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface ----------------------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: d2,d3,drdon
! -------------------------------------------------------
  cup=cos(up)
  sup=sin(up)
! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  CALL prodmv(xp,Hp,yp)
   
  d2 = DOT_PRODUCT(xp-x,xp-x)
  d3 = d2**(1.5d0)

! (-1/2)*der(D^2) = <Xpl-Xast, dR_on*R_i*R_om*Y>
  drdon=DOT_PRODUCT(xp-x,donx)


! ======= DERIVATIVE with respect to Omnod ===========
  fon = (gmp/d3)*(1.d0-ea*cu)*(1.d0-ep*cup)*drdon


! **************************************************
!  CALL check_deriv(cup,xp,d2,fon)
! **************************************************
  RETURN 
END FUNCTION fon


! *********************************************************************
! subroutine computing Omnod.nod derivative (d(1/d5-(r*r5)/|r5|^3)/donod)*(sin oppure cos)(critical angle)
! This derivative is zero when the orbit of the planet is circular
DOUBLE PRECISION FUNCTION fon_res(up) 
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface ----------------------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: d2,d3,drdon,dindon
  REAL(KIND=dkind) :: auxres,vsize
! -------------------------------------------------------
  cup=cos(up)
  sup=sin(up)

! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  CALL prodmv(xp,Hp,yp)
   
  d2 = DOT_PRODUCT(xp-x,xp-x)
  d3 = d2**(1.5d0)

! (-1/2)*der(D^2) = <Xpl-Xast, dR_on*R_i*R_om*Y>
  drdon=DOT_PRODUCT(xp-x,donx)

! d/don((r*r5)/|r5|^3*(1-e*cosu))
  dindon=DOT_PRODUCT(donx,xp)/vsize(xp)**3


! ======= DERIVATIVE with respect to Omnod ===========
! CORRECT VALUE
  auxres = 2.d0*gmp*(1.d0-ea*cu)*(1.d0-ep*cup)* (drdon/d3 - dindon)

! DUMMY VALUE
!  auxres = 2.d0*gmp*(1.d0-ea*cu)*(1.d0-ep*cup)*(- dindon)

  IF (sin_flag.eq.0) THEN    ! resonant added
     ckl=cos(kappaa*(usel-ea*su)+kappap*(up-ep*sup))    !resonant_added
     fon_res=auxres*ckl       ! resonant added
  ELSE                     ! resonant added
     skl=sin(kappaa*(usel-ea*su)+kappap*(up-ep*sup))    !resonant_added
     fon_res=auxres*skl        ! resonant added
  ENDIF
     ! **************************************************
     !  CALL check_deriv(cup,xp,d2,fon)
! **************************************************
  RETURN 
END FUNCTION fon_res


! *********************************************************************
! subroutine computing (1/d5-(r*r5)/|r5|^3)*(sin or cos)(critical angle) 
! 
DOUBLE PRECISION FUNCTION fsig_res(up)          ! resonant added
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface ----------------------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: d2,dist,drdon
  REAL(KIND=dkind) :: vsize,auxres
!  REAL(KIND=dkind) :: omega

! -------------------------------------------------------
  cup=cos(up)
  sup=sin(up)
!  CALL rot_matrix(elplr,om,omnod,g,zl,a)
! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  
  CALL prodmv(xp,Hp,yp)

  d2 = DOT_PRODUCT(xp-x,xp-x)
  dist = sqrt(d2)

! CORRECT FORMULA
  auxres = 2.d0*gmp*(1/dist - DOT_PRODUCT(x,xp)/vsize(xp)**3)* &
       & (1.d0-ea*cu)*(1.d0-ep*cup)

! DUMMY VALUE
!  auxres = 2.d0*gmp*1/dist*(1.d0-ea*cu)*(1.d0-ep*cup)
!  auxres = 2.d0*gmp*(- DOT_PRODUCT(x,xp)/vsize(xp)**3)* &
!       & (1.d0-ea*cu)*(1.d0-ep*cup)
!  auxres = (1.d0-ea*cu)*(1.d0-ep*cup)
!  auxres = 2.d0*gmp*(1/dist)*(1.d0-ea*cu)*(1.d0-ep*cup)

  IF (sin_flag.eq.0) THEN   
     ckl=cos(kappaa*(usel-ea*su)+kappap*(up-ep*sup))  
     fsig_res=auxres*ckl    
  ELSE                     
     skl=sin(kappaa*(usel-ea*su)+kappap*(up-ep*sup)) 
     fsig_res=auxres*skl   
  ENDIF

! DUMMY VALUE
!  omega = datan2(R_om(2,1),R_om(1,1))


  RETURN 
END FUNCTION fsig_res



! *********************************************************************
! subroutine computing S derivative (d(1/d)/dS)      
! This derivative is zero when the orbit of the planet is circular
DOUBLE PRECISION FUNCTION fs(up)     ! resonant added
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface ----------------------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: d2,d3,drda,drde
! -------------------------------------------------------
  cup=cos(up)
  sup=sin(up)

! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  CALL prodmv(xp,Hp,yp)

  d2 = DOT_PRODUCT(xp-x,xp-x)
  d3 = d2**(1.5d0)

! (-1/2)*der(D^2) = <Xpl-Xast, dR_on*R_i*R_om*Y>
!  drdon=DOT_PRODUCT(xp-x,donx)
  drda=DOT_PRODUCT(xp-x,dax)*(1.d0-ea*cu)

  drde = DOT_PRODUCT(xp-x,dex)*(1.d0-ea*cu)-cu*d2

! ======= DERIVATIVE with respect to S ===========
fs = -(gmp/d3)*(1.d0-ep*cup)*((beta**2/(ky*sqrt(aa)*ea))*drde + &
          &   2.d0*sqrt(a)/ky*drda)
 

! **************************************************
!  CALL check_deriv(cup,xp,d2,fon)
! **************************************************
  RETURN 
END FUNCTION fs


! subroutine computing S derivative (d(1/d)/dS-(r*r5)/|r5|^3)*(sin oppure cos)(critical angle)      
! This derivative is zero when the orbit of the planet is circular
DOUBLE PRECISION FUNCTION fs_res(up)     ! resonant added
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface ----------------------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: dist,d2,d3,drda,drde,dindde,dindda
  REAL(KIND=dkind) :: auxres,derauxres
  REAL(KIND=dkind) :: vsize
  REAL(KIND=dkind) :: derckl,derskl
! -------------------------------------------------------
  cup=cos(up)
  sup=sin(up)

! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  CALL prodmv(xp,Hp,yp)

  d2 = DOT_PRODUCT(xp-x,xp-x)
  d3 = d2**(1.5d0)
  dist = sqrt(d2)

! (-1/2)*der(D^2) = <Xpl-Xast, dR_on*R_i*R_om*Y>
!  drdon=DOT_PRODUCT(xp-x,donx)
  drda=DOT_PRODUCT(xp-x,dax)*(1.d0-ea*cu)

! d/da((r*r5)/|r5|^3*(1-e*cosu))
  dindda=DOT_PRODUCT(dax,xp)*(1.d0-ea*cu)/vsize(xp)**3

  drde = DOT_PRODUCT(xp-x,dex)*(1.d0-ea*cu)-cu*d2

! d/de((r*r5)/|r5|^3*(1-e*cosu))
  dindde = ( DOT_PRODUCT(dex,xp)*(1.d0-ea*cu) - &
       & DOT_PRODUCT(x,xp)*cu )/vsize(xp)**3


! CORRECT VALUE
  auxres = -2.d0*gmp*(1/dist - DOT_PRODUCT(x,xp)/vsize(xp)**3)* &
       & (1.d0-ea*cu)*(1.d0-ep*cup)

! ======= DERIVATIVE with respect to S ===========

  derauxres = -2.d0*gmp*(1.d0-ep*cup) * ( &
       & (beta**2/(ky*dsqrt(aa)*ea)) * (drde/d3 - dindde) + &
       & (2*sqrt(a)/ky) * (drda/d3 - dindda) ) * kappaa !(=dL/dS)   

  IF(sin_flag.eq.0) THEN   
     ckl=cos(kappaa*(usel-ea*su)+kappap*(up-ep*sup)) 
     derckl = sin(kappaa*(usel-ea*su)+kappap*(up-ep*sup))*kappaa*su
     fs_res= derauxres*ckl + auxres*derckl*(beta**2/(ky*dsqrt(aa)*ea)*kappaa)
  ELSE                     
     skl=sin(kappaa*(usel-ea*su)+kappap*(up-ep*sup))   
     derskl = -cos(kappaa*(usel-ea*su)+kappap*(up-ep*sup))*kappaa*su
     fs_res= derauxres*skl + auxres*derskl*(beta**2/(ky*dsqrt(aa)*ea)*kappaa)  
  ENDIF

! **************************************************
!  CALL check_deriv(cup,xp,d2,fon)
! **************************************************
  RETURN 
END FUNCTION fs_res









! ********************************************************************
! subroutine computing rotational matrix (with derivatives) 
! and coordinate vectors
SUBROUTINE rot_matrix(elpl,om,omnod,g,zl,a)
! ----------------- interface ------------------
  TYPE(orbit_elem),INTENT(IN) :: elpl !planet elems 
  REAL(KIND=dkind),INTENT(IN) :: om,omnod,g,zl,a
! ---------------- end interface ---------------
  REAL(KIND=dkind) :: omp,ip,onod,comp,somp,cip,sip,con,son,com,som
  REAL(KIND=dkind),DIMENSION(3,3) :: maux
! ----------------------------------------------
! planet elements
  apl = elpl%coord(1)
  ep = elpl%coord(2)
  omp = elpl%coord(5) 
  ip = elpl%coord(3) 
  betap = sqrt(1.d0-ep**2)
! Rotation matrix of angle omega for the planet
  comp = cos(omp)
  somp = sin(omp)
  R_omp = 0.d0
  R_omp(1,1) = comp
  R_omp(1,2) =-somp
  R_omp(2,1) = somp
  R_omp(2,2) = comp
  R_omp(3,3) = 1.d0
! Rotation matrix of angle Inclination for the planet
  cip = cos(ip)
  sip = sin(ip)
  R_ip = 0.d0
  R_ip(1,1) = 1.d0
  R_ip(2,2) = cip
  R_ip(2,3) =-sip
  R_ip(3,2) = sip
  R_ip(3,3) = cip
! asteroid elements
  aa = a
  ea = dsqrt(1.d0-(g/ky)**2/aa)
  beta = sqrt(1.d0-ea**2) 
! Rotation matrix of angle omega for the asteroid
  com = cos(om)
  som = sin(om)
  R_om = 0.d0
  R_om(1,1) = com
  R_om(1,2) =-som
  R_om(2,1) = som
  R_om(2,2) = com
  R_om(3,3) = 1.d0
! Derivative of the rotation matrix of angle omega for the asteroid
  dR_om = 0.d0
  dR_om(1,1) =-som
  dR_om(1,2) =-com
  dR_om(2,1) = com
  dR_om(2,2) =-som
 ! Compute cosine of inclination from Z/G for the asteroid
  ci = zl/g
  IF(ci.lt.1.d0)THEN 
     si = dsqrt(1.d0-ci**2) 
  ELSE 
     si = 0.d0 
  ENDIF
 ! Rotation matrix of angle Inclination for the asteroid
  R_i = 0.d0
  R_i(1,1) = 1.d0
  R_i(2,2) = ci
  R_i(2,3) =-si
  R_i(3,2) = si
  R_i(3,3) = ci
! Derivative of the rotation matrix of angle Inclination for the asteroid
  dR_i = 0.d0
  dR_i(2,2) =-si
  dR_i(2,3) =-ci
  dR_i(3,2) = ci
  dR_i(3,3) =-si
! Angle difference between asteroid Omnod and planet Omnod
  onod = omnod-elpl%coord(4)
  con = cos(onod)
  son = sin(onod)
! Rotation matrix of angle Omega.nodal for the asteroid
  R_on = 0.d0
  R_on(1,1) = con
  R_on(1,2) =-son
  R_on(2,1) = son
  R_on(2,2) = con
  R_on(3,3) = 1.d0
! Derivative of the rotation matrix of angle Omega.nodal for the asteroid
  dR_on = 0.d0
  dR_on(1,1) =-son
  dR_on(1,2) =-con
  dR_on(2,1) = con
  dR_on(2,2) =-son

!  CALL test_orbits(0.d0)

  CALL prodmm(maux,R_on,R_i)
  CALL prodmm(Hast,maux,R_om) 

  CALL prodmm(Hp,R_ip,R_omp)

  CALL prodmm(maux,R_on,dR_i)
  CALL prodmm(dHdI,maux,R_om)

  CALL prodmm(maux,R_on,R_i)
  CALL prodmm(dHdom,maux,dR_om)

  CALL prodmm(maux,dR_on,R_i)
  CALL prodmm(dHdon,maux,R_om)

END SUBROUTINE rot_matrix

! *******************************************************
SUBROUTINE test_orbits(u)
! --------------- interface -------------
  REAL(KIND=dkind),INTENT(IN) :: u
! ------------ end interface ------------  
  REAL(KIND=dkind),DIMENSION(3,3) :: maux
  REAL(KIND=dkind) :: s
  INTEGER :: i,j
! ---------------------------------
  cu=cos(u)
  su=sin(u)

! Asteroid coordinates in the orbital plane
  y(1)=aa*(cu-ea)
  y(2)=aa*su*beta
  y(3)=0.d0

! Asteroid coordinates in the rotated system
!  CALL prodmm(maux,R_on,R_i)
!  CALL prodmm(H2,maux,R_om) 
  CALL prodmv(x,Hast,y)

  WRITE(*,*)'u',u
  WRITE(*,*)'aa,ea,beta',aa,ea,beta
!  WRITE(*,*)'co,so',co,so
!  WRITE(*,*)'ci,si',ci,si
  WRITE(*,*)'ast position', x(1:3)

END SUBROUTINE test_orbits

! ==================================================================
! COMPUTING THE AVERAGING PERTURBING FUNCTION VALUE
! AT EACH STEP OF THE EVOLUTION TIME
! ==================================================================
! INTEGRATION DOMAIN: T = [0,DPIG] x [0,DPIG]                       
! y(1)=omega, y(2)=G, y(3)=Omega y(4)=Z y(5)=S, y(6)=sigma
SUBROUTINE secpert(numpla,elpl,om,G,omnod,zl,bigelle,sigma,ddd,eee,nnn) 
! **** controllare se passare n serve oppure no *****
! ***************************************************
! --------------- interface ---------------------
  INTEGER,INTENT(IN) :: numpla !planets number
  TYPE(orbit_elem),INTENT(IN) :: elpl !planet elems 
  REAL(KIND=dkind),INTENT(IN) :: om,omnod,G,zl,bigelle,sigma
  REAL(KIND=dkind),INTENT(OUT) :: ddd,eee !ddd=hamiltonian value
  INTEGER,INTENT(OUT) :: nnn
! -------------- end interface -------------------
! for dqags                                                         
  INTEGER :: limx,limx4 
  PARAMETER (limx=500,limx4=4*limx)
  INTEGER :: ier,limit,iwork(limx),lenw,last 
! function evaluations                                              
  INTEGER :: neval
  REAL(KIND=dkind) :: abserr,work(limx4) 
! output of dqags                                                 
  REAL(KIND=dkind) :: rm,rm_res_c,rm_res_s
  REAL(KIND=dkind) :: np,a
! ======================================================
  gmp = gm(numpla)
  a = (bigelle/ky)**2
  CALL rot_matrix(elpl,om,omnod,g,zl,a) 
! initialization
  ddd=0.d0 
  eee=0.d0 
  nnn=0 
! preparing dqags call
!  epsabs=1.d-9 
!  epsrel=1.d-7 
  limit=limx 
  lenw=limx4
!     ========= perturbing function ==========
! per il momento uso l'integratore dqags anziche' dqagp
!  CALL dqags(ff,0.d0,dpig,npts2,points,epsabs,epsrel,rm,abserr,  &
!       &        neval,ier,leniw,lenw,last,iwork,work)                     
!  ddd=ddd+rm/(dpig**2) 
!  eee=eee+abserr/(dpig**2) 
!  nnn=nnn+neval
  
  rm_res_c=0.d0
  rm_res_c=0.d0
  CALL dqags(ffd_pert,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)                               

! ********************************
!  rm =0.d0 ! DUMMY VALUE!!!!  

  IF (numpla.eq.n_res) THEN
!  write(*,*)'compute resonant term'
!  ider =5 ! forse non serve
!  res_flag=1

     sin_flag=0
     CALL dqags(ffd_pert_res,0.d0,dpig,epsabs,epsrel,rm_res_c,abserr,neval,ier, &
          &        limit,lenw,last,iwork,work)
     sin_flag= 1 
     CALL dqags(ffd_pert_res,0.d0,dpig,epsabs,epsrel,rm_res_s,abserr,neval,ier,&
          & limit,lenw,last,iwork,work)
     
     np = ky/elpl%coord(1)**(3.d0/2.d0)

! CORRECT FORMULA  
  ddd = ddd - ky**4/(2.d0*(kappaa*bigelle)**2) + np*kappap*bigelle - &
       & (rm_res_c*cos(sigma)+rm_res_s*sin(sigma) )/(dpig**2)

! GFGFG DUMMY VALUE!!!
!     ddd = -ky**4/(2.d0*(kappaa*bigelle)**2) + np*kappap*bigelle
!  ddd = ddd - (rm_res_c*cos(sigma)+rm_res_s*sin(sigma) )/(dpig**2)

!     write(*,133)'secpert: ky,kappaa,kappap,bigelle=',ky,kappaa,kappap,bigelle
!133  format(a33,f10.5,2x,(2i2,2x),f10.5)

!  ddd = -(rm_res_c*cos(sigma)+rm_res_s*sin(sigma) )/(dpig**2)
!    ddd = 0.d0   ! to take only the int(1/d) part
!  ddd =  -rm_res_c
!  ddd =  -rm_res_s
  ENDIF

  
! HINT! temporary switch off 
  ddd=ddd+rm/(dpig**2) 
  eee=eee+abserr/(dpig**2) 
  nnn=nnn+neval
  
!  write(*,*)'ddd=',ddd


END SUBROUTINE secpert

! ********************************************************************     
! subroutine computing unidim.average of the unidim.average         
DOUBLE PRECISION FUNCTION ffd_pert(u) 
! --------------- interface ---------------------
  REAL(KIND=dkind) :: u !eccentric anomaly of the asteroid
! ------------ end interface ------------
! for dqagsc                                                        
  INTEGER :: limx,limx4 
  PARAMETER (limx=500,limx4=4*limx) 
  INTEGER :: neval,ier,limit,iwork(limx),lenw,last 
  REAL(KIND=dkind) :: abserr,work(limx4),risult 
! ----------------------------------------------
  cu=cos(u)
  su=sin(u)
! Asteroid coord
  y(1) = aa*(cu-ea)
  y(2) = aa*su*beta
  y(3) = 0.d0
  CALL prodmv(x,Hast,y)
! ------------ preparing dqagsc call ------------
!  epsabs=1.d-8 
!  epsrel=1.d-5 
  limit=limx 
  lenw=limx4
! -------------------------------------------------------
  CALL dqagsc(fun_pert,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
          &        ier,limit,lenw,last,iwork,work) 
!     if(ier.gt.1)write(*,*)'pert fun ',ier,u

  ffd_pert=risult  

  RETURN 
END FUNCTION ffd_pert

! ********************************************************************     
! subroutine computing unidim.average of the unidim.average         
DOUBLE PRECISION FUNCTION ffd_pert_res(u) 
! --------------- interface ---------------------
  REAL(KIND=dkind) :: u !eccentric anomaly of the asteroid
! ------------ end interface ------------
! for dqagsc                                                        
  INTEGER :: limx,limx4 
  PARAMETER (limx=500,limx4=4*limx) 
  INTEGER :: neval,ier,limit,iwork(limx),lenw,last 
  REAL(KIND=dkind) :: abserr,work(limx4),risult 
! ----------------------------------------------

  cu=cos(u)
  su=sin(u)
  usel = u
! Asteroid coord
  y(1) = aa*(cu-ea)
  y(2) = aa*su*beta
  y(3) = 0.d0
  CALL prodmv(x,Hast,y)

!  write(*,*)'y=',y
!  write(*,*)''
!  write(*,*)'x=',x


! ------------ preparing dqagsc call ------------
!  epsabs=1.d-8 
!  epsrel=1.d-5 
  limit=limx 
  lenw=limx4
! -------------------------------------------------------

! ****** GFG 27/6/2016 ******

  risult = 0.d0
  abserr = 0.d0
  neval = 0
!  write(*,*)'x in ffd_pert_res:',x

  CALL dqagsc(fsig_res,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
       &        ier,limit,lenw,last,iwork,work)
  
   ffd_pert_res=risult  

  RETURN 
END FUNCTION ffd_pert_res


! ******************************************************
! subroutine computing perturbing function
DOUBLE PRECISION FUNCTION fun_pert(up)
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface ------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: d2
! --------------------------------------
  cup=cos(up)
  sup=sin(up)
! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  CALL prodmv(xp,Hp,yp)

  d2 = DOT_PRODUCT(xp-x,xp-x)

! ----- perturbing function -----
  fun_pert = -(gmp/sqrt(d2))*(1.d0-ea*cu)*(1.d0-ep*cup)

  RETURN 
END FUNCTION fun_pert

!! ******************************************************
!! subroutine computing perturbing function
!DOUBLE PRECISION FUNCTION fun_pert_res(up)
!! --------------- interface ---------------------
!  REAL(KIND=dkind) :: up  !planet eccentric anomaly
!! ------------ end interface ------------
!  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
!  REAL(KIND=dkind) :: d2
!! --------------------------------------
!  cup=cos(up)
!  sup=sin(up)
!! Planet coordinates
!  yp(1)=apl*(cup-ep)
!  yp(2)=apl*sup*betap
!  yp(3)=0.d0
!  CALL prodmv(xp,Hp,yp)
!
!  d2 = DOT_PRODUCT(xp-x,xp-x)
!
!!! ----- perturbing function -----
 ! fun_pert_res = -(gmp/sqrt(d2))*(1.d0-ea*cu)*(1.d0-ep*cup)
!
 ! RETURN 
!END FUNCTION fun_pert_res

! ==================================================================
! COMPUTING SECOND DERIVATIVES OF R
! ==================================================================
! INTEGRATION DOMAIN: T = [0,DPIG] x [0,DPIG]                       
SUBROUTINE derrhs(n,elpl,om,omnod,g,zl,a,ddd_ii,eee_ii,nnn_ii) 
! --------------- interface ---------------------
  INTEGER,INTENT(IN) :: n !planets number
  TYPE(orbit_elem),INTENT(IN) :: elpl !planet elems 
  REAL(KIND=dkind),INTENT(IN) :: om,omnod,g,zl,a
  REAL(KIND=dkind),INTENT(OUT) :: ddd_ii(4),eee_ii(4) 
  INTEGER,INTENT(OUT) :: nnn_ii(4)
! -------------- end interface -------------------
! for dqags                                                         
  INTEGER :: limx,limx4 
  PARAMETER (limx=500,limx4=4*limx)
  INTEGER :: ier,limit,iwork(limx),lenw,last 
! function evaluations                                              
  INTEGER :: neval 
  REAL(KIND=dkind) :: abserr,work(limx4) 
! output of dqags                                                   
  REAL(KIND=dkind) :: rm
! loop indexes                                                      
  INTEGER :: k,i
! for derivatives wrt I
  REAL(KIND=dkind) :: ddR_I(3,3),maux(3,3)
! ======================================================
  gmp = gm(n)
  CALL rot_matrix(elpl,om,omnod,g,zl,a) 
! matrix for 2nd-derivative of ast coord wrt Inc
  ddR_i = 0.d0
  ddR_i(2:3,2:3) = -R_i(2:3,2:3)
  CALL prodmm(maux,R_on,ddR_i)
  CALL prodmm(ddHdI,maux,R_om)
! initialization
  ddd_ii=0.d0 
  eee_ii=0.d0 
  nnn_ii=0 
! preparing dqags call
!  epsabs=1.d-8 
!  epsrel=1.d-5 
  limit=limx 
  lenw=limx4 
! ---------------------------------------------
! first computing integral with respect to u
! ---------------------------------------------
! 2nd-derivative with respect to G
  ider=5
  CALL dqags(ffd_ii,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)
  ddd_ii(1)=ddd_ii(1)+rm/(dpig**2)
  eee_ii(1)=eee_ii(1)+abserr/(dpig**2)
  nnn_ii(1)=nnn_ii(1)+neval
! 2nd-derivative with respect to omega                                  
  ider=6 
  CALL dqags(ffd_ii,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)                               
  ddd_ii(2)=ddd_ii(2)+rm/(dpig**2) 
  eee_ii(2)=eee_ii(2)+abserr/(dpig**2) 
  nnn_ii(2)=nnn_ii(2)+neval
! 2nd-derivative with respect to Z
  ider=7
  CALL dqags(ffd_ii,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)
  ddd_ii(3)=ddd_ii(3)+rm/(dpig**2)
  eee_ii(3)=eee_ii(3)+abserr/(dpig**2)
  nnn_ii(3)=nnn_ii(3)+neval
! 2nd-derivative with respect to onod
  ider=8
  CALL dqags(ffd_ii,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)
  ddd_ii(4)=ddd_ii(4)+rm/(dpig**2)
  eee_ii(4)=eee_ii(4)+abserr/(dpig**2)
  nnn_ii(4)=nnn_ii(4)+neval

END SUBROUTINE derrhs

! subroutine computing unidim.average of the unidim.average         
DOUBLE PRECISION FUNCTION ffd_ii(u) 
! --------------- interface ---------------------
  REAL(KIND=dkind) :: u !asteroid eccentric anomaly
! ------------ end interface ------------
  REAL(KIND=dkind) :: dy(3)!,dR(3,3),ddH(3,3),maux(3,3)
! for dqagsc
  INTEGER :: limx,limx4 
  PARAMETER (limx=500,limx4=4*limx) 
  INTEGER :: neval,ier,limit,iwork(limx),lenw,last 
  REAL(KIND=dkind) :: abserr,work(limx4),risult 
! ------------------------------------------------------------
  cu=cos(u)
  su=sin(u)
! Asteroid coord
  y(1) = aa*(cu-ea)
  y(2) = aa*su*beta
  y(3) = 0.d0
  CALL prodmv(x,Hast,y)
! Asteroid coord vectors : 1th-derivatives w.r.t. ecc,Inc,omega,Omnod
  dy(1) =-aa
  dy(2) =-aa*su*ea/beta
  dy(3) = 0.d0
  CALL prodmv(dex,Hast,dy)
  CALL prodmv(dIx,dHdI,y)
  CALL prodmv(domx,dHdom,y) 
  CALL prodmv(donx,dHdon,y)
! Asteroid coord vectors : 2th-derivatives w.r.t. ecc,Inc,omega,Omnod
  CALL prodmv(ddeIx,dHdI,dy)
  dy(1) = 0.d0
  dy(2) =-aa*su/beta*(1.d0+(ea/beta)**2)
  dy(3) = 0.d0
  CALL prodmv(ddex,Hast,dy)
  CALL prodmv(ddIx,ddHdI,y)
  ddomx = -x
  ddonx(1:2) = -x(1:2)
  ddonx(3) = 0.d0
! ------------ preparing dqagsc call ------------
!  epsabs=1.d-8 
!  epsrel=1.d-5 
  limit=limx 
  lenw=limx4
!  ------------ selecting derivative  ------------
  IF(ider.eq.5)THEN 
     call dqagsc(fg_ii,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
          &        ier,limit,lenw,last,iwork,work) 
!     if(ier.gt.1)write(*,*)'d/dom ',ier,u
  ELSEIF(ider.eq.6)THEN
     call dqagsc(fom_ii,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
          &        ier,limit,lenw,last,iwork,work)
! if(ier.gt.1)write(*,*)'d/dom ',ier,u 
  ELSEIF(ider.eq.7)THEN
     call dqagsc(fz_ii,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
          &        ier,limit,lenw,last,iwork,work)
! if(ier.gt.1)write(*,*)'d/dZ ',ier,u   
  ELSEIF(ider.eq.8)THEN
     call dqagsc(fon_ii,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
          &        ier,limit,lenw,last,iwork,work)
! if(ier.gt.1)write(*,*)'d/donod ',ier,u
  ELSE 
     WRITE(*,*)' ffd_ii: ider=',ider 
     STOP 
  ENDIF
 
  ffd_ii=risult
  
  RETURN 
END FUNCTION ffd_ii


! ********************************************************************* 
! subroutine computing G derivative (ddR/ddG)
DOUBLE PRECISION FUNCTION fg_ii(up)
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly 
! ------------ end interface ----------------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: d2,d3,de,dI,dde,ddI,dfdg,ddfdg
  REAL(KIND=dkind),DIMENSION(3) :: dgx,ddgx
! --------------------------------------------------
  cup=cos(up)
  sup=sin(up)
! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  CALL prodmv(xp,Hp,yp)

  d2 = DOT_PRODUCT(xp-x,xp-x)
  d3 = d2**(1.5d0)

  de = -beta/(ky*dsqrt(aa)*ea)   !de/dG
  dI = ci/(si*ky*dsqrt(aa)*beta) !dI/dG
  dgx = dex*de+dIx*dI

  dde = (de**2)*(-ea/beta**2-1.d0/ea)            !dde/ddG   
  ddI = (dI**2)*(-1.d0/(ci*si))+de*dI*ea/beta**2 !ddI/ddG
  ddgx = ddex*(de**2) + dex*dde + ddIx*(dI**2) + dIx*ddI + 2.d0*ddeIx*dI*de

  dfdg = DOT_PRODUCT(xp-x,dgx)
  ddfdg = 3.d0/d2*(dfdg**2) + DOT_PRODUCT(xp-x,ddgx) - DOT_PRODUCT(dgx,dgx)

! ========= SECOND DERIVATIVE with respect to G ============
  fg_ii = (gmp/d3)*(1.d0-ep*cup)*((1.d0-ea*cu)*ddfdg-2.d0*cu*de*dfdg-cu*dde*d2)

! **************************************************
!  CALL check_secderiv(cup,xp,d3,fg_ii)
! **************************************************

  RETURN
END FUNCTION fg_ii

! ********************************************************************* 
! subroutine computing omega derivative (ddR/ddom)
DOUBLE PRECISION FUNCTION fom_ii(up)
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface ------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: d2,d3,dfdom,ddfdom
! --------------------------------------
! check on derivate
  REAL(KIND=dkind) :: maux(3,3),vaux(3),dm(3,3),ddHdom(3,3)
! --------------------------------------
  cup=cos(up)
  sup=sin(up)
! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  CALL prodmv(xp,Hp,yp)

  d2 = DOT_PRODUCT(xp-x,xp-x)
  d3 = d2**(1.5d0)

  dfdom = DOT_PRODUCT(xp-x,domx)
  ddfdom = 3.d0/d2*(dfdom**2) + DOT_PRODUCT(xp-x,ddomx) - DOT_PRODUCT(domx,domx)
  
! ========= SECOND DERIVATIVE with respect to omega ============
  fom_ii = (gmp/d3)*(1.d0-ea*cu)*(1.d0-ep*cup)*ddfdom

! **************************************************
!  CALL check_secderiv(cup,xp,d3,fom_ii)
! **************************************************
  RETURN
END FUNCTION fom_ii


! ********************************************************************
! subroutine computing Z derivative (ddR/ddZ) 
DOUBLE PRECISION FUNCTION fz_ii(up) 
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface ------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: d2,d3,dfdI,ddfdI,dI
! ----------------------------------------
  cup=cos(up)
  sup=sin(up)
! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  CALL prodmv(xp,Hp,yp)

  d2=DOT_PRODUCT(xp-x,xp-x)
  d3=d2**(1.5d0)

  dfdI = DOT_PRODUCT(xp-x,dIx)
  ddfdI = 3.d0/d2*(dfdI**2) + DOT_PRODUCT(xp-x,ddIx) - DOT_PRODUCT(dIx,dIx)

  dI = -1.d0/(ky*sqrt(aa)*beta*si)

! ========== SECOND DERIVATIVE with respect to Z ==================
  fz_ii = (gmp/d3)*(1.d0-ea*cu)*(1.d0-ep*cup)*(dI**2)*(ddfdI - ci/si*dfdI)

! **************************************************
!  CALL check_secderiv(cup,xp,d3,fz_ii)
! **************************************************
  RETURN
END FUNCTION fz_ii


! ********************************************************************
! subroutine computing Omnod derivative (ddR/ddOmnod) 
DOUBLE PRECISION FUNCTION fon_ii(up) 
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface -------------------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: d2,d3,dfdon,ddfdon
! ----------------------------------------------------
  cup=cos(up)
  sup=sin(up)
! Planet coordinates
  yp(1)=apl*(cup-ep)
  yp(2)=apl*sup*betap
  yp(3)=0.d0
  CALL prodmv(xp,Hp,yp)

  d2 = DOT_PRODUCT(xp-x,xp-x)
  d3 = d2**(1.5d0)

  dfdon = DOT_PRODUCT(xp-x,donx)
  ddfdon = 3.d0/d2*(dfdon**2) + DOT_PRODUCT(xp-x,ddonx) - DOT_PRODUCT(donx,donx) 

! ========== SECOND DERIVATIVE with respect to Omnod ==================
  fon_ii = (gmp/d3)*(1.d0-ea*cu)*(1.d0-ep*cup)*ddfdon

! **************************************************
!  CALL check_secderiv(cup,xp,d3,fon_ii)
! **************************************************
  RETURN
END FUNCTION fon_ii

! ********************************************************
! CHECK on FIRST derivatives of R w.r.t. G,om,Z,Omnod
! ********************************************************
SUBROUTINE check_deriv(cup,xp,d2,der)
! --------------- interface ------------
  REAL(KIND=dkind) :: cup,xp(3),d2,der
! -------- end interface ---------------
  REAL(KIND=dkind) :: inc,f0,f0_inc
  REAL(KIND=dkind) :: g,g_inc,zl,zl_inc,co,so,co_inc,so_inc
  REAL(KIND=dkind) :: ea_inc,ci_inc,si_inc,d2_inc
  REAL(KIND=dkind),DIMENSION(3,3) :: R_inc,maux,Hast_inc
  REAL(KIND=dkind),DIMENSION(3) :: y_inc,x_inc
! -------------------------------------

!  write(*,*)'ider=',ider

  f0 = (gmp/dsqrt(d2))*(1.d0-ep*cup)*(1.d0-ea*cu)
  inc = 1.d-6

  IF(ider.eq.1)THEN
     write(*,*)'**************************'
     write(*,*)'CHECK derivative w.r.t. G'
     g = ky*sqrt(aa)*beta
     g_inc = g+inc
     ea_inc = dsqrt(1.d0-(g_inc/ky)**2/aa)
     ci_inc = g*ci/g_inc
     si_inc = dsqrt(1.d0-ci_inc**2)
     R_inc = 0.d0
     R_inc(1,1) = 1.d0
     R_inc(2,2) = ci_inc
     R_inc(2,3) =-si_inc
     R_inc(3,2) = si_inc
     R_inc(3,3) = ci_inc
     CALL prodmm(maux,R_on,R_inc)
     CALL prodmm(Hast_inc,maux,R_om)
     y_inc(1) = aa*(cu-ea_inc)
     y_inc(2) = aa*su*dsqrt(1.d0-ea_inc**2)
     y_inc(3) = 0.d0
     CALL prodmv(x_inc,Hast_inc,y_inc)
     d2_inc = DOT_PRODUCT(xp-x_inc,xp-x_inc)
     f0_inc = (gmp/dsqrt(d2_inc))*(1.d0-ep*cup)*(1.d0-ea_inc*cu)

  ELSEIF(ider.eq.2)THEN
     write(*,*)'**************************'
     write(*,*)'CHECK derivative w.r.t. omega'
     co=R_om(1,1)
     so=R_om(2,1)
     co_inc = co*cos(inc)-so*sin(inc)
     so_inc = so*cos(inc)+co*sin(inc)
     R_inc = 0.d0
     R_inc(1,1) = co_inc
     R_inc(1,2) =-so_inc
     R_inc(2,1) = so_inc
     R_inc(2,2) = co_inc
     R_inc(3,3) = 1.d0
     CALL prodmm(maux,R_on,R_i)
     CALL prodmm(Hast_inc,maux,R_inc) 
     CALL prodmv(x_inc,Hast_inc,y)
     d2_inc = DOT_PRODUCT(xp-x_inc,xp-x_inc)
     f0_inc = (gmp/dsqrt(d2_inc))*(1.d0-ep*cup)*(1.d0-ea*cu)

  ELSEIF(ider.eq.3)THEN
     write(*,*)'**************************'
     write(*,*)'CHECK derivative w.r.t. Z'
     g = ky*sqrt(aa)*beta
     zl=g*ci
     zl_inc = zl+inc
     ci_inc = zl_inc/g
     si_inc = sqrt(1.d0-ci_inc**2)
     R_inc = 0.d0
     R_inc(1,1) = 1.d0
     R_inc(2,2) = ci_inc
     R_inc(2,3) =-si_inc
     R_inc(3,2) = si_inc
     R_inc(3,3) = ci_inc
     CALL prodmm(maux,R_on,R_inc)
     CALL prodmm(Hast_inc,maux,R_om) 
     CALL prodmv(x_inc,Hast_inc,y)
     d2_inc = DOT_PRODUCT(xp-x_inc,xp-x_inc)
     f0_inc = (gmp/sqrt(d2_inc))*(1.d0-ep*cup)*(1.d0-ea*cu)

  ELSEIF(ider.eq.4)THEN
     write(*,*)'**************************'
     write(*,*)'CHECK derivative w.r.t. Omnod'
     co=R_on(1,1)
     so=R_on(2,1)
     co_inc = co*cos(inc)-so*sin(inc)
     so_inc = so*cos(inc)+co*sin(inc)
     R_inc = 0.d0
     R_inc(1,1) = co_inc
     R_inc(1,2) =-so_inc
     R_inc(2,1) = so_inc
     R_inc(2,2) = co_inc
     R_inc(3,3) = 1.d0
     CALL prodmm(maux,R_inc,R_i)
     CALL prodmm(Hast_inc,maux,R_om) 
     CALL prodmv(x_inc,Hast_inc,y)
     d2_inc = DOT_PRODUCT(xp-x_inc,xp-x_inc)
     f0_inc = (gmp/dsqrt(d2_inc))*(1.d0-ep*cup)*(1.d0-ea*cu)

  ELSE
     write(*,*)'check_deriv : no derivative selected ',ider
  ENDIF
 
  write(*,*)'incremental ratio:',(f0_inc-f0)/inc
  write(*,*)'derivative:       ',der
  write(*,*)'********************************************'


  STOP

END SUBROUTINE check_deriv

! ********************************************************
! CHECK on SECOND derivatives of R w.r.t. G,om,Z,Omnod
! ********************************************************
SUBROUTINE check_secderiv(cup,xp,d3,der)
! --------------- interface -----------------------------
  REAL(KIND=dkind) :: cup,xp(3),d3,der
! -------- end interface --------------------------------
  REAL(KIND=dkind) :: inc,df,df_inc,d3_inc,d2_inc,d2
  REAL(KIND=dkind) :: g,g_inc,zl,zl_inc,co,so,co_inc,so_inc
  REAL(KIND=dkind) :: ea_inc,beta_inc,ci_inc,si_inc,de,de_inc,dI,dI_inc
  REAL(KIND=dkind),DIMENSION(3,3) :: R_inc,maux,Hast_inc,dR_inc,dH_inc
  REAL(KIND=dkind),DIMENSION(3) :: y_inc,x_inc,dx_inc,dy_inc,dx,dex_inc,dIx_inc
! -------------------------------------------------------

  write(*,*)'ider=',ider

  inc = 1.d-5

  IF(ider.eq.5)THEN
     write(*,*)'**************************'
     write(*,*)'CHECK second derivative w.r.t. G'
     g = ky*sqrt(aa)*beta
     dI = -1.d0/(g*si)
     de = -beta**2/(g*ea)
     dI = ci/(si*g)
     dx = dex*de + dIx*dI
     d2 = DOT_PRODUCT(xp-x,xp-x)
     df = (gmp/d3)*(1.d0-ep*cup)*((1.d0-ea*cu)*DOT_PRODUCT(xp-x,dx)-cu*de*d2)
     g_inc = g+inc
     ea_inc = dsqrt(1.d0-(g_inc/ky)**2/aa)
     beta_inc = dsqrt(1.d0-ea_inc**2)
     y_inc(1) = aa*(cu-ea_inc)
     y_inc(2) = aa*su*beta_inc
     y_inc(3) = 0.d0
     ci_inc = g*ci/g_inc
     si_inc = dsqrt(1.d0-ci_inc**2)
     R_inc = 0.d0
     R_inc(1,1) = 1.d0
     R_inc(2,2) = ci_inc
     R_inc(2,3) =-si_inc
     R_inc(3,2) = si_inc
     R_inc(3,3) = ci_inc
     CALL prodmm(maux,R_on,R_inc)
     CALL prodmm(Hast_inc,maux,R_om)
     CALL prodmv(x_inc,Hast_inc,y_inc)
     dy_inc(1) =-aa
     dy_inc(2) =-aa*su*(ea_inc/beta_inc)
     dy_inc(3) = 0.d0
     CALL prodmv(dex_inc,Hast_inc,dy_inc)
     de_inc = -beta_inc**2/(g_inc*ea_inc)
     dR_inc = 0.d0
     dR_inc(2,2) =-si_inc
     dR_inc(2,3) =-ci_inc
     dR_inc(3,2) = ci_inc
     dR_inc(3,3) =-si_inc
     CALL prodmm(maux,R_on,dR_inc)
     CALL prodmm(dH_inc,maux,R_om)
     CALL prodmv(dIx_inc,dH_inc,y_inc)
     dI_inc = ci_inc/(si_inc*g_inc)
     dx_inc = dex_inc*de_inc + dIx_inc*dI_inc
     d2_inc = DOT_PRODUCT(xp-x_inc,xp-x_inc)
     d3_inc = d2_inc**(1.5d0)
     df_inc = (gmp/d3_inc)*(1.d0-ep*cup)*((1.d0-ea_inc*cu)* &
          & DOT_PRODUCT(xp-x_inc,dx_inc) - cu*de_inc*d2_inc)

  ELSEIF(ider.eq.6)THEN
     write(*,*)'**************************'
     write(*,*)'CHECK second derivative w.r.t. omega'
     df = (gmp/d3)*(1.d0-ep*cup)*(1.d0-ea*cu)*DOT_PRODUCT(xp-x,domx)
     co = R_om(1,1)
     so = R_om(2,1)
     co_inc = co*cos(inc)-so*sin(inc)
     so_inc = so*cos(inc)+co*sin(inc)
     R_inc = 0.d0
     R_inc(1,1) = co_inc
     R_inc(1,2) =-so_inc
     R_inc(2,1) = so_inc
     R_inc(2,2) = co_inc
     R_inc(3,3) = 1.d0
     CALL prodmm(maux,R_on,R_i)
     CALL prodmm(Hast_inc,maux,R_inc)
     CALL prodmv(x_inc,Hast_inc,y)
     dR_inc = 0.d0
     dR_inc(1,1) =-so_inc
     dR_inc(1,2) =-co_inc
     dR_inc(2,1) = co_inc
     dR_inc(2,2) =-so_inc
     CALL prodmm(maux,R_on,R_i)
     CALL prodmm(dH_inc,maux,dR_inc)
     CALL prodmv(dx_inc,dH_inc,y)
     d3_inc = DOT_PRODUCT(xp-x_inc,xp-x_inc)**(1.5d0)
     df_inc = (gmp/d3_inc)*(1.d0-ea*cu)*(1.d0-ep*cup)* &
          & DOT_PRODUCT(xp-x_inc,dx_inc)

  ELSEIF(ider.eq.7)THEN
     write(*,*)'**************************'
     write(*,*)'CHECK second derivative w.r.t. Z'
     g = ky*sqrt(aa)*beta
     dI = -1.d0/(g*si)
     df = (gmp/d3)*(1.d0-ep*cup)*(1.d0-ea*cu)*DOT_PRODUCT(xp-x,dIx)*dI
     zl=g*ci
     zl_inc = zl+inc
     ci_inc = zl_inc/g
     si_inc = sqrt(1.d0-ci_inc**2)
     R_inc = 0.d0
     R_inc(1,1) = 1.d0
     R_inc(2,2) = ci_inc
     R_inc(2,3) =-si_inc
     R_inc(3,2) = si_inc
     R_inc(3,3) = ci_inc
     CALL prodmm(maux,R_on,R_inc)
     CALL prodmm(Hast_inc,maux,R_om) 
     CALL prodmv(x_inc,Hast_inc,y)     
     dR_inc = 0.d0
     dR_inc(2,2) =-si_inc
     dR_inc(2,3) =-ci_inc
     dR_inc(3,2) = ci_inc
     dR_inc(3,3) =-si_inc
     CALL prodmm(maux,R_on,dR_inc)
     CALL prodmm(dH_inc,maux,R_om)
     CALL prodmv(dx_inc,dH_inc,y)
     d3_inc = DOT_PRODUCT(xp-x_inc,xp-x_inc)**(1.5d0)
     dI_inc = -1.d0/(g*si_inc)
     df_inc = (gmp/d3_inc)*(1.d0-ea*cu)*(1.d0-ep*cup)* &
          & DOT_PRODUCT(xp-x_inc,dx_inc)*dI_inc

  ELSEIF(ider.eq.8)THEN
     write(*,*)'**************************'
     write(*,*)'CHECK second derivative w.r.t. Omnod'
     df = (gmp/d3)*(1.d0-ep*cup)*(1.d0-ea*cu)*DOT_PRODUCT(xp-x,donx)
     co = R_on(1,1)
     so = R_on(2,1)
     co_inc = co*cos(inc)-so*sin(inc)
     so_inc = so*cos(inc)+co*sin(inc)
     R_inc = 0.d0
     R_inc(1,1) = co_inc
     R_inc(1,2) =-so_inc
     R_inc(2,1) = so_inc
     R_inc(2,2) = co_inc
     R_inc(3,3) = 1.d0
     CALL prodmm(maux,R_inc,R_i)
     CALL prodmm(Hast_inc,maux,R_om)
     CALL prodmv(x_inc,Hast_inc,y)
     dR_inc = 0.d0
     dR_inc(1,1) =-so_inc
     dR_inc(1,2) =-co_inc
     dR_inc(2,1) = co_inc
     dR_inc(2,2) =-so_inc
     CALL prodmm(maux,dR_inc,R_i)
     CALL prodmm(dH_inc,maux,R_om)
     CALL prodmv(dx_inc,dH_inc,y)
     d3_inc = DOT_PRODUCT(xp-x_inc,xp-x_inc)**(1.5d0)
     df_inc = (gmp/d3_inc)*(1.d0-ea*cu)*(1.d0-ep*cup)* &
          & DOT_PRODUCT(xp-x_inc,dx_inc)
     
  ELSE
     write(*,*)'check_secder : no derivative selected ',ider
  ENDIF

  write(*,*)'incremental ratio:',(df_inc-df)/inc
  write(*,*)'derivative:       ',der
  write(*,*)'********************************************'

  STOP

END SUBROUTINE check_secderiv


END MODULE right_hand_side
