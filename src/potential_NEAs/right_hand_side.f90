MODULE right_hand_side
  USE fund_const
  USE orbit_elements
  IMPLICIT NONE
  PRIVATE
! asteroid/planet orbital elements
  REAL(KIND=dkind) :: aa,ea,beta,ci,si,apl,ep,betap,gmp
! rotational matrix
  REAL(KIND=dkind),DIMENSION(3,3) :: R_omp,R_ip,Hp,R_om,R_i,R_on,Hast
  REAL(KIND=dkind),DIMENSION(3,3) :: dR_om,dR_i,dR_on,dHdom,dHdI,dHdon
  REAL(KIND=dkind),DIMENSION(3,3) :: ddHdI
! coord vectors and derivatives
  REAL(KIND=dkind) :: x(3),y(3),cu,su
  REAL(KIND=dkind),DIMENSION(3) :: dex,dIx,domx,donx
  REAL(KIND=dkind),DIMENSION(3) :: ddex,ddeIx,ddIx,ddomx,ddonx
! for selecting derivatives
  INTEGER :: ider
! planet data                                                       
!  INCLUDE 'pldata.h90'
! ap = semimajor axis of the planet, ky = Gauss constant 
! gm = (masse pianeti)/(massa Sole)
       INTEGER nplax,nplax2,npl,inpl,ioupl,ndum
       parameter (nplax=10)
       parameter (nplax2=20)
       DOUBLE PRECISION ap(nplax),gm(nplax)
       DOUBLE PRECISION ky,bigg 
!       COMMON/pldata/gm,ap,ky,ndum,npl,inpl,ioupl,bigg



! routines
  PUBLIC :: rhs2,secpert,pladat!derrhs

! variables
  PUBLIC :: inpl,ioupl,ky,ap

CONTAINS
! =================================================================
! RIGHT SIDE COMPUTATION                                       
! INTEGRATION DOMAIN: T = [0,DPIG] x [0,DPIG]                       
SUBROUTINE rhs2(n,elpl,om,omnod,g,zl,a,ddd,eee,nnn) 
! --------------- interface ---------------------
  INTEGER,INTENT(IN) :: n !planets number
  TYPE(orbit_elem),INTENT(IN) :: elpl !planet elems 
  REAL(KIND=dkind),INTENT(IN) :: om,omnod,g,zl,a
  REAL(KIND=dkind),INTENT(OUT) :: ddd(4),eee(4) 
  INTEGER,INTENT(OUT) :: nnn(4)
! -------------- end interface -------------------
! for dqags                                                         
  INTEGER :: limx,limx4 
  PARAMETER (limx=500,limx4=4*limx)
  INTEGER :: ier,limit,iwork(limx),lenw,last 
! function evaluations                                              
  INTEGER :: neval 
  REAL(KIND=dkind) :: epsabs,epsrel,abserr,work(limx4) 
! output of dqags                                                   
  REAL(KIND=dkind) :: rm
! loop indexes                                                      
  INTEGER :: k,i
! ======================================================
  gmp = gm(n)
! compute rotational matrix and Cartesian coords for planet and asteroid
  CALL rot_matrix(elpl,om,omnod,g,zl,a)
! ---- initialization -----
  DO k=1,4
     ddd(k)=0.d0 
     eee(k)=0.d0 
     nnn(k)=0 
  ENDDO
! ---- preparing dqags call ---
  epsabs=1.d-8 
  epsrel=1.d-5 
  limit=limx 
  lenw=limx4 
! ********************************************************************     
!                    COMPUTING DERIVATIVES 
! ********************************************************************     
! first computing integral with respect to u

! derivative with respect to G
  ider=1 
  CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)                               
  ddd(1)=ddd(1)+rm/(dpig**2) 
  eee(1)=eee(1)+abserr/(dpig**2) 
  nnn(1)=nnn(1)+neval 
  
! derivative with respect to omega                                  
  ider=2 
  CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)                               
  ddd(2)=ddd(2)+rm/(dpig**2) 
  eee(2)=eee(2)+abserr/(dpig**2) 
  nnn(2)=nnn(2)+neval 
                                                                        
! derivative with respect to Z                                      
  ider=3 
  CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)                               
  ddd(3)=ddd(3)+rm/(dpig**2) 
  eee(3)=eee(3)+abserr/(dpig**2) 
  nnn(3)=nnn(3)+neval 

! derivative with respect to onod
  ider=4 
  CALL dqags(ffd,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)                               
  ddd(4)=ddd(4)+rm/(dpig**2) 
  eee(4)=eee(4)+abserr/(dpig**2) 
  nnn(4)=nnn(4)+neval

END SUBROUTINE rhs2
                                                                        
! ********************************************************************     
! subroutine computing unidim.average of the unidim.average         
DOUBLE PRECISION FUNCTION ffd(u) 
! --------------- interface ---------------------
  REAL(KIND=dkind) :: u !eccentric anomaly of the asteroid
! ------------ end interface ----------------------
  REAL(KIND=dkind) :: dy(3)
! for dqags                                                         
  INTEGER :: limx,limx4 
  PARAMETER (limx=500,limx4=4*limx) 
  INTEGER :: neval,ier,limit,iwork(limx),lenw,last 
  REAL(KIND=dkind) :: epsabs,epsrel,abserr,work(limx4),result 
! --------------------------------------------------
  cu=cos(u)
  su=sin(u)
! Asteroid coordinates
  y(1) = aa*(cu-ea)
  y(2) = aa*su*beta
  y(3) = 0.d0
  CALL prodmv(x,Hast,y)
! 1th-derivatives w.r.t. ecc,Inc,omega,Omnod for asteroid coord
  dy(1) =-aa
  dy(2) =-aa*su*ea/beta
  dy(3) = 0.d0
  CALL prodmv(dex,Hast,dy)
  CALL prodmv(dIx,dHdI,y)
  CALL prodmv(domx,dHdom,y) 
  CALL prodmv(donx,dHdon,y)
! ------------ preparing dqags call ------------
  epsabs=1.d-8 
  epsrel=1.d-5 
  limit=limx 
  lenw=limx4
!  ------------ selecting derivative  ------------
  IF(ider.eq.1)THEN 
     call dqagsc(fg,0.d0,dpig,epsabs,epsrel,result,abserr,neval,    &
          &        ier,limit,lenw,last,iwork,work)                           
! if(ier.gt.1)write(*,*)'d/G ',ier,u 
  ELSEIF(ider.eq.2)THEN 
     call dqagsc(fom,0.d0,dpig,epsabs,epsrel,result,abserr,neval,    &
          &        ier,limit,lenw,last,iwork,work)                           
! if(ier.gt.1)write(*,*)'d/dom ',ier,u                              
  ELSEIF(ider.eq.3)THEN 
     call dqagsc(fz,0.d0,dpig,epsabs,epsrel,result,abserr,neval,    &
          &        ier,limit,lenw,last,iwork,work)                           
! if(ier.gt.1)write(*,*)'d/dZ ',ier,u
  ELSEIF(ider.eq.4)THEN 
     call dqagsc(fon,0.d0,dpig,epsabs,epsrel,result,abserr,neval,    &
          &        ier,limit,lenw,last,iwork,work)                           
! if(ier.gt.1)write(*,*)'d/donod ',ier,u
  ELSE 
     WRITE(*,*)' ffd: ider=',ider 
     STOP 
  ENDIF
 
  ffd=result

  RETURN 
END FUNCTION ffd
                
! ********************************************************************* 
! subroutine computing G derivative (dR/dG)
DOUBLE PRECISION FUNCTION fg(up)
! --------------- interface ---------------------
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------ end interface -----------------------
  REAL(KIND=dkind) :: cup,sup,xp(3),yp(3)
  REAL(KIND=dkind) :: d2,d3,drde,drdi
  REAL(KIND=dkind) :: inc,f0,f0_inc,g,g_inc,maux(3,3)
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
  fg = (gmp/d3)*(1.d0-ep*cup)*((-beta/(ky*dsqrt(aa)*ea))*drde + &
       &   (ci/(si*ky*dsqrt(aa)*beta))*drdi)

! **************************************************
!  CALL check_deriv(cup,xp,d2,fg)
! *************************************************
  RETURN 
END FUNCTION fg

! ***********************************************************************
! subroutine computing omega derivative (dR/dg)
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

! =========== DERIVATIVE with respect to omega ============== 
  fom = (gmp/d3)*(1.d0-ea*cu)*(1.d0-ep*cup)*drdom

! **************************************************
!  CALL check_deriv(cup,xp,d2,fom)
! **************************************************
  RETURN 
END FUNCTION fom
                                                                        
! ********************************************************************      
! subroutine computing Z derivative (dR/dZ)
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
  fz = (gmp/d3)*(-1.d0/(ky*beta*sqrt(aa)*si))*(1.d0-ea*cu)*(1.d0-ep*cup)*drdi

! **************************************************
!  CALL check_deriv(cup,xp,d2,fz)
! **************************************************
  RETURN 
END FUNCTION fz

! *********************************************************************
! subroutine computing Omnod.nod derivative (dR/donod)
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
SUBROUTINE secpert(n,elpl,om,omnod,g,zl,a,ddd,eee,nnn) 
! --------------- interface ---------------------
  INTEGER,INTENT(IN) :: n !planets number
  TYPE(orbit_elem),INTENT(IN) :: elpl !planet elems 
  REAL(KIND=dkind),INTENT(IN) :: om,omnod,g,zl,a
  REAL(KIND=dkind),INTENT(OUT) :: ddd,eee !ddd=hamiltonian value
  INTEGER,INTENT(OUT) :: nnn
! -------------- end interface -------------------
! for dqags                                                         
  INTEGER :: limx,limx4 
  PARAMETER (limx=500,limx4=4*limx)
  INTEGER :: ier,limit,iwork(limx),lenw,last 
! function evaluations                                              
  INTEGER :: neval 
  REAL(KIND=dkind) :: epsabs,epsrel,abserr,work(limx4) 
! output of dqags                                                 
  REAL(KIND=dkind) :: rm
! ======================================================
  gmp = gm(n)
  CALL rot_matrix(elpl,om,omnod,g,zl,a) 
! initialization
  ddd=0.d0 
  eee=0.d0 
  nnn=0 
! preparing dqags call
  epsabs=1.d-9 
  epsrel=1.d-7 
  limit=limx 
  lenw=limx4
!     ========= perturbing function ==========
! per il momento uso l'integratore dqags anziche' dqagp
!  CALL dqags(ff,0.d0,dpig,npts2,points,epsabs,epsrel,rm,abserr,  &
!       &        neval,ier,leniw,lenw,last,iwork,work)                     
!  ddd=ddd+rm/(dpig**2) 
!  eee=eee+abserr/(dpig**2) 
!  nnn=nnn+neval
  CALL dqags(ffd_pert,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
       &        limit,lenw,last,iwork,work)                               
  ddd=ddd+rm/(dpig**2) 
  eee=eee+abserr/(dpig**2) 
  nnn=nnn+neval
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
  REAL(KIND=dkind) :: epsabs,epsrel,abserr,work(limx4),result 
! ----------------------------------------------
  cu=cos(u)
  su=sin(u)
! Asteroid coord
  y(1) = aa*(cu-ea)
  y(2) = aa*su*beta
  y(3) = 0.d0
  CALL prodmv(x,Hast,y)
! ------------ preparing dqagsc call ------------
  epsabs=1.d-8 
  epsrel=1.d-5 
  limit=limx 
  lenw=limx4
! -------------------------------------------------------
  CALL dqagsc(fun_pert,0.d0,dpig,epsabs,epsrel,result,abserr,neval,    &
          &        ier,limit,lenw,last,iwork,work) 
!     if(ier.gt.1)write(*,*)'pert fun ',ier,u

  ffd_pert=result  

  RETURN 
END FUNCTION ffd_pert

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
  fun_pert = (gmp/sqrt(d2))*(1.d0-ea*cu)*(1.d0-ep*cup)

  RETURN 
END FUNCTION fun_pert

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
  REAL(KIND=dkind) :: epsabs,epsrel,abserr,work(limx4) 
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
  epsabs=1.d-8 
  epsrel=1.d-5 
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
  REAL(KIND=dkind) :: epsabs,epsrel,abserr,work(limx4),result 
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
  epsabs=1.d-8 
  epsrel=1.d-5 
  limit=limx 
  lenw=limx4
!  ------------ selecting derivative  ------------
  IF(ider.eq.5)THEN 
     call dqagsc(fg_ii,0.d0,dpig,epsabs,epsrel,result,abserr,neval,    &
          &        ier,limit,lenw,last,iwork,work) 
!     if(ier.gt.1)write(*,*)'d/dom ',ier,u
  ELSEIF(ider.eq.6)THEN
     call dqagsc(fom_ii,0.d0,dpig,epsabs,epsrel,result,abserr,neval,    &
          &        ier,limit,lenw,last,iwork,work)
! if(ier.gt.1)write(*,*)'d/dom ',ier,u 
  ELSEIF(ider.eq.7)THEN
     call dqagsc(fz_ii,0.d0,dpig,epsabs,epsrel,result,abserr,neval,    &
          &        ier,limit,lenw,last,iwork,work)
! if(ier.gt.1)write(*,*)'d/dZ ',ier,u   
  ELSEIF(ider.eq.8)THEN
     call dqagsc(fon_ii,0.d0,dpig,epsabs,epsrel,result,abserr,neval,    &
          &        ier,limit,lenw,last,iwork,work)
! if(ier.gt.1)write(*,*)'d/donod ',ier,u
  ELSE 
     WRITE(*,*)' ffd_ii: ider=',ider 
     STOP 
  ENDIF
 
  ffd_ii=result
  
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


SUBROUTINE pladat 
!  IMPLICIT NONE 
! planet data                                                           
!      INCLUDE 'pldata.h90' 
  DOUBLE PRECISION gk,gk2,gjyr,gjyr2 
  INTEGER i 
! ======== PIANETI presi in considerazione ==============               
  inpl=2 
  ioupl=6 
  npl=ioupl-inpl+1 
! ============ planet semimajor axis ==================                 
  ap(1)= .3870992058d0 
  ap(2)= .7233274811d0 
  ap(3)= 1.0000036214d0 
  ap(4)= 1.5235973464d0 
  ap(5)= 5.2024107723d0 
  ap(6)= 9.5575876779d0 
  ap(7)= 19.3008879212d0 
  ap(8)= 30.2722024706d0 
  ap(9)= 39.7533710065d0 
! ============ masse dei pianeti =======================                
  gm(1)= 0.00000016601d0 
  gm(2)= 0.00000244781d0 
  gm(3)= 0.0000030404d0 
  gm(4)= 0.00000032272d0 
  gm(5)= 0.00095479d0 
  gm(6)= 0.00028589d0 
  gm(7)= 0.000043662d0 
  gm(8)= 0.000051514d0 
  gm(9)= 0.0000000073964d0 
! =========== conversioni di unita' =====================               
!  Gauss constant                                                       
  gk=0.01720209895d0 
  gk2=gk*gk 
!  conversion to internal units: 1AU, 1JYR=365.25 d(Julian year)        
  gjyr=365.25d0 
  gjyr2=gjyr*gjyr 
  bigg=gk2*gjyr2 
  ky=gk*gjyr 
  DO i =1,9 
     gm(i)=gm(i)*bigg 
  ENDDO
!  RETURN 
END SUBROUTINE pladat


END MODULE right_hand_side
