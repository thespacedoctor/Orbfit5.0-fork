! =========================================================         
! RIGHT HAND SIDE COMPUTATION
! INTEGRATION DOMAIN: T = [-PIG,PIG] x [0,DPIG]
! REFERENCE SYSTEM: x-axis towards ascending node                   
SUBROUTINE rhsloop(elem,ddd,eee,nnn) 
  USE fund_const
  USE orbit_elements
  USE planet_orbits 
 IMPLICIT NONE 
! ===================== INTERFACE =================================
! input: asteroid orbital elements
  TYPE(orbit_elem),INTENT(IN) :: elem ! asteroid orbital element
! output: right hand side,errors,number eval                        
  REAL(KIND=dkind),INTENT(OUT) :: ddd(4),eee(4) 
  INTEGER,INTENT(OUT) :: nnn(4) 
! ======================= END INTERFACE ============================
! Asteroid orbital elements
  REAL(KIND=dkind) :: om,omnod
! Planet orbital elements
  TYPE(orbit_elem) :: pla
  REAL(KIND=dkind) :: ompl,ipl
! for dqags                                                         
  INTEGER :: limx,limx4 
  PARAMETER (limx=500,limx4=4*limx)
  INTEGER :: ier,limit,iwork(limx),lenw,last 
! function evaluations                                              
  INTEGER :: neval
  REAL(KIND=dkind) :: epsabs,epsrel,abserr,work(limx4) 
! output of dqags                                                   
  REAL(KIND=dkind) :: rm
! for placar
  INTEGER :: fail_flag
  REAL(KIND=dkind),DIMENSION(3,3) :: maux
! =======================================================================
! loop index                                                      
  INTEGER :: n,k
! =======================================================================
! external function                                                 
  EXTERNAL loopffd 
! =======================================================================
! Common
! functions of the asteroid elements
  REAL(KIND=dkind) :: aa,ea,beta,ci,si
  REAL(KIND=dkind),DIMENSION(3,3) :: R_om,R_i,R_on,dR_om,dR_i,dR_on,&
       & H,dHdI,dHdom,dHdon
  COMMON/elems/aa,ea,beta,ci,si,R_om,R_i,R_on,dR_om,dR_i,dR_on,&
       & H,dHdI,dHdom,dHdon
! functions of the planet elements                                 
  REAL(KIND=dkind) :: apl,epl,betap,gmp
  REAL(KIND=dkind),DIMENSION(3,3) :: R_omp,R_ip,Hp
  COMMON/elepla/apl,epl,betap,gmp,R_omp,R_ip,Hp
! for selecting derivatives
  INTEGER :: ider
  COMMON/derfla/ider
! =======================================================================
! planet data                                                       
  INCLUDE 'pldata.h90'
! ***********************************************************************
! Asteroid elements
  aa = elem%coord(1)
  ea = elem%coord(2)
  beta = sqrt(1.d0-ea**2) 
  om = elem%coord(5)
! Rotation matrix of angle omega for the asteroid
  R_om = 0.d0
  R_om(1,1) = cos(om)
  R_om(1,2) =-sin(om)
  R_om(2,1) = sin(om)
  R_om(2,2) = cos(om)
  R_om(3,3) = 1.d0
! Derivative of the rotation matrix of angle omega for the asteroid
  dR_om = 0.d0
  dR_om(1,1) =-sin(om)
  dR_om(1,2) =-cos(om)
  dR_om(2,1) = cos(om)
  dR_om(2,2) =-sin(om)
! Rotation matrix of angle Inclination for the asteroid
  ci = cos(elem%coord(3))
  si = sin(elem%coord(3))
  R_i = 0.d0
  R_i(1,1) = 1.d0
  R_i(2,2) = ci
  R_i(2,3) =-si
  R_i(3,2) =si
  R_i(3,3) =ci
! Derivative of the rotation matrix of angle Inclination for the asteroid
  dR_i = 0.d0
  dR_i(2,2) =-si
  dR_i(2,3) =-ci
  dR_i(3,2) = ci
  dR_i(3,3) =-si
! ================= initialization ====================================
  DO k=1,4
     ddd(k) = 0.d0 
     eee(k) = 0.d0 
     nnn(k) = 0 
  ENDDO
! ===================== loop on number of planets =====================
  DO 10 n=inpl,ioupl
     pla = el_pla(n)
! planet elements
     apl = pla%coord(1)
     epl = pla%coord(2)
     ipl = pla%coord(3)
     ompl = pla%coord(5) 
     betap = sqrt(1.d0-epl**2)
     gmp = gm(n) 
! Rotation matrix of angle omega for the planet
     R_omp = 0.d0
     R_omp(1,1) = cos(ompl)
     R_omp(1,2) =-sin(ompl)
     R_omp(2,1) = sin(ompl)
     R_omp(2,2) = cos(ompl)
     R_omp(3,3) = 1.d0
! Rotation matrix of angle Inclination for the planet
     R_ip = 0.d0
     R_ip(1,1) = 1.d0
     R_ip(2,2) = cos(ipl)
     R_ip(2,3) =-sin(ipl)
     R_ip(3,2) = sin(ipl)
     R_ip(3,3) = cos(ipl)
 
    CALL prodmm(Hp,R_ip,R_omp)
     
! Rotation matrix of angle Omega.nodal for the asteroid
     omnod = elem%coord(4)-pla%coord(4)
     R_on = 0.d0
     R_on(1,1) = cos(omnod)
     R_on(1,2) =-sin(omnod)
     R_on(2,1) = sin(omnod)
     R_on(2,2) = cos(omnod)
     R_on(3,3) = 1.d0
! Derivative of the rotation matrix of angle Omega.nodal for the asteroid
     dR_on = 0.d0
     dR_on(1,1) =-sin(omnod)
     dR_on(1,2) =-cos(omnod)
     dR_on(2,1) = cos(omnod)
     dR_on(2,2) =-sin(omnod)
                  
  CALL prodmm(maux,R_on,R_i)
  CALL prodmm(H,maux,R_om) 

  CALL prodmm(maux,R_on,dR_i)
  CALL prodmm(dHdI,maux,R_om)

  CALL prodmm(maux,R_on,R_i)
  CALL prodmm(dHdom,maux,dR_om)

  CALL prodmm(maux,dR_on,R_i)
  CALL prodmm(dHdon,maux,R_om)
                            
! =================== preparing dqags call ===========================
     epsabs=1.d-8 
     epsrel=1.d-5 
     limit=limx 
     lenw=limx4 

! ********************************************************************     
!                    COMPUTING DERIVATIVES 
! ********************************************************************     
! first computing integral with respect to u
                                                                        
! Derivative with respect to G                                      
     ider=1 
     CALL dqags(loopffd,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,&
          &        limit,lenw,last,iwork,work)                               
     ddd(1)=ddd(1)+rm/(dpig**2) 
     eee(1)=eee(1)+abserr/(dpig**2) 
     nnn(1)=nnn(1)+neval 
                                                                        
! Derivative with respect to omega                                  
     ider=2 
     CALL dqags(loopffd,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,&
          &        limit,lenw,last,iwork,work)                               
     ddd(2)=ddd(2)+rm/(dpig**2) 
     eee(2)=eee(2)+abserr/(dpig**2) 
     nnn(2)=nnn(2)+neval 
                                                                        
! Derivative with respect to Z                                      
     ider=3 
     CALL dqags(loopffd,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,&
          &        limit,lenw,last,iwork,work)                               
     ddd(3)=ddd(3)+rm/(dpig**2) 
     eee(3)=eee(3)+abserr/(dpig**2) 
     nnn(3)=nnn(3)+neval 

! derivative with respect to Omega_nod
     ider=4 
     CALL dqags(loopffd,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
          &        limit,lenw,last,iwork,work)                               
     ddd(4)=ddd(4)+rm/(dpig**2) 
     eee(4)=eee(4)+abserr/(dpig**2) 
     nnn(4)=nnn(4)+neval 

10 END DO
                                                                        
  RETURN 
END SUBROUTINE rhsloop
                                                                        
! ********************************************************************     
! subroutine computing unidim.average of the unidim.average         
DOUBLE PRECISION FUNCTION loopffd(u)  
  USE fund_const
  IMPLICIT NONE 
  REAL(KIND=dkind) :: u ! eccentric anomaly of the asteroid
  REAL(KIND=dkind),DIMENSION(3,3) :: maux
  REAL(KIND=dkind) :: dy(3)
! for dqags                                                         
  INTEGER :: limx,limx4 
  PARAMETER (limx=500,limx4=4*limx) 
  INTEGER :: neval,ier,limit,iwork(limx),lenw,last
  REAL(KIND=dkind) :: epsabs,epsrel,abserr,work(limx4),risult 
! =====================================================================   
! external function                                                 
  EXTERNAL loopfg,loopfom,loopfz,loopfon 
! =====================================================================
! common for ff
! functions of the asteroid elements
  REAL(KIND=dkind) :: aa,ea,beta,ci,si
  REAL(KIND=dkind),DIMENSION(3,3) :: R_om,R_i,R_on,dR_om,dR_i,dR_on, &
       & H,dHdI,dHdom,dHdon
  COMMON/elems/aa,ea,beta,ci,si,R_om,R_i,R_on,dR_om,dR_i,dR_on,& 
       H,dHdI,dHdom,dHdon
! functions of the planet elements
  REAL(KIND=dkind) :: apl,epl,betap,gmp 
  REAL(KIND=dkind),DIMENSION(3,3) :: R_omp,R_ip,Hp
  COMMON/elepla/apl,epl,betap,gmp,R_omp,R_ip,Hp
! common per passare ad f
  REAL(KIND=dkind) :: x(3),y(3),cu,su,dx(3)                   
  COMMON/elele/x,y,cu,su,dx
! for selecting derivatives
  INTEGER :: ider                          
  COMMON/derfla/ider 
! ===================== u depending computations ======================    
  cu = cos(u)
  su = sin(u)

! Asteroid coordinate
  y(1) = aa*(cu-ea)
  y(2) = aa*su*beta
  y(3) = 0.d0
  CALL prodmv(x,H,y)

! Asteroid derivative coord w.r.t. ecc
  dy(1) =-aa
  dy(2) =-aa*su*ea/beta
  dy(3) = 0.d0
  CALL prodmv(dx,H,dy)

! ================= preparing dqags call ===============================
  epsabs=1.d-8 
  epsrel=1.d-5
  limit=limx 
  lenw=limx4 
! ====================== selecting derivative =========================   
  IF(ider.eq.1)THEN 
     call dqagsc(loopfg,-pig,pig,epsabs,epsrel,risult,abserr,neval, &
          &        ier,limit,lenw,last,iwork,work)                           
!     if(ier.gt.1)write(*,*)'d/G ',ier,u                                
  ELSEIF(ider.eq.2)THEN 
     call dqagsc(loopfom,-pig,pig,epsabs,epsrel,risult,abserr,neval, &
          &        ier,limit,lenw,last,iwork,work)                           
!     if(ier.gt.1)write(*,*)'d/dom ',ier,u                              
  ELSEIF(ider.eq.3)THEN 
     call dqagsc(loopfz,-pig,pig,epsabs,epsrel,risult,abserr,neval, &
          &        ier,limit,lenw,last,iwork,work)                           
!     if(ier.gt.1)write(*,*)'d/dZ ',ier,u  
  ELSEIF(ider.eq.4)THEN 
     call dqagsc(loopfon,0.d0,dpig,epsabs,epsrel,risult,abserr,neval,    &
          &        ier,limit,lenw,last,iwork,work)                           
!     if(ier.gt.1)write(*,*)'d/domnod ',ier,u                           
  ELSE 
     WRITE(*,*)' loopffd: ider=',ider 
     STOP 
  ENDIF
                                                                        
  loopffd=risult 
  
  RETURN 
END FUNCTION loopffd
                                                                        
! ********************************************************************* 
! subroutine computing G derivative (dR/dG)
DOUBLE PRECISION FUNCTION loopfg(up) 
  USE fund_const
  IMPLICIT NONE 
  REAL(KIND=dkind) :: up  ! eccentric anomaly of the planet
  REAL(KIND=dkind),DIMENSION(3,3) :: H1,H2,maux,mmaux
  REAL(KIND=dkind) :: xp(3),yp(3)  ! Planet coordinate vectors
  REAL(KIND=dkind) :: d2,d3,drde,drdi,vaux(3),cup,sup,dy(3),s
! =====================================================================
! functions of the asteroid elements
  REAL(KIND=dkind) :: aa,ea,beta,ci,si
  REAL(KIND=dkind),DIMENSION(3,3) :: R_om,R_i,R_on,dR_om,dR_i,dR_on, &
       & H,dHdI,dHdom,dHdon
  COMMON/elems/aa,ea,beta,ci,si,R_om,R_i,R_on,dR_om,dR_i,dR_on, &
       & H,dHdI,dHdom,dHdon
! functions of the planet elements
  REAL(KIND=dkind) :: apl,epl,betap,gmp
  REAL(KIND=dkind),DIMENSION(3,3) :: R_omp,R_ip,Hp
  COMMON/elepla/apl,epl,betap,gmp,R_omp,R_ip,Hp
! common per passare ad f
  REAL(KIND=dkind) :: x(3),y(3),cu,su,dx(3)                  
  COMMON/elele/x,y,cu,su,dx
! ======================================================================
  INCLUDE 'pldata.h90' 
! ================= up depending computations ==========================
  cup = cos(up)
  sup = sin(up)

! Planet coordinate in the orbital plane
  yp(1) = apl*(cup-epl)
  yp(2) = apl*sup*betap
  yp(3) = 0.d0

! Planet coordinate in the rotated system
  CALL prodmv(xp,Hp,yp)  

  d2=(xp(1)-x(1))**2+(xp(2)-x(2))**2+(xp(3)-x(3))**2
  d3=d2**(1.5d0)       

! =================== DERIVATA with respect to e ========================
! ((-1/2)*der(D^2)*(1-e2*cosu)+D^2*der(1-e2*cosu))*(1-e1*cos(up))
! s = (-1/2)*der(D^2) = <Xpl-Xast, R_on*R_i*R_om*dx>

  s=DOT_PRODUCT((xp-x),dx)
  drde=s*(1.d0-ea*cu)-cu*d2

! ===================DERIVATA with respect to I =========================
! der(D^2)*(1-e2*cosu)*(1-e1*cos(up))
! -2*s = der(D^2) = -2*<Xpl-Xast, R_on*dR_i*R_om*Y>

  CALL prodmv(vaux,dHdI,y)
  s=DOT_PRODUCT((xp-x),vaux)
  drdi=s*(1.d0-ea*cu)
      
! ==================DERIVATA with respect to G ==========================

  loopfg = (gmp/d3)*(1.d0-epl*cup)*((-beta/(ky*sqrt(aa)*ea))*drde+&
       &   (ci/(si*ky*sqrt(aa)*beta))*drdi) 
  
  RETURN 
END FUNCTION loopfg

! ***********************************************************************
! subroutine computing omega derivative (dR/dg)
DOUBLE PRECISION FUNCTION loopfom(up)
  USE fund_const
  IMPLICIT NONE                                   
  REAL(KIND=dkind) :: up  ! eccentric anomaly of the planet
  REAL(KIND=dkind) :: xp(3),yp(3)  ! Planet coordinate vectors
  REAL(KIND=dkind) :: d3,drdom,vaux(3),cup,sup,s
! ========================================================================
! common asteroid elements
  REAL(KIND=dkind) :: aa,ea,beta,ci,si
  REAL(KIND=dkind),DIMENSION(3,3) :: R_om,R_i,R_on,dR_om,dR_i,dR_on, &
       & H,dHdI,dHdom,dHdon
  COMMON/elems/aa,ea,beta,ci,si,R_om,R_i,R_on,dR_om,dR_i,dR_on, &
       & H,dHdI,dHdom,dHdon 
! common planet elements
  REAL(KIND=dkind) :: apl,epl,betap,gmp
  REAL(KIND=dkind),DIMENSION(3,3) :: R_omp,R_ip,Hp
  COMMON/elepla/apl,epl,betap,gmp,R_omp,R_ip,Hp
! common per passare ad f
  REAL(KIND=dkind) :: x(3),y(3),cu,su,dx(3)
  COMMON/elele/x,y,cu,su,dx
! ================== up depending computations ========================
  cup = cos(up)
  sup = sin(up)

! Planet coordinate in the orbital plane
  yp(1) = apl*(cup-epl)
  yp(2) = apl*sup*betap
  yp(3) = 0.d0

! Planet coordinate in the rotate system
  CALL prodmv(xp,Hp,yp)
  d3=((xp(1)-x(1))**2+(xp(2)-x(2))**2+(xp(3)-x(3))**2)**(1.5d0)
                                                                    
! =============== DERIVATIVE with respect to omega ====================
! s = (-1/2)*der(D^2) = -2*<Xpl-Xast, R_on*R_i*dR_om*Y>

  CALL prodmv(vaux,dHdom,y)
  drdom=DOT_PRODUCT((xp-x),vaux)

  loopfom = (gmp/d3)*(1.d0-ea*cu)*(1.d0-epl*cup)*drdom

  RETURN 
END FUNCTION loopfom
                                                                        
! ********************************************************************      
! subroutine computing Z derivative (dR/dZ)
DOUBLE PRECISION FUNCTION loopfz(up) 
  USE fund_const
  IMPLICIT NONE                                   
  REAL(KIND=dkind) :: up  ! eccentric anomaly of the planet
  REAL(KIND=dkind) :: xp(3),yp(3)  ! Planet coordinate vectors
  REAL(KIND=dkind) :: d3,drdi,vaux(3),s,cup,sup
! ===================================================================    
! common asteroid elements
  REAL(KIND=dkind) :: aa,ea,beta,ci,si
  REAL(KIND=dkind),DIMENSION(3,3) :: R_om,R_i,R_on,dR_om,dR_i,dR_on, &
       & H,dHdI,dHdom,dHdon
  COMMON/elems/aa,ea,beta,ci,si,R_om,R_i,R_on,dR_om,dR_i,dR_on,& 
       & H,dHdI,dHdom,dHdon
! common planet elements
  REAL(KIND=dkind) :: apl,epl,betap,gmp
  REAL(KIND=dkind),DIMENSION(3,3) :: R_omp,R_ip,Hp
  COMMON/elepla/apl,epl,betap,gmp,R_omp,R_ip,Hp
! common per passare ad f
  REAL(KIND=dkind) :: x(3),y(3),cu,su,dx(3)                     
  COMMON/elele/x,y,cu,su,dx
! =====================================================================
  INCLUDE 'pldata.h90'
! ================== up depending computations ========================
  cup = cos(up)
  sup = sin(up)

! Planet coordinate in the orbital plane
  yp(1) = apl*(cup-epl)
  yp(2) = apl*sup*betap
  yp(3) = 0.d0
! Planet coordinate in the rotated system
  CALL prodmv(xp,Hp,yp)
 
  d3=((xp(1)-x(1))**2+(xp(2)-x(2))**2+(xp(3)-x(3))**2)**(1.5d0)
                                                              
! ================= DERIVATIVE with respect to I ======================
! s = (-1/2)*der(D^2) = <Xpl-Xast, R_on*dR_i*R_om*Y>

  CALL prodmv(vaux,dHdI,y)
  drdi=DOT_PRODUCT((xp-x),vaux)

! ================= DERIVATIVE with respect to Z ======================

  loopfz=(gmp/d3)*(-1.d0/(ky*beta*sqrt(aa)*si))*(1.d0-ea*cu)*(1.d0-epl*cup)*drdi

  RETURN 
END FUNCTION loopfz

! *********************************************************************
! subroutine computing Omnod.nod derivative (dR/donod)
! This derivative is zero when the orbit of the planet is circular
DOUBLE PRECISION FUNCTION loopfon(up) 
  USE fund_const
  IMPLICIT NONE                                   
  REAL(KIND=dkind) :: up  ! eccentric anomaly of the planet
  REAL(KIND=dkind) :: xp(3),yp(3)  ! Vectors planet coordinate
  REAL(KIND=dkind) :: d3,drdon,vaux(3),s,cup,sup
! =====================================================================
! common asteroid elements
  REAL(KIND=dkind) :: aa,ea,beta,ci,si
  REAL(KIND=dkind),DIMENSION(3,3) :: R_om,R_i,R_on,dR_om,dR_i,dR_on, &
       & H,dHdI,dHdom,dHdon
  COMMON/elems/aa,ea,beta,ci,si,R_om,R_i,R_on,dR_om,dR_i,dR_on,&
       & H,dHdI,dHdom,dHdon
! common planet elements
  REAL(KIND=dkind) :: apl,epl,betap,gmp
  REAL(KIND=dkind),DIMENSION(3,3) :: R_omp,R_ip,Hp
  COMMON/elepla/apl,epl,betap,gmp,R_omp,R_ip,Hp
! common per passare ad f
  REAL(KIND=dkind) :: x(3),y(3),cu,su,dx(3)                     
  COMMON/elele/x,y,cu,su,dx
! ================== up depending computations ========================
  cup = cos(up)
  sup = sin(up)

! Planet coordinate in the orbital plane
  yp(1) = apl*(cup-epl)
  yp(2) = apl*sup*betap
  yp(3) = 0.d0
! Planet coordinate in the rotated system
  CALL prodmv(xp,Hp,yp)

  d3=((xp(1)-x(1))**2+(xp(2)-x(2))**2+(xp(3)-x(3))**2)**(1.5d0)

! =============== DERIVATIVE with respect to Omega.nod ================= 
! s = (-1/2)*der(D^2) = <Xpl-Xast, dR_on*R_i*R_om*Y>

  CALL prodmv(vaux,dHdon,y)
  drdon=DOT_PRODUCT((xp-x),vaux)
  
  loopfon = (gmp/d3)*(1.d0-ea*cu)*(1.d0-epl*cup)*drdon
 
   RETURN 
END FUNCTION loopfon
