MODULE kantorovich
  USE fund_const
  USE orbit_elements
  USE critical_points
  IMPLICIT NONE
  PRIVATE
! orbital elems and their functions
  TYPE(orbit_elem) :: elem,elpl !planet and asteroid elems
  REAL(KIND=dkind) :: aa,ea,beta,beta2,ci,si
  REAL(KIND=dkind) :: apl,epl,betap,betap2,gmp
  REAL(KIND=dkind),DIMENSION(3,3) :: R_om,R_i,R_on,dR_om,dR_i,dR_on, &
	& Hast,dHdI,dHdom,dHdon,R_omp,R_ip,Hp
! dmintil functions
  INTEGER :: nmin  !number of minimal points
  REAL(KIND=dkind),DIMENSION(nminx) :: dmintil
  REAL(KIND=dkind),DIMENSION(4,nminx) :: dDeldmint2 !derivatives w.r.t. G,Z,g,z
  REAL(KIND=dkind),DIMENSION(nminx) :: u1sgn,u2sgn
  REAL(KIND=dkind),DIMENSION(2,nminx) :: Vh
! functions of min loc points
  REAL(KIND=dkind) :: p1,p2 !angle translation parameters 
  REAL(KIND=dkind) ::pinf,psup,ainf,asup !integration intervals
  REAL(KIND=dkind),DIMENSION(3) :: x1,x2
  REAL(KIND=dkind),DIMENSION(3) :: y2,duy2,dux2,dux1
  REAL(KIND=dkind),DIMENSION(2,2,nminx) :: Ah
! deltarhs
  REAL(KIND=dkind),DIMENSION(2,2,4,nminx) :: dDelA !DEL-derivs of Ah matrix
  REAL(KIND=dkind),DIMENSION(2,4,nminx) :: dDelVh !dervivatives w.r.t. G,Z,g,z
  REAL(KIND=dkind),DIMENSION(4,4):: dKEP_dDEL
  REAL(KIND=dkind),DIMENSION(3) :: xa,ya,dexa,dIxa,domxa,donxa
  REAL(KIND=dkind),DIMENSION(2,nminx) :: kappa
  REAL(KIND=dkind) :: cu,su,dl
  INTEGER :: ider
! approxrhs
  REAL(KIND=dkind) :: Delta,rho,sigma,tau,lbar1,lbar2!,detAh
  REAL(KIND=dkind),DIMENSION(4) :: dDeldetA,dDelDelta,dDelrho,dDelsigma,dDeltau
  REAL(KIND=dkind),DIMENSION(2,4) :: dlbar
  REAL(KIND=dkind) :: dmint2 !dmintil and its square
  REAL(KIND=dkind),DIMENSION(4) :: ddmint2 !dmintil DEL-derivs 
  INTEGER :: ddel !DEL index indicating the variable of integration
  INCLUDE 'pldata.h90'

! variables
!  PUBLIC :: elpl,elem
! routines
  PUBLIC :: split,coord_tau_Ah ,rot_matrices

CONTAINS
! -----------------------------------------------------------------
! MOID DERIVATIVES, ITS MINIMUM AND ITS LOC MIN POINTS
! -----------------------------------------------------------------
SUBROUTINE split(n,elp,om,omnod,g,zl,aa,ddd,eee,nnn)
! ===================== INTERFACE ================================= 
  INTEGER, INTENT(IN) :: n  !planet number
  TYPE(orbit_elem),INTENT(IN) :: elp !planet elems
  REAL(KIND=dkind),INTENT(IN) :: om,omnod,g,zl,aa
! right hand side,errors,number eval                        
  REAL(KIND=dkind),INTENT(OUT) :: ddd(4),eee(4) 
  INTEGER,INTENT(OUT) :: nnn(4)
! ======================= END INTERFACE ===========================
  REAL(KIND=dkind),DIMENSION(4) :: dd,ee,dd1,ee1,dd_1,ee_1
  REAL(KIND=dkind) :: ddf,eef,ddf_rem
  INTEGER,DIMENSION(4) :: nn,nn1,nn_1
  INTEGER :: nnf
! Varibles in the argument of dmintil_rms
  INTEGER :: nummin  !number of minimal points 
  REAL(KIND=dkind),DIMENSION(5,nminx) :: ddmintdel2  !COM-derivatives
  REAL(KIND=dkind),DIMENSION(nminx) :: v1min,v2min,detH!,sint1t2
!  LOGICAL,PARAMETER :: chk_der=.true.
! angle tanslation parameters vectors
  REAL(KIND=dkind),DIMENSION(nminx) :: t1,t2
! =================================================================
!  LOGICAL :: check_der
!  REAL(KIND=dkind) :: incr,g1,g2,zeta
  INTEGER :: j
! =================================================================

!  check_der=.false.
!  check_der=.true.

! Planet KEP elements
  elpl = elp
  gmp = gm(n)

! Asteroid KEP elements at the intermediate time steps
  elem = undefined_orbit_elem
  elem%coo = 'KEP'
  elem%coord(1) = aa
  elem%coord(2) = sqrt(1.d0-(g/ky)**2/aa)
  elem%coord(3) = ACOS(zl/g)
  elem%coord(4) = omnod
  elem%coord(5) = om
  elem%coord(6) = 0.d0 
!  elem%t=elpl%t 

  CALL rot_matrices(elpl,elem)

  CALL dmintil_rms(elpl,elem,nummin,dmintil,DDMINTDEL2=ddmintdel2, &
       & V1MIN=v1min,V2MIN=v2min)

!  IF(nummin.gt.2)THEN
!     write(*,*)'number of minimal points: ',nummin
!     STOP
!  ELSE

  nmin = nummin
! ***** dummy value ****
!     nmin = 1
! *********************
  IF(nmin.gt.1)THEN
     IF(abs(dmintil(2)).lt.sqrt(gmp/39.4769))THEN
        nmin=2 ! no more than 2 singularity extraction
        WRITE(*,*)'DOUBLE CROSSING!!:',dmintil(1:2)
     ELSE
        nmin=1
     ENDIF
  ENDIF
  

  CALL delaunay_ddmintil(ddmintdel2) !output:dDeldmint2(1:4,1:nmin)

  CALL anomalies_at_min(v1min,v2min,u1sgn,u2sgn,t1,t2)

  CALL int_interval(t1,t2,nmin,pinf,psup,ainf,asup,p1,p2)

! ------------ initialization --------------
  dd1 = 0.d0
  ee1 = 0.d0
  nn1 = 0.d0
! ------------------------------------------

  DO j=1,nmin

! Mean anomalies vector at min loc points
     Vh(1,j) = u1sgn(j)-epl*sin(u1sgn(j))
     Vh(2,j) = u2sgn(j)- ea*sin(u2sgn(j))
! dmintil for approx apptoxrhs
     dmint2 = dmintil(j)**2
     ddmint2 = dDeldmint2(:,j)

     CALL coord_tau_Ah(u1sgn(j),u2sgn(j),Ah(:,:,j))

     CALL compute_derivs(u1sgn(j),u2sgn(j),dDelVh(:,:,j),dDelA(:,:,:,j))

! +++++++++++++++++++++ APPROXRHS2 +++++++++++++++++++++++++
     CALL approxrhs2(j,dd_1,ee_1,nn_1) !dd1=G,Z,om,omnod-derivs
     CALL riordina(dd_1) !change dd_1 in G,g,Z,z-derivs 
     CALL riordina(ee_1)
     CALL riordina_int(nn_1)
     dd1 = dd1 + dd_1
     ee1 = ee1 + ee_1
     nn1 = nn1 + nn_1

  ENDDO

! ********************************
! uncomment to check the derivatives of the average of 1/delta_h
!  CALL check_deriv_approx(dd_1)
!  write(*,*)'dDelA(2,2,1:4,1)',dDelA(2,2,1:4,1)
!  CALL checkder(dDelA(2,2,1:4,1),4)
! ********************************

! +++++++++++ DELTARHS2 +++++++++++++
  CALL deltarhs2(dd,ee,nn)
! dd = (G,g,Z,z)
  DO j = 1,4 
     ddd(j) = dd(j) + dd1(j)
     eee(j) = ee(j) + ee1(j)
     nnn(j) = nn(j) + nn1(j)
  ENDDO

END SUBROUTINE split

! =================================================================         
! APPROXIMATE PERTURBING FUNCTION (first loc min point)            
! INTEGRATION DOMAIN: T = [-pi,pi] x [-pi,pi]                        
SUBROUTINE app_pertfun(ddd,eee,nnn)
! ===================== INTERFACE =================================
  REAL(KIND=dkind),INTENT(OUT) :: ddd,eee !right hand side,errors
  INTEGER,INTENT(OUT) :: nnn !number eval
! ======================= END INTERFACE ===========================
! for dqagp
!  REAL(KIND=dkind) :: points
!  INTEGER :: npts2
! for dqags                                                         
  INTEGER,PARAMETER :: limx=300,limx4=4*limx
  INTEGER :: ier,leniw,iwork(limx),lenw,last 
! function evaluations                                              
  INTEGER :: neval
  REAL(KIND=dkind) :: epsabs,epsrel,abserr,work(limx4) 
! output of dqags                                                   
  REAL(KIND=dkind) :: rm
! =================== preparing dqags call ==========================
  epsabs = 1.d-8 
  epsrel = 1.d-5 
  leniw = limx 
  lenw = limx4
!  npts2 = 3
!  points = u2sgn
! ================= initialization ==================================
  ddd=0.d0 
  eee=0.d0 
  nnn=0 
! ********************************************************************     
! first computing integral with respect to l
!     ASTEROID INTEGRATION RANGE = [AINF,ASUP]                        
!  write(*,*)'integration extrema: ',ainf,asup
!  CALL dqagp(ffun,ainf,asup,npts2,points,epsabs,epsrel,rm, &
!	& abserr,neval,ier,leniw,lenw,last,iwork,work)

  CALL dqags(ffun,ainf,asup,epsabs,epsrel,rm, &
       & abserr,neval,ier,leniw,lenw,last,iwork,work)
  ddd = ddd + rm*(gmp/dpig**2) 
  eee = eee + abserr*(gmp/dpig**2) 
  nnn = nnn + neval 

END SUBROUTINE app_pertfun
 
! subroutine computing unidim.average of the unidim.average         
DOUBLE PRECISION FUNCTION ffun(u)
  REAL(KIND=dkind) :: u !asteroid eccentric anomaly
! ---------- end interface ----------
! for dqagp
!  REAL(KIND=dkind) :: points(3),eta
!  INTEGER :: npts2
! for dqags                                                         
  INTEGER, PARAMETER :: limx=500,limx4=4*limx
  INTEGER :: neval,ier,leniw,iwork(limx),lenw,last 
  REAL(KIND=dkind) :: epsabs,epsrel,abserr,work(limx4),risult
  INTEGER :: j !loop index 
! ===================== u depending computations ======================     
  cu = cos(u)
  su = sin(u)
  DO j=1,nmin
     kappa(2,j) = u-ea*su-Vh(2,j)
  ENDDO
! ================= preparing dqags call ===============================
  epsabs = 1.d-8 
  epsrel = 1.d-5
  leniw = limx
  lenw = limx4 
! deciding parameter                                                
!  eta = 1.d-7

  CALL dqagsc(fun,pinf,psup,epsabs,epsrel,risult,abserr, &
       & neval,ier,leniw,lenw,last,iwork,work)
  ffun = risult

  RETURN 
END FUNCTION ffun
  
! ********************************************************************* 
! Perturbative function
DOUBLE PRECISION FUNCTION fun(up)
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ----------------- end interface ------------------
  REAL(KIND=dkind) :: cup,sup
  REAL(KIND=dkind) :: xp(3),yp(3)  !planet coordinate vectors
  REAL(KIND=dkind) :: d2,d3,drde,drdi,delta2
  REAL(KIND=dkind) :: s,vaux(3),maux(2,2),vvaux(2),dGdelta,fg,fgd
  INTEGER :: j !loop index
! ==================================================
  cup = cos(up)
  sup = sin(up)

  fgd = 0.d0
  DO j=1,nmin
     kappa(1,j) = up-epl*sup - Vh(1,j)
! squared approximated distance
     CALL prodmv22(vvaux,Ah(:,:,j),kappa(:,j))
     delta2 = DOT_PRODUCT(kappa(:,j),vvaux) + dmintil(j)**2 
     fgd = fgd + 1.d0/sqrt(delta2)
  ENDDO
  fun = fgd*(1.d0-ea*cu)*(1.d0-epl*cup)
  
  RETURN 
END FUNCTION fun


! =================================================================         
! DERIVATIVE OF REMAINDER FUNCTION (first loc min point)            
! INTEGRATION DOMAIN: T = [-pi,pi] x [-pi,pi]                        
SUBROUTINE deltarhs2(ddd,eee,nnn)
! ------------------ inteface ---------------------
! right hand side,errors,number eval                        
  REAL(KIND=dkind),INTENT(OUT) :: ddd(4),eee(4) 
  INTEGER,INTENT(OUT) :: nnn(4)
! ----------------- end inteface ------------------
! for dqagp
  REAL(KIND=dkind) :: points(4)
  INTEGER :: npts2
! for dqags                                                         
  INTEGER,PARAMETER :: limx=300,limx4=4*limx
  INTEGER :: ier,leniw,iwork(limx),lenw,last 
! function evaluations                                              
  INTEGER :: neval
  REAL(KIND=dkind) :: epsabs,epsrel,abserr,work(limx4) 
! output of dqags                                                   
  REAL(KIND=dkind) :: rm
! loop indexes                                                      
  INTEGER :: j,k
! --- preparing dqags call ---
  epsabs = 1.d-8 
  epsrel = 1.d-5 
  leniw = limx 
  lenw = limx4
  npts2 = 3
  DO j = 1,nmin
     points(j) = u2sgn(j)
  ENDDO
! ------ initialization -------
  DO k=1,4
     ddd(k)=0.d0 
     eee(k)=0.d0 
     nnn(k)=0 
  ENDDO
! ********************************************************************     
!                    COMPUTING DERIVATIVES 
! ********************************************************************     
! first computing integral with respect to l
!     ASTEROID INTEGRATION RANGE = [AINF,ASUP]                        
                                                                       
! ------------------ Derivative w.r.t. G --------------------
  ider = 1 
  CALL dqagp(deltaff,ainf,asup,npts2,points,epsabs,epsrel,rm, &
	& abserr,neval,ier,leniw,lenw,last,iwork,work)
  ddd(1) = ddd(1) + rm*(gmp/dpig**2) 
  eee(1) = eee(1) + abserr*(gmp/dpig**2)
  nnn(1) = nnn(1) + neval 
! ------------------ Derivative w.r.t. g=omega --------------------
  ider = 2 
  CALL dqagp(deltaff,ainf,asup,npts2,points,epsabs,epsrel,rm, &
	& abserr,neval,ier,leniw,lenw,last,iwork,work)
  ddd(2) = ddd(2) + rm*(gmp/dpig**2) 
  eee(2) = eee(2) + abserr*(gmp/dpig**2) 
  nnn(2) = nnn(2) + neval 
! ----------------- Derivative w.r.t. Z ----------------
  ider = 3 
  CALL dqagp(deltaff,ainf,asup,npts2,points,epsabs,epsrel,rm, &
	& abserr,neval,ier,leniw,lenw,last,iwork,work)
  ddd(3) = ddd(3) + rm*(gmp/dpig**2) 
  eee(3) = eee(3) + abserr*(gmp/dpig**2) 
  nnn(3) = nnn(3) + neval 
! ---------------- Derivative w.r.t. z=Omnod -----------------
  ider = 4 
  CALL dqagp(deltaff,ainf,asup,npts2,points,epsabs,epsrel,rm, &
	& abserr,neval,ier,leniw,lenw,last,iwork,work)
  ddd(4) = ddd(4) + rm*(gmp/dpig**2) 
  eee(4) = eee(4) + abserr*(gmp/dpig**2) 
  nnn(4) = nnn(4) + neval 

END SUBROUTINE deltarhs2
                                                                        
! ********************************************************************     
! subroutine computing unidim.average of the unidim.average         
DOUBLE PRECISION FUNCTION deltaff(u) 
  REAL(KIND=dkind) :: u !asteroid eccentric anomaly
! ----------------- end interface -------------------
  REAL(KIND=dkind) :: dy(3)
! for dqapc
  REAL(KIND=dkind) :: points(3),eta
  INTEGER :: npts2
! for dqags                                                         
  INTEGER, PARAMETER :: limx=500,limx4=4*limx
  INTEGER :: neval,ier,leniw,iwork(limx),lenw,last 
  REAL(KIND=dkind) :: epsabs,epsrel,abserr,work(limx4),risult
  INTEGER :: j !loop index 
! ===============================================================     
  cu = cos(u)
  su = sin(u)
  DO j=1,nmin
     kappa(2,j)  = u-ea*su-Vh(2,j)
  ENDDO
!  dl = 1.d0-ea*cu !dl/du, l=u-ea*sin(u)
  
! Asteroid coordinate in the orbital plane
  ya(1) = aa*(cu-ea)
  ya(2) = aa*beta*su
  ya(3) = 0.d0

! Asteroid coordinate in the rotated system
  CALL prodmv(xa,Hast,ya)

! Asteroid coord derivative w.r.t. e in the orbital plane
  dy(1) =-aa
  dy(2) =-aa*ea/beta*su
  dy(3) = 0.d0

! Asteroid coord derivatives in the rotated system
  CALL prodmv(dexa,Hast,dy)   !e-deriv
  CALL prodmv(dIxa,dHdI,ya)   !I-deriv
  CALL prodmv(domxa,dHdom,ya) !omega-deriv
  CALL prodmv(donxa,dHdon,ya) !Omnod-deriv

! preparing dqags call
  epsabs = 1.d-8 
  epsrel = 1.d-5
  leniw = limx
  lenw = limx4 
  npts2 = 3
  DO j=1,nmin
     points(j) = u1sgn(j)
  ENDDO
! deciding parameter
  eta = 1.d-7 
! ====================== selecting derivative =========================
!     PLANET INTEGRATION RANGE = [PINF,PSUP]

  IF(ider.eq.1) THEN
     IF(ABS(u2sgn(1) - u).le.eta)THEN
        CALL dqagpc(deltafg,pinf,psup,npts2,points,epsabs,epsrel, &
             & risult,abserr,neval,ier,leniw,lenw,last,iwork,work)
! if(ier.gt.1) write(*,*)'d/G ',ier,f
     ELSE !IF(ABS(u2sgn(1) - u).gt.eta)THEN 
        CALL dqagsc(deltafg,pinf,psup,epsabs,epsrel,risult,abserr, &
             & neval,ier,leniw,lenw,last,iwork,work)
! if(ier.gt.1) write(*,*)'d/G ',ier,f
     ENDIF

  ELSEIF(ider.eq.2)THEN
     IF(ABS(u2sgn(1) - u).le.eta)THEN 
     	CALL dqagpc(deltafom,pinf,psup,npts2,points,epsabs,epsrel, &
           & risult,abserr,neval,ier,leniw,lenw,last,iwork,work)
! if(ier.gt.1) write(*,*)'d/dom ',ier,f
     ELSE !IF(ABS(u2sgn(1) - u).ge.eta)THEN
        CALL dqagsc(deltafom,pinf,psup,epsabs,epsrel,risult,abserr, &
             & neval,ier,leniw,lenw,last,iwork,work)
! if(ier.gt.1) write(*,*)'d/dom ',ier,f                                
     ENDIF

  ELSEIF(ider.eq.3) THEN
     IF(ABS(u2sgn(1) - u).le.eta)THEN 
        CALL dqagpc(deltafz,pinf,psup,npts2,points,epsabs,epsrel, &
             & risult,abserr,neval,ier,leniw,lenw,last,iwork,work)
! if(ier.gt.1) write(*,*)'d/dZ ',ier,f
     ELSE !IF(ABS(u2sgn(1) -u).ge.eta)THEN
        CALL dqagsc(deltafz,pinf,psup,epsabs,epsrel,risult,abserr, &
             & neval,ier,leniw,lenw,last,iwork,work)
! if(ier.gt.1) write(*,*)'d/dZ ',ier,f
     ENDIF

  ELSEIF(ider.eq.4)THEN
     IF(ABS(u2sgn(1) - u).le.eta)THEN 
        CALL dqagpc(deltafon,pinf,psup,npts2,points,epsabs,epsrel, &
             & risult,abserr,neval,ier,leniw,lenw,last,iwork,work)
! if(ier.gt.1) write(*,*)'d/donod ',ier,f
     ELSE !IF(ABS(u2sgn(1) - u).ge.eta)THEN                      
        CALL dqagsc(deltafon,pinf,psup,epsabs,epsrel, &
             & risult,abserr,neval,ier,leniw,lenw,last,iwork,work)
!  if(ier.gt.1) write(*,*)'d/donod ',ier,f
     ENDIF
     
  ELSE 
     WRITE(*,*)' deltaffd: ider=',ider 
     STOP 
  ENDIF
  
  deltaff = risult
  RETURN
END FUNCTION deltaff
                
! ********************************************************************* 
! subroutine computing G derivative (dR/dG)
DOUBLE PRECISION FUNCTION deltafg(up)
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------- end interface ---------------
  REAL(KIND=dkind) :: cup,sup!,dlp
  REAL(KIND=dkind) :: xp(3),yp(3)  !planet coord vectors
  REAL(KIND=dkind) :: d2,d3,drde,drdi,delta2,dGdelta,fg,fgd
  REAL(KIND=dkind) :: s,vaux(3),maux(2,2),vvaux(2)
  INTEGER :: j !lopp index
! =============================================================
  cup = cos(up)
  sup = sin(up)
  DO j=1,nmin
     kappa(1,j) = up-epl*sup - Vh(1,j)
  ENDDO
!  dlp = 1.d0-epl*cup !dlp/dup, lp=1-epl*sin(upl)

! Planet coordinate
  yp(1) = apl*(cup-epl)
  yp(2) = apl*betap*sup
  yp(3) = 0.d0
  CALL prodmv(xp,Hp,yp) !coord in the rotated system

! ======================== INTEGRAND FUNCTION =========================
! keplerian square and cube distance
  d2 = (xp(1)-xa(1))**2+(xp(2)-xa(2))**2+(xp(3)-xa(3))**2
  d3 = d2**(1.5d0)      
! derivative of D**2*dl/du w.r.t. e
  s = DOT_PRODUCT((xp-xa),dexa)
  drde = s*(1.d0-ea*cu) - d2*cu
! derivative of D**2*dl/du w.r.t. I
  s = DOT_PRODUCT((xp-xa),dIxa)
  drdi = s*(1.d0-ea*cu)
! derivative of 1/D*dl/du w.r.t. G
  fg = 1.d0/d3*(drde*dKEP_dDEL(1,1) + drdi*dKEP_dDEL(2,1))
! -------------------------------------------------------------
! squared approximated distance at each min loc point
  fgd = 0.d0
  DO j=1,nmin
     CALL prodmv22(vvaux,Ah(:,:,j),kappa(:,j))
     delta2 = DOT_PRODUCT(kappa(:,j),vvaux) + dmintil(j)**2  
! ---- derivative of d**2*dl/du w.r.t. G ----
     dGdelta = su*vvaux(2)*dKEP_dDEL(1,1) + DOT_PRODUCT(dDelVh(:,1,j),vvaux)
     maux = dDelA(:,:,1,j)
     CALL prodmv22(vvaux,maux,kappa(:,j))
     dGdelta = (dGdelta - 0.5d0*DOT_PRODUCT(kappa(:,j),vvaux) - &
          & dDeldmint2(1,j))*(1.d0-ea*cu)
! ---- derivative of 1/d*dl/du w.r.t. G ---- 
     fgd = fgd + 1.d0/delta2**(1.5d0)*(dGdelta-delta2*cu*dKEP_dDEL(1,1))
  ENDDO

  deltafg = (fg-fgd)*(1.d0-epl*cup)

  RETURN 
END FUNCTION deltafg
                                                                        
! ***********************************************************************
! subroutine computing omega derivative (dR/dg)
DOUBLE PRECISION FUNCTION deltafom(up)
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------- end interface ---------------
  REAL(KIND=dkind) :: cup,sup!,dlp
  REAL(KIND=dkind) :: xp(3),yp(3)  !planet coord vectors
  REAL(KIND=dkind) :: d3,drdom,delta3,domdelta,fom,fomd
  REAL(KIND=dkind) :: vaux(3),vvaux(2),maux(2,2)
  INTEGER :: j !loop index
! ============================================================
  cup = cos(up)
  sup = sin(up)
  DO j=1,nmin
     kappa(1,j) = up-epl*sup - Vh(1,j)
  ENDDO
!  dlp = 1.d0-epl*cup

! Planet coordinate
  yp(1) = apl*(cup-epl)
  yp(2) = apl*betap*sup
  yp(3) = 0.d0
  CALL prodmv(xp,Hp,yp)

! cube keplerian distance
  d3 = ((xp(1)-xa(1))**2+(xp(2)-xa(2))**2+(xp(3)-xa(3))**2)**(1.5d0)
! derivative of D**2 w.r.t. omega
  drdom = DOT_PRODUCT((xp-xa),domxa)
  fom = 1.d0/d3*drdom
! ------------------------------------------------------------------
  fomd = 0.d0
  DO j=1,nmin
! cube approximated distance
     CALL prodmv22(vvaux,Ah(:,:,j),kappa(:,j))
     delta3 = (DOT_PRODUCT(kappa(:,j),vvaux) + dmintil(j)**2)**(1.5d0) 
! --- derivative of d**2 w.r.t. omega ---
     domdelta = DOT_PRODUCT(dDelVh(:,3,j),vvaux)
     maux = dDelA(:,:,3,j)
     CALL prodmv22(vvaux,maux,kappa(:,j))
     domdelta = domdelta-0.5d0*DOT_PRODUCT(kappa(:,j),vvaux)-dDeldmint2(3,j)
! --- derivative od 1/d w.r.t omega --- 
     fomd = fomd + 1.d0/delta3*domdelta
  ENDDO

  deltafom = (fom-fomd)*(1.d0-ea*cu)*(1.d0-epl*cup)

  RETURN 
END FUNCTION deltafom
                                                                        
! ********************************************************************      
! subroutine computing Z derivative (dR/dZ)
DOUBLE PRECISION FUNCTION deltafz(up) 
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! -------------- end interface -----------------
  REAL(KIND=dkind) :: cup,sup!,dlp
  REAL(KIND=dkind) :: xp(3),yp(3)  !planet coord vectors
  REAL(KIND=dkind) :: d3,drdi,s,delta3,dZdelta,fz,fzd
  REAL(KIND=dkind) :: vaux(3),vvaux(2),maux(2,2)
  INTEGER :: j !loop index
! ==========================================================
  cup = cos(up)
  sup = sin(up)
  DO j=1,nmin
     kappa(1,j) = up-epl*sup - Vh(1,j)
  ENDDO

! Planet coordinate
  yp(1) = apl*(cup-epl)
  yp(2) = apl*betap*sup
  yp(3) = 0.d0
  CALL prodmv(xp,Hp,yp)

! cube keplerian distance
  d3 = ((xp(1)-xa(1))**2+(xp(2)-xa(2))**2+(xp(3)-xa(3))**2)**(1.5d0)
! derivative of D**2 w.r.t. I
  drdi = DOT_PRODUCT((xp-xa),dIxa)
! derivative of 1/D w.r.t. Z
  fz = (1.d0/d3)*drdi*dKEP_dDEL(2,2)
! -----------------------------------------------------------------
  fzd =0.d0
  DO j=1,nmin
! cube approximated distance
     CALL prodmv22(vvaux,Ah(:,:,j),kappa(:,j))
     delta3 = (DOT_PRODUCT(kappa(:,j),vvaux) + dmintil(j)**2)**(1.5d0)
! --- derivative of d**2 w.r.t. Z -----
     dZdelta = DOT_PRODUCT(dDelVh(:,2,j),vvaux)
     maux = dDelA(:,:,2,j)
     CALL prodmv22(vvaux,maux,kappa(:,j))
     dZdelta = dZdelta-0.5d0*DOT_PRODUCT(kappa(:,j),vvaux)-dDeldmint2(2,j)
! --- derivative of 1/d w.r.t. Z -----
     fzd = fzd + 1.d0/delta3*dZdelta
  ENDDO

  deltafz = (fz-fzd)*(1.d0-ea*cu)*(1.d0-epl*cup)

  RETURN 
END FUNCTION deltafz

! *********************************************************************
! subroutine computing Omnod.nod derivative (dR/donod)
DOUBLE PRECISION FUNCTION deltafon(up) 
  REAL(KIND=dkind) :: up  !planet eccentric anomaly
! ------------- end interface ------------------
  REAL(KIND=dkind) :: cup,sup,dlp
  REAL(KIND=dkind) :: xp(3),yp(3)  !planet coord vectors
  REAL(KIND=dkind) :: d3,drdon,delta3,dondelta,fon,fond
  REAL(KIND=dkind) :: vaux(3),maux(2,2),vvaux(2)
  INTEGER :: j !lopp index
! ==========================================================
  cup = cos(up)
  sup = sin(up)
  DO j=1,nmin
     kappa(1,j) = up-epl*sup - Vh(1,j)
  ENDDO

! Planet coordinate
  yp(1) = apl*(cup-epl)
  yp(2) = apl*betap*sup
  yp(3) = 0.d0
  CALL prodmv(xp,Hp,yp)

! cube keplerian distance
  d3 = ((xp(1)-xa(1))**2+(xp(2)-xa(2))**2+(xp(3)-xa(3))**2)**(1.5d0)
! derivative of 1/D w.r.t. Omega.nod
  drdon = DOT_PRODUCT((xp-xa),donxa)
  fon = (1.d0/d3)*drdon
! -----------------------------------------------------------------
  fond = 0.d0
  DO j=1,nmin
! cube approximated distance
     CALL prodmv22(vvaux,Ah(:,:,j),kappa(:,j))
     delta3 = (DOT_PRODUCT(kappa(:,j),vvaux) + dmintil(j)**2)**(1.5d0)
! --- derivative of d**2 w.r.t. Omega.nod ----
     dondelta = DOT_PRODUCT(dDelVh(:,4,j),vvaux)
     maux = dDelA(:,:,4,j)
     CALL prodmv22(vvaux,maux,kappa(:,j))
     dondelta = dondelta-0.5d0*DOT_PRODUCT(kappa(:,j),vvaux)-dDeldmint2(4,j)
     fond = fond + (1.d0/delta3)*dondelta
  ENDDO

  deltafon = (fon-fond)*(1.d0-ea*cu)*(1.d0-epl*cup)

   RETURN 
END FUNCTION deltafon

! ******************************************************************
! APPROXIMATE RIGHT-HAND SIDE (using polar coordinates)
! INTEGRATION DOMAIN: T = [PINF,PSUP] x [AINF,ASUP]
SUBROUTINE approxrhs2(jh,ddd,eee,nnn)
  INTEGER,INTENT(IN) :: jh !index of current min point                      
  REAL(KIND=dkind),INTENT(OUT) :: ddd(4),eee(4) !right hand side,errors
  INTEGER,INTENT(OUT) :: nnn(4) !number eval  
! --------------------- end interface ----------------------
  REAL(KIND=dkind) :: vsize,vtau1,dv,detAh!,lbar1,lbar2
! integration extrema                                               
  REAL(KIND=dkind),DIMENSION(6) :: ct,st,t 
! for dqags                                                         
  INTEGER,PARAMETER :: limx=300,limx4=4*limx
  INTEGER :: ier(5),limit,iwork(limx),lenw,last 
! function evaluations                                              
  INTEGER,DIMENSION(5) :: neval
  REAL(KIND=dkind) :: epsabs,epsrel,abserr(5),work(limx4) 
! output of dqags                                                   
  REAL(KIND=dkind),DIMENSION(5) :: rm
! loop indexes                                                      
  INTEGER :: j,k
! ==================================================================
! Determinant of matrix Ah, Delta
  detAh = Ah(1,1,jh)*Ah(2,2,jh)-Ah(1,2,jh)**2	
  IF(detAh.le.0.d0)THEN 
     write(*,*)' APPROXRHS2, detA=',detAh 
     STOP 
  ENDIF
  Delta = 1.d0/(sqrt(detAh))
  vtau1 = sqrt(Ah(1,1,jh))
! --------------------------------------------------------------------
! rho = rho4 ; sigma = rho3 ; tau = rho1
  rho = 1.d0/(Delta*vtau1)
  sigma = Ah(1,2,jh)/vtau1
  tau = vtau1
! --------------------------------------------------------------------
  DO k=1,4
! Derivative of detAh,Delta w.r.t. DEL
     dDeldetA(k) = dDelA(1,1,k,jh)*Ah(2,2,jh) + Ah(1,1,jh)*dDelA(2,2,k,jh)- &
          & 2.d0*Ah(1,2,jh)*dDelA(1,2,k,jh)
     dDelDelta(k) = -0.5d0*Delta**3*dDeldetA(k)
! Derivatives of rho,sigma,tau w.r.t. DEL
     dv = -0.5d0*dDelA(1,1,k,jh)/vtau1**3 !dv=d(1/vtau1))
     dDelrho(k) = 0.5d0*Delta/vtau1*dDeldetA(k) + dv/Delta
     dDelsigma(k) = dDelA(1,2,k,jh)/vtau1 + Ah(1,2,jh)*dv
     dDeltau(k) = 0.5d0/vtau1*dDelA(1,1,k,jh)
  ENDDO
! --------------------------------------------------------------------
  lbar1 = Vh(1,jh)
  lbar2 = Vh(2,jh)
  dlbar = dDelVh(:,:,jh)
! -------------------------------------------------------------------- 
!       INTEGRATION EXTREMA COMPUTATION
! ------------------- t(1) ----------------------                         
  t(1) = 0.d0 
! ------------------- t(2) ---------------------- 
  ct(2) = sigma*(pig+p2-lbar2) + tau*(pig+p1-lbar1) 
  st(2) = rho*(pig+p2-lbar2) 
  t(2)  = datan2(st(2),ct(2)) 
  IF (t(2).lt.0.d0) THEN 
     WRITE(*,*)'ERROR in t(2)!' 
  ENDIF
! ------------------- t(3) ---------------------- 
  ct(3) = sigma*(pig+p2-lbar2) - tau*(pig-p1+lbar1) 
  st(3) = rho*(pig+p2-lbar2) 
  t(3)  = datan2(st(3),ct(3)) 
  IF (t(3).lt.0.d0) THEN 
     WRITE(*,*)'ERROR in t(3)!' 
  ENDIF
! ------------------- t(4) ---------------------- 
  ct(4) = -sigma*(pig-p2+lbar2) - tau*(pig-p1+lbar1) 
  st(4) = -rho*(pig-p2+lbar2) 
  t(4)  = datan2(st(4),ct(4)) 
  IF (t(4).gt.0.d0) THEN 
     WRITE(*,*)'ERROR in t(4)!' 
  ENDIF
  t(4) = t(4) + dpig 
! checking angle between t(3) and t(4)                              
  IF((t(4)-t(3)).gt.pig) THEN 
     WRITE(*,*)'ERROR in DIFFERENCE t(4)-t(3)!' 
  ENDIF
! ------------------- t(5) ---------------------- 
  ct(5) = -sigma*(pig-p2+lbar2) + tau*(pig+p1-lbar1) 
  st(5) = -rho*(pig-p2+lbar2)
  t(5)  = datan2(st(5),ct(5)) 
  IF (t(5).gt.0.d0) THEN 
     WRITE(*,*)'ERROR in t(5)!' 
  ENDIF
  t(5) = t(5) + dpig 
! checking angle between t(5) and t(1)                              
  IF((t(5)-t(1)).lt.pig) THEN 
     WRITE(*,*)'ERROR in DIFFERENCE t(1)-t(5)!' 
  ENDIF
! ------------------- t(6) ---------------------- 
  t(6)= dpig 
! -----------------------------------------------------------------
! initialization
  DO k=1,4
     ddd(k)=0.d0 
     eee(k)=0.d0 
  nnn(k)=0 
  ENDDO
! preparing dqags call
  epsabs=1.d-8 
  epsrel=1.d-5 
  limit=limx 
  lenw=limx4 
! ----------------------------------------------------------------- 
!          APPROXIMATE RIGHT HAND SIDE COMPUTATION                           
! 1=G-deriv, 2=Z-deriv, 3=g-deriv, 4=z-deriv

  DO j=1,4
     ddel = j
                                                
     CALL dqags(ff1,t(1),t(2),epsabs,epsrel,rm(1),abserr(1), &
          & neval(1),ier(1),limit,lenw,last,iwork,work)        
!     ======================================================
     CALL dqags(ff2,t(2),t(3),epsabs,epsrel,rm(2),abserr(2), &
          & neval(2),ier(2),limit,lenw,last,iwork,work)        
!     ======================================================
      CALL dqags(ff3,t(3),t(4),epsabs,epsrel,rm(3),abserr(3), &
          & neval(3),ier(3),limit,lenw,last,iwork,work)        
!     ======================================================   
     CALL dqags(ff4,t(4),t(5),epsabs,epsrel,rm(4),abserr(4), &
          & neval(4),ier(4),limit,lenw,last,iwork,work)        
!     ====================================================== 
     CALL dqags(ff5,t(5),t(6),epsabs,epsrel,rm(5),abserr(5), &
          & neval(5),ier(5),limit,lenw,last,iwork,work)        
! ====================================================================
! SUM of the RESULTS for each call                                  
     DO  k = 1,5 
        ddd(j) = ddd(j) + rm(k)*(gmp/dpig**2) 
        eee(j) = eee(j) + abserr(k)*(gmp/dpig**2) 
        nnn(j) = nnn(j) + neval(k) 
     ENDDO

! COMMENTED ONLY FOR A TEST !!! (please UNCOMMENT!!!!)
! Derivative with respect to a generic Delaunay element
     ddd(j)=-(gmp/dpig)*(Delta*ddmint2(j)/sqrt(dmint2) + &
          & dDelDelta(j)*sqrt(dmint2)) + ddd(j)

  ENDDO

END SUBROUTINE approxrhs2

! ===========================================================       
!  FUNCTIONS defining the right-hand side of the equations
! ===========================================================       
  DOUBLE PRECISION FUNCTION ff1(t) 
    REAL(KIND=dkind),INTENT(IN) :: t
! ------------------ end interface ---------------------------
    REAL(KIND=dkind) :: r1,dr1 

! INTEGRAND FUNCTION                                                
    r1 = rho*tau*(pig+p1-lbar1)/(rho*cos(t)-sigma*sin(t)) 
                                                                        
    dr1 = ((dDeltau(ddel)*rho + tau*dDelrho(ddel))*(pig+p1-lbar1) - &
         & tau*rho*dlbar(1,ddel))/(rho*cos(t)-sigma*sin(t)) - &
         & rho*tau*(pig+p1-lbar1)*(dDelrho(ddel)*cos(t) - &
         & dDelsigma(ddel)*sin(t))/((rho*cos(t)-sigma*sin(t))**2) 

    ff1 = dDelDelta(ddel)*sqrt(dmint2+r1**2) + &
         & Delta/sqrt(dmint2+r1**2)*(ddmint2(ddel)+r1*dr1)

    RETURN 
  END FUNCTION ff1
                                                                        
! ******************************************************************
  DOUBLE PRECISION FUNCTION ff2(t) 
    REAL(KIND=dkind),INTENT(IN) :: t
! ------------------ end interface --------------------------------
    REAL(KIND=dkind) :: r2,dr2 
! ================================================================= 
! INTEGRAND FUNCTION                                                
    r2 = rho*(pig+p2-lbar2)/sin(t) 
                                                                        
    dr2 = (dDelrho(ddel)*(pig+p2-lbar2) - rho*dlbar(2,ddel))/sin(t)
    
    ff2 = dDelDelta(ddel)*dsqrt(dmint2+r2**2) + &
         & Delta/sqrt(dmint2+r2**2)*(ddmint2(ddel)+r2*dr2)

    RETURN 
  END FUNCTION ff2
      
! ******************************************************************
  DOUBLE PRECISION FUNCTION ff3(t) 
    REAL(KIND=dkind),INTENT(IN) :: t
! ------------------ end interface --------------------------------
    REAL(KIND=dkind) :: r3,dr3 
! ================================================================= 
! INTEGRAND FUNCTION                                                
    r3 = -rho*tau*(pig-p1+lbar1)/(rho*cos(t)-sigma*sin(t)) 
                                                                        
    dr3 = -((dDeltau(ddel)*rho + tau*dDelrho(ddel))*(pig-p1+lbar1) + &
         & tau*rho*dlbar(1,ddel))/(rho*cos(t)-sigma*sin(t)) + &
         & tau*rho*(pig-p1+lbar1)*(dDelrho(ddel)*cos(t) - &
         & dDelsigma(ddel)*sin(t))/((rho*cos(t)-sigma*sin(t))**2)          
    
    ff3 = dDelDelta(ddel)*sqrt(dmint2+r3**2) + &
         & Delta/sqrt(dmint2+r3**2)*(ddmint2(ddel)+r3*dr3)   
    
    RETURN 
  END FUNCTION ff3
                                                                        
! ******************************************************************
  DOUBLE PRECISION FUNCTION ff4(t) 
    REAL(KIND=dkind),INTENT(IN) :: t
! ------------------ end interface -------------------------------- 
    REAL(KIND=dkind) :: r4,dr4 
! ================================================================= 
! INTEGRAND FUNCTION                                                
    r4 = -rho*(pig-p2+lbar2)/sin(t) 
    
    dr4 = -(dDelrho(ddel)*(pig-p2+lbar2) + rho*dlbar(2,ddel))/sin(t)
     
    ff4 = dDelDelta(ddel)*sqrt(dmint2+r4**2) + &
         & Delta/sqrt(dmint2+r4**2)*(ddmint2(ddel)+r4*dr4)   
 
    RETURN 
  END FUNCTION ff4
                                                                        
! ******************************************************************
  DOUBLE PRECISION FUNCTION ff5(t) 
    REAL(KIND=dkind),INTENT(IN) :: t
! ------------------ end interface -------------------------------- 
    REAL(KIND=dkind) :: r5,dr5 
! ================================================================= 
! INTEGRAND FUNCTION                                                
    r5 = tau*rho*(pig+p1-lbar1)/(rho*cos(t)-sigma*sin(t)) 
                                                                        
    dr5 = ((dDeltau(ddel)*rho + tau*dDelrho(ddel))*(pig+p1-lbar1) - &
         & tau*rho*dlbar(1,ddel))/(rho*cos(t)-sigma*sin(t)) - &
         & rho*tau*(pig+p1-lbar1)*(dDelrho(ddel)*cos(t) - &
         & dDelsigma(ddel)*sin(t))/((rho*cos(t)-sigma*sin(t))**2)          
                                                                        
    ff5 = dDelDelta(ddel)*sqrt(dmint2+r5**2) + &
         & Delta/sqrt(dmint2+r5**2)*(ddmint2(ddel)+r5*dr5)   

    RETURN 
  END FUNCTION ff5

! =================================================================
! COMPUTE ASTEROID AND PLANET ROTATIONAL MATRICES
! input : asteroid elems (elem), plantet elems (elpl)
! output: aa,ea,beta,beta2,ci,si;apl,epl,betap,betap2;
!         R_om,R_i,R_on,dR_om,dR_i,dR_on,Hast,dHdI,dHdom,dHdon,R_omp,R_ip,Hp
! =================================================================
SUBROUTINE rot_matrices(elpl,elem)
  TYPE(orbit_elem),INTENT(IN) :: elem,elpl !planet and asteroid elems
! end interface
  REAL(KIND=dkind) :: com,som,con,son,omaux
  REAL(KIND=dkind) :: omp,ip,onod,comp,somp,cip,sip
  REAL(KIND=dkind),DIMENSION(3,3) :: maux

! Planet elements
  apl = elpl%coord(1)
  epl = elpl%coord(2)
  omp = elpl%coord(5) 
  ip  = elpl%coord(3)
  betap2 = 1.d0-epl**2
  betap  = sqrt(betap2)
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
! Asteroid elements
  aa = elem%coord(1)
  ea = elem%coord(2)
  beta2 = 1.d0-ea**2 
  beta = sqrt(beta2)
! Rotation matrix of angle omega for the asteroid
  com = cos(elem%coord(5))
  som = sin(elem%coord(5))
  R_om = 0.d0
  R_om(1,1) = com
  R_om(1,2) =-som
  R_om(2,1) = som
  R_om(2,2) = com
  R_om(3,3) = 1.d0
! Derivatives of the rotation matrix of angle omega for the asteroid
  dR_om = 0.d0
  dR_om(1,1) =-som
  dR_om(1,2) =-com
  dR_om(2,1) = com
  dR_om(2,2) =-som
! Compute cosine of inclination from Z/G for the asteroid
  ci = cos(elem%coord(3)) !ci = zl/g  !old ci=zl/(ky*sqrt(a)*beta) 
  IF(ci.lt.1.d0) THEN 
     si = sqrt(1.d0-ci**2) 
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
! Derivatives of the rotation matrix of angle Inclination for the asteroid
  dR_i = 0.d0
  dR_i(2,2) =-si
  dR_i(2,3) =-ci
  dR_i(3,2) = ci
  dR_i(3,3) =-si
! Angle difference between Omega.nod of the asteroid and Omega.nod of the planet
  onod = elem%coord(4)-elpl%coord(4)
  con = cos(onod)
  son = sin(onod)
! Rotation matrix of angle Omnod for the asteroid
  R_on = 0.d0
  R_on(1,1) = con
  R_on(1,2) =-son
  R_on(2,1) = son
  R_on(2,2) = con
  R_on(3,3) = 1.d0
! Derivatives of the rotation matrix of angle Omnod for the asteroid
  dR_on = 0.d0
  dR_on(1,1) =-son
  dR_on(1,2) =-con
  dR_on(2,1) = con
  dR_on(2,2) =-son

  CALL prodmm(maux,R_on,R_i)
  CALL prodmm(Hast,maux,R_om)  !Hast = R_on x R_i x R_om

  CALL prodmm(Hp,R_ip,R_omp)   !Hp = R_ip x R_om

  CALL prodmm(maux,R_on,dR_i)
  CALL prodmm(dHdI,maux,R_om)  !dHdI: derivative of Hast w.r.t. inclination

  CALL prodmm(maux,R_on,R_i)
  CALL prodmm(dHdom,maux,dR_om)!dHdom: derivative of Hast w.r.t. omega

  CALL prodmm(maux,dR_on,R_i)
  CALL prodmm(dHdon,maux,R_om) !dHdon: derivative of Hast w.r.t. Omega.nod

END SUBROUTINE rot_matrices

! =================================================================
! COMPUTE  SINGULAR ECCENTRIC ANOMALY VALUES
! input : pla and ast singular true anomalies (f1h,f2h)
!         pla and ast eccentricities (epl,ea)
! output: u1sgn,u2sgn and their cosine and sine
!         lbar1,lbar2,dl1,dl2
! =================================================================
SUBROUTINE anomalies_at_min(v1min,v2min,u1sgn,u2sgn,t1,t2)
! ------------- interface ------------
  REAL(KIND=dkind),DIMENSION(nminx),INTENT(IN) :: v1min,v2min
  REAL(KIND=dkind),DIMENSION(nminx),INTENT(OUT) :: u1sgn,u2sgn
  REAL(KIND=dkind),DIMENSION(nminx),INTENT(OUT) :: t1,t2 !angle traslations
! ----------- end interface -----------
  REAL(KIND=dkind) :: f1,f2,sf1,sf2,cf1,cf2,cu1,su1,cu2,su2,chusgn
  INTEGER :: j
! -------------------------------------

  IF(nmin.eq.1)THEN
     
     f1 = v1min(1)
     f2 = v2min(1)

     cf1 = cos(f1)
     sf1 = sin(f1)
     cf2 = cos(f2)
     sf2 = sin(f2)

! Eccentric and mean anomalies values at loc min points
     cu1 = (epl+cf1)/(1.d0+epl*cf1)
     su1 = sqrt(1.d0-(cu1**2))*SIGN(1.d0,sf1)
     u1sgn(1) = ATAN2(su1,cu1)
     IF(ABS(u1sgn(1)).gt.3.d0) THEN
!     WRITE(*,*) 'change planet ecc anomaly in [O,dpig] range'
        CALL rad02pi(u1sgn(1),chusgn) 
        u1sgn(1) = chusgn
        t1(1) = pig
!        write(*,*) cu1,cos(u1sgn),u1sgn
     ELSE
        t1(1) = 0.d0
     ENDIF

     cu2 = (ea+cf2)/(1+ea*cf2)
     su2 = sqrt(1.d0-(cu2**2))*SIGN(1.d0,sf2)
     u2sgn(1) = ATAN2(su2,cu2)
     IF(ABS(u2sgn(1)).gt.3.d0) THEN
!     WRITE(*,*) 'change asteroid ecc anomaly in [0,dpig] range'
        CALL rad02pi(u2sgn(1),chusgn) 
        u2sgn(1) = chusgn
        t2(1) = pig
!        write(*,*) cu2,cos(u2sgn),u2sgn
     ELSE
        t2(1) = 0.d0
     ENDIF

  ELSE
     DO j=1,nmin

        f1 = v1min(j)
        f2 = v2min(j)

        cf1 = cos(f1)
        sf1 = sin(f1)
        cf2 = cos(f2)
        sf2 = sin(f2)

! Eccentric and mean anomalies values at loc min points
        cu1 = (epl+cf1)/(1.d0+epl*cf1)
        su1 = sqrt(1.d0-(cu1**2))*SIGN(1.d0,sf1)
        u1sgn(j) = ATAN2(su1,cu1)
!        IF(ABS(u1sgn(j)).gt.3.d0) THEN
!     WRITE(*,*) 'change planet ecc anomaly in [O,dpig] range'
!        CALL rad02pi(u1sgn(j),chusgn) 
!        u1sgn(j) = chusgn
!        t1(j) = pig
!        write(*,*) cu1,cos(u1sgn),u1sgn
!     ELSE
        t1(j) = 0.d0
!     ENDIF

        cu2 = (ea+cf2)/(1+ea*cf2)
        su2 = sqrt(1.d0-(cu2**2))*SIGN(1.d0,sf2)
        u2sgn(j) = ATAN2(su2,cu2)
!        IF(ABS(u2sgn(j)).gt.3.d0) THEN
!     WRITE(*,*) 'change asteroid ecc anomaly in [0,dpig] range'
        CALL rad02pi(u2sgn(j),chusgn) 
        u2sgn(j) = chusgn
        t2(j) = pig
!        write(*,*) cu2,cos(u2sgn),u2sgn
!     ELSE
!        t2(1) = 0.d0
!     ENDIF
     ENDDO
  ENDIF

END SUBROUTINE anomalies_at_min

! ======================================
! Check extrema of integration interval
! ======================================
SUBROUTINE int_interval(t1,t2,nummin,pinf,psup,ainf,asup,p1,p2)
  USE fund_const
  USE critical_points
  IMPLICIT NONE
! ----------------- interface ------------------
  REAL(KIND=dkind),DIMENSION(nminx),INTENT(IN) :: t1,t2
  INTEGER,INTENT(IN) :: nummin
  REAL(KIND=dkind),INTENT(OUT) :: pinf,psup,ainf,asup,p1,p2
! --------------- end interface -----------------
  IF(nummin.eq.1)THEN
     pinf =-pig+t1(1)
     psup = pig+t1(1)
     p1 = t1(1)
     ainf =-pig+t2(1)
     asup = pig+t2(1)
     p2 = t2(1)
  ELSE
     pinf =-pig
     psup = pig
     p1 = t1(1)
     ainf = 0.d0
     asup = dpig
     p2 = t2(1)
  ENDIF

END SUBROUTINE int_interval

! =================================================================
! COMPUTE PLANET AND ASTEROID COORD AT SING POINT, 
! TANGENT VECTORS AND MATRIX A_h
! =================================================================
SUBROUTINE coord_tau_Ah(u1,u2,As)

  REAL(KIND=dkind),INTENT(IN) :: u1,u2 !ecc anomalies
  REAL(KIND=dkind),DIMENSION(2,2),INTENT(OUT) :: As !matrix A_h in delta_h
! --------------- end interface -----------------
  REAL(KIND=dkind) :: cu1,su1,cu2,su2
  REAL(KIND=dkind) :: dl1,dl2,du1,ddu1,du2,ddu2
  REAL(KIND=dkind) :: Pi1,Pi2
  REAL(KIND=dkind),DIMENSION(3) :: ddux1,ddux2,ddlx1,ddlx2
  REAL(KIND=dkind),DIMENSION(3) :: y1,duy1,dduy1,dduy2
  REAL(KIND=dkind),DIMENSION(3) :: tang1,tang2 !tangent vectors
! ------------------------------------------------
  cu1 = cos(u1)
  su1 = sin(u1)
  cu2 = cos(u2)
  su2 = sin(u2)

! Planet functions
  y1(1) = apl*(cu1-epl)
  y1(2) = apl*betap*su1
  y1(3) = 0.d0
  CALL prodmv(x1,Hp,y1) !planet coord (in the rotated system)
!  write(*,*)'planet coord: ',x1

  duy1(1) =-apl*su1
  duy1(2) = apl*betap*cu1
  duy1(3) = 0.d0
  CALL prodmv(dux1,Hp,duy1) !dXpl/du in u=u1sgn

! second derivatives 
  dduy1(1) =-apl*cu1
  dduy1(2) =-apl*betap*su1
  dduy1(3) = 0.d0
  CALL prodmv(ddux1,Hp,dduy1) !ddXpl/ddu

! Asteroid functions
  y2(1) = aa*(cu2-ea)
  y2(2) = aa*beta*su2
  y2(3) = 0.d0
  CALL prodmv(x2,Hast,y2) !asteroid coord (in the rotated system)

  duy2(1) =-aa*su2 
  duy2(2) = aa*beta*cu2
  duy2(3) = 0.d0
  CALL prodmv(dux2,Hast,duy2) !dXa/du in u=u2sgn

! second derivatives
  dduy2(1) =-aa*cu2 
  dduy2(2) =-aa*beta*su2 
  dduy2(3) = 0.d0
  CALL prodmv(ddux2,Hast,dduy2) !ddXa/ddu

! --------------------------------------------------------------------
! Tangent vectors at lbar1,lbar2
  tang1(1:3) = dux1/(1.d0-epl*cu1) 
  tang2(1:3) = dux2/(1.d0-ea*cu2)

! --------------------------------------------------------------------
! dot products giving \Pi_j(\calE)
  dl1 = 1.d0-epl*cu1
  dl2 = 1.d0-ea*cu2
  du1   = 1.d0/dl1
  ddu1  = -epl*su1*du1**3
  ddlx1 = ddux1*du1**2 + dux1*ddu1 
  du2 = 1.d0/dl2
  ddu2 = -ea*su2*du2**3
  ddlx2 = ddux2*du2**2 + dux2*ddu2 

  Pi1 = DOT_PRODUCT(ddlx1,x1-x2)
  Pi2 = DOT_PRODUCT(ddlx2,x1-x2)
! -------------------------------------------------------------------
! Matrix A_h at singular  points (lbar1,lbar2)
  As(1,1) = DOT_PRODUCT(tang1,tang1) + Pi1
  As(2,2) = DOT_PRODUCT(tang2,tang2) - Pi2
  As(1,2) =-DOT_PRODUCT(tang1,tang2)
  As(2,1) = As(1,2)

END SUBROUTINE coord_tau_Ah

! =================================================================
! COMPUTE DERIVATIVES W.R.T. DEL-ELEMS AT SING POINTS OF
! X1,X2,DMINTIL^2,TAU1,TAU2,A_h,V_h
! =================================================================
SUBROUTINE compute_derivs(u1,u2,dV,dA)
! ------------------- interface ----------------------
  REAL(KIND=dkind),INTENT(IN) :: u1,u2
  REAL(KIND=dkind),DIMENSION(2,4),INTENT(OUT) :: dV
  REAL(KIND=dkind),DIMENSION(2,2,4),INTENT(OUT) :: dA
! ------------------ end interface --------------------
  REAL(KIND=dkind) :: cu1,su1,cu2,su2,dl1,dl2
  REAL(KIND=dkind) :: du1,du2,ddu1,ddu2,dddu1,dddu2
! Variables for derivatives
  REAL(KIND=dkind),DIMENSION(3) :: tang1,tang2
  REAL(KIND=dkind),DIMENSION(3) :: dduy1,dduy2,dey2,deduy2
  REAL(KIND=dkind),DIMENSION(3) :: ddux1,ddux2,dex2,dedux2
  REAL(KIND=dkind),DIMENSION(3) :: ddlx1,ddlx2
  REAL(KIND=dkind),DIMENSION(3) :: ddduy1,ddduy2,dddux1,dddux2
  REAL(KIND=dkind),DIMENSION(2,2) :: Hess,InvHess
  REAL(KIND=dkind) :: detH
  REAL(KIND=dkind),DIMENSION(3) :: dIx2,domx2,donx2,dIdux2,domdux2,dondux2
  REAL(KIND=dkind),DIMENSION(3) :: dedduy2,deddux2
  REAL(KIND=dkind),DIMENSION(3) :: dIddux2,domddux2,donddux2
  REAL(KIND=dkind),DIMENSION(2) :: dedud2,dIdud2,domdud2,dondud2
  REAL(KIND=dkind),DIMENSION(2,4) :: dKepUh,dKepVh
! DEL derivatives
  REAL(KIND=dkind),DIMENSION(3,4) :: dKeptang1,dKeptang2,dtang1,dtang2
  REAL(KIND=dkind),DIMENSION(3,4) :: dKepx2,dKepdux2,dKepddux2
  REAL(KIND=dkind),DIMENSION(3,4) :: dKepsum1,dKepsum2
  REAL(KIND=dkind),DIMENSION(3,4) :: dKepDeltah
  REAL(KIND=dkind),DIMENSION(1,4) :: dKepPi1,dKepPi2
  REAL(KIND=dkind),DIMENSION(1,4) :: dPi1,dPi2

  INTEGER :: k
! auxliar variables
  REAL(KIND=dkind),DIMENSION(3) :: aux1,aux2
! ======================================================================
  cu1 = cos(u1)
  su1 = sin(u1)
  cu2 = cos(u2)
  su2 = sin(u2)

! dell/du
  dl1 = 1.d0-epl*cu1
  dl2 = 1.d0-ea*cu2

! --------------------------------------------------------------------
! Tangent vectors at lbar1,lbar2
  tang1(1:3) = dux1/dl1
  tang2(1:3) = dux2/dl2

! Derivatives of coord vectors (at loc min) w.r.t. KEP,eccentric anomalies
! planet
  dduy1(1) =-apl*cu1
  dduy1(2) =-apl*betap*su1
  dduy1(3) = 0.d0
  CALL prodmv(ddux1,Hp,dduy1) !ddXpl/ddu

  ddduy1(1) = apl*su1
  ddduy1(2) =-apl*betap*cu1
  ddduy1(3) = 0.d0
  CALL prodmv(dddux1,Hp,ddduy1) !dddXpl/dddu

! asteroid
  dduy2(1) =-aa*cu2 
  dduy2(2) =-aa*beta*su2 
  dduy2(3) = 0.d0
  CALL prodmv(ddux2,Hast,dduy2) !ddXa/ddu

  ddduy2(1) = aa*su2 
  ddduy2(2) =-aa*beta*cu2 
  ddduy2(3) = 0.d0
  CALL prodmv(dddux2,Hast,ddduy2) !dddXa/dddu

  dey2(1) =-aa 
  dey2(2) =-aa*ea/beta*su2
  dey2(3) = 0.d0
  CALL prodmv(dex2,Hast,dey2) !dXa/de
  deduy2(1) = 0.d0
  deduy2(2) =-aa*ea/beta*cu2
  deduy2(3) = 0.d0
  CALL prodmv(dedux2,Hast,deduy2) !ddXa/dedu
  dedduy2(1) = 0.d0
  dedduy2(2) = -dey2(2)
  dedduy2(3) = 0.d0
  CALL prodmv(deddux2,Hast,dedduy2) !dddXa/deddu

  CALL prodmv(dIx2,dHdI,y2)       !dXa/dI
  CALL prodmv(dIdux2,dHdI,duy2)   !ddXa/dIdu
  CALL prodmv(dIddux2,dHdI,dduy2) !dddXa/dIddu

  CALL prodmv(donx2,dHdon,y2)       !dXa/don
  CALL prodmv(dondux2,dHdon,duy2)   !ddXa/dondu
  CALL prodmv(donddux2,dHdon,dduy2) !dddXa/donddu

  CALL prodmv(domx2,dHdom,y2)       !dXa/dom
  CALL prodmv(domdux2,dHdom,duy2)   !ddXa/domdu
  CALL prodmv(domddux2,dHdom,dduy2) !dddXa/domddu

  dKepx2(1:3,1) = dex2
  dKepx2(1:3,2) = dIx2
  dKepx2(1:3,3) = donx2
  dKepx2(1:3,4) = domx2

  dKepdux2(1:3,1) = dedux2
  dKepdux2(1:3,2) = dIdux2
  dKepdux2(1:3,3) = dondux2
  dKepdux2(1:3,4) = domdux2

  dKepddux2(1:3,1) = deddux2
  dKepddux2(1:3,2) = dIddux2
  dKepddux2(1:3,3) = donddux2
  dKepddux2(1:3,4) = domddux2

! --------------------------------------------------------------------
! Hessian matrix of d2=<Xp-Xa,Xp-Xa> and its inverse at u1sgn,u2sgn
  Hess(1,1) = 2.d0*(DOT_PRODUCT(dux1,dux1)+DOT_PRODUCT((x1-x2),ddux1))
  Hess(2,2) = 2.d0*(DOT_PRODUCT(dux2,dux2)-DOT_PRODUCT((x1-x2),ddux2))
  Hess(1,2) = -2.d0*DOT_PRODUCT(dux1,dux2)
  Hess(2,1) = Hess(1,2)
  CALL inv22(Hess,InvHess,detH)
  IF(detH.le.0.d0)THEN 
     write(*,*)' Kantorovich, Hessian determinat: detH = ',detH
     write(*,*)'InvHess',InvHess
     STOP 
  ENDIF
! ====================================================================
! Derivatives of u1sgn,u2sgn w.r.t. KEP
! d(u1,u2)/de = -InvHess*d(grad d2)/de

! Derivatives of gradient of d2 w.r.t. ecc
  dedud2(1) =-2.d0*DOT_PRODUCT(dex2,dux1)
  dedud2(2) = 2.d0*DOT_PRODUCT(dex2,dux2) - &
       & 2.d0*DOT_PRODUCT((x1-x2),dedux2)

  CALL prodmv22(dKepUh(1:2,1),-InvHess,dedud2)  !ecc-deriv
! --------------------------------------------------------------------
! Derivatives of gradient of d2 w.r.t. I
  dIdud2(1) =-2.d0*DOT_PRODUCT(dIx2,dux1)
  dIdud2(2) = 2.d0*DOT_PRODUCT(dIx2,dux2)-2.d0*DOT_PRODUCT((x1-x2),dIdux2)

  CALL prodmv22(dKepUh(1:2,2),-InvHess,dIdud2)  !Incl-deriv
! --------------------------------------------------------------------
! Derivatives of gradient of d2 w.r.t. Omega.nod
  dondud2(1) =-2.d0*DOT_PRODUCT(donx2,dux1)
  dondud2(2) = 2.d0*DOT_PRODUCT(donx2,dux2) - &
       & 2.d0*DOT_PRODUCT((x1-x2),dondux2)

  CALL prodmv22(dKepUh(1:2,3),-InvHess,dondud2) !Omnod-deriv
! --------------------------------------------------------------------
! Derivatives of gradient of d2 w.r.t. omega
  domdud2(1) =-2.d0*DOT_PRODUCT(domx2,dux1)
  domdud2(2) = 2.d0*DOT_PRODUCT(domx2,dux2) - &
       & 2.d0*DOT_PRODUCT((x1-x2),domdux2)

  CALL prodmv22(dKepUh(1:2,4),-InvHess,domdud2) !omega-deriv
! ====================================================================
! KEP-derivatives of vector Vh=(lbar1,lbar2)
  dKepVh(1,1:4) = dl1*dKepUh(1,1:4)
  dKepVh(2,1:4) = dl2*dKepUh(2,1:4)
  dkepVh(2,1) = dKepVh(2,1)-su2 !e-deriv
! G,Z,g,z-derivatives of vector Vh=(lbar1,lbar2)
  dV = MATMUL(dKepVh(1:2,1:4),dKEP_dDEL(1:4,1:4))
! ====================================================================
! Derivative of tang1,tang2 w.r.t. KEP
! here tang1 = dX2/dl2, tang2 = dX2/dl2 are derived at Vh=(l1sgn,l2sgn)
  aux1 = ddux1 - dux1*epl*su1/dl1
  aux2 = ddux2 - dux2*ea*su2/dl2
! e-derivs
  dKeptang1(:,1) = aux1*dKepUh(1,1)/dl1
  dKeptang2(:,1) = (aux2*dKepUh(2,1) + dedux2+dux2*cu2/dl2)/dl2
! I-derivs
  dKeptang1(:,2) = aux1*dKepUh(1,2)/dl1
  dKeptang2(:,2) = (aux2*dKepUh(2,2) + dIdux2)/dl2
! Omnod-derivs
  dKeptang1(:,3) = aux1*dKepUh(1,3)/dl1
  dKeptang2(:,3) = (aux2*dKepUh(2,3) + dondux2)/dl2
! omega-derivs
  dKeptang1(:,4) = aux1*dKepUh(1,4)/dl1
  dKeptang2(:,4) = (aux2*dKepUh(2,4) + domdux2)/dl2
! --------------------------------------------------------------------  
! Derivative of tang1,tang2 w.r.t. DEL
  dtang1 = MATMUL(dKeptang1,dKEP_dDEL)
  dtang2 = MATMUL(dKeptang2,dKEP_dDEL)

! ====================================================================
! Derivative of Pi1,Pi2 w.r.t. KEP

! 1st, 2nd and 3rd derivatives of u wrt ell
  du1 = 1.d0/dl1
  ddu1 = -epl*su1*du1**3
! d3u/dl3 = 3*(du/dl)^5*(d2l/du2)^2 - (du/dl)^4*(d3l/du3)
  dddu1 = 3.d0*du1**5*(epl*su1)**2 - du1**4*(epl*cu1) 

  du2 = 1.d0/dl2
  ddu2 = -ea*su2*du2**3
  dddu2 = 3.d0*du2**5*(ea*su2)**2 - du2**4*(ea*cu2)

  ddlx1 = ddux1*du1**2 + dux1*ddu1
  ddlx2 = ddux2*du2**2 + dux2*ddu2

! derivatives of sum_j = \partial^2\mathcal{X}_j/\partial ell_j^2 
  DO k=1,4
     dKepsum1(1:3,k) = (dddux1(1:3)*du1**3 + 3.d0*ddux1(1:3)*du1*ddu1 + &
          & dux1(1:3)*dddu1)*dKepVh(1,k)
     dKepsum2(1:3,k) = dKepddux2(1:3,k)*du2**2 + dKepdux2(1:3,k)*ddu2 +&
          &            (dddux2(1:3)*du2**3 + 3.d0*ddux2(1:3)*du2*ddu2 + &
          & dux2(1:3)*dddu2)*dKepVh(2,k)
     IF(k.eq.1)THEN
        dKepsum2(1:3,k) = dKepsum2(1:3,k) + 2.d0*ddux2(1:3)*du2*cu2*du2**2  &
             &-dux2(1:3)*su2*(1.d0+2.d0*ea*cu2)*du2**4 
     ENDIF
     dKepDeltah(1:3,k) = -dKepx2(1:3,k) + tang1(1:3)*dKepVh(1,k) - &
          & tang2(1:3)*dKepVh(2,k)

     dKepPi1(1,k) = DOT_PRODUCT(dKepsum1(:,k),x1-x2) + &
          & DOT_PRODUCT(ddlx1,dKepDeltah(:,k))

     dKepPi2(1,k) = DOT_PRODUCT(dKepsum2(:,k),x1-x2) + &
          & DOT_PRODUCT(ddlx2,dKepDeltah(:,k))
  ENDDO

  dPi1(1,:) = MATMUL(dKepPi1(1,:),dKEP_dDEL(:,:))
  dPi2(1,:) = MATMUL(dKepPi2(1,:),dKEP_dDEL(:,:))
! ====================================================================
! Derivatives of Ah(1,1), Ah(2,2) and Ah(1,2) w.r.t. DEL
  DO k=1,4
     dA(1,1,k) = 2.d0*DOT_PRODUCT(dux1/dl1,dtang1(:,k)) + dPi1(1,k) 
     dA(2,2,k) = 2.d0*DOT_PRODUCT(dux2/dl2,dtang2(:,k)) - dPi2(1,k)
     dA(1,2,k) = - DOT_PRODUCT(dux1/dl1,dtang2(:,k)) - &
          & DOT_PRODUCT(dtang1(:,k),dux2/dl2)
     dA(2,1,k) = dA(1,2,k)
  ENDDO

END SUBROUTINE compute_derivs

! ============================================================
! Compute dmintil^2 DEL-derivatives
SUBROUTINE delaunay_ddmintil(dcomdmint)
! ----------------------- inteface ------------------------
  REAL(KIND=dkind),DIMENSION(5,nminx),INTENT(IN) :: dComdmint
! -------------------- end interface ----------------------
  TYPE(orbit_elem) :: com2 ! cometary asteroid elements
  REAL(KIND=dkind),DIMENSION(6,6) :: dCOM_dKEP
  REAL(KIND=dkind),DIMENSION(1,5) :: grad_tmp,vectmp
  REAL(KIND=dkind),DIMENSION(4) :: dDeldmint
  INTEGER :: fail_flag,j
! ---------------------------------------------------------
! Jacobian matrix dKEP_dDEL
  dKEP_dDEL=0.d0
  dKEP_dDEL(4,3) =  1.d0 !domega/dg
  dKEP_dDEL(3,4) =  1.d0 !dOmnod/dz
  dKEP_dDEL(1,1) = -beta/(ky*sqrt(aa)*ea)    !de/dG
  dKEP_dDEL(2,1) =  ci/(si*ky*sqrt(aa)*beta) !dI/dG
  dKEP_dDEL(2,2) = -1.d0/(ky*sqrt(aa)*beta*si) !dI/dZ
!  write(*,*)'dKEP_dDEL :',dKEP_dDEL
! ---------------------------------------------------------
! Jacobian matrix dCOM_dKEP
  CALL coo_cha(elem,'COM',com2,fail_flag,dCOM_dKEP)
  IF(fail_flag.ge.5)THEN
     WRITE(*,*)'error in coo_cha! fail_flag=',fail_flag
  ENDIF
! ---------------------------------------------------------
! Derivatives of dmintil**2 w.r.t. DEL
  DO j=1,nmin
     vectmp(1,1:5) = dcomdmint(1:5,j)
     grad_tmp(1,1:5) = MATMUL(vectmp(1,1:5),dCOM_dKEP(1:5,1:5))
     dDeldmint  = MATMUL(grad_tmp(1,2:5),dKEP_dDEL(1:4,1:4))
     dDeldmint2(:,j) = dmintil(j)*dDeldmint !actually it is half the derivative
  ENDDO
END SUBROUTINE delaunay_ddmintil


! *********************************************************
! CHECK on derivatives of approximated function
! change the value of i to select G,Z,g=omega,z=Omnod
! *********************************************************
SUBROUTINE check_deriv_approx(deriv)
! ----------------- interface ------------------
  REAL(KIND=dkind),DIMENSION(4),INTENT(IN) :: deriv
! ----------------- end interface --------------
! for dmintil_rms
  INTEGER :: nummin
  REAL(KIND=dkind),DIMENSION(5,nminx) :: ddmintdel2
  REAL(KIND=dkind),DIMENSION(nminx) :: v1min,v2min
  REAL(KIND=dkind),DIMENSION(nminx) :: t1,t2
  REAL(KIND=dkind) :: incr,ddf,ddf_rem,eef
  REAL(KIND=dkind) :: g1,g2,z1,z2
  INTEGER :: nnf,i
! ----------------------------------------------
  write(*,*)'**********************'
  CALL app_pertfun(ddf,eef,nnf)
!  write(*,*)'double integral app pert function:',ddf
!  CALL app_pertfun_bis(1,ddf,eef,nnf)
!  write(*,*)'double integral app pert function (int):',ddf
  ddf_rem=ddf
  incr=1.d-5
  i=1  !no. to select the variable for increment 
  IF(i.eq.1)THEN
     write(*,*)'  check for G-derivative  '
     g1 = ky*sqrt(aa)*beta
     g2 = g1 + incr
     elem%coord(2) = sqrt(1-(g2/ky)**2/aa)
     elem%coord(3) = ACOS(g1*ci/g2)*SIGN(1.d0,si)
  ELSEIF(i.eq.2)THEN
     write(*,*)'  check for Z-derivative  '
     g1 = ky*sqrt(aa)*beta
     z1 = G1*ci
     z2 = z1+incr
     elem%coord(3) = ACOS(z2/g1)*SIGN(1.d0,si)
  ELSEIF(i.eq.3)THEN
     write(*,*)'  check for g-derivative  '
     elem%coord(5) = elem%coord(5) + incr
  ELSEIF(i.eq.4)THEN
     write(*,*)' check for z-derivative  '
     elem%coord(4) = elem%coord(4) + incr
  ENDIF
  CALL rot_matrices(elpl,elem)
  CALL dmintil_rms(elpl,elem,nummin,dmintil,DDMINTDEL2=ddmintdel2, &
       & V1MIN=v1min,V2MIN=v2min)
  CALL anomalies_at_min(v1min,v2min,u1sgn,u2sgn,t1,t2)
  CALL int_interval(t1,t2,nummin,pinf,psup,ainf,asup,p1,p2)
  Vh(1,1) = u1sgn(1)-epl*sin(u1sgn(1))
  Vh(2,1) = u2sgn(1)- ea*sin(u2sgn(1))
  CALL coord_tau_Ah(u1sgn(1),u2sgn(1),Ah(:,:,1))

     CALL app_pertfun(ddf,eef,nnf)
!     write(*,*)'double integral with incremented g:',ddf
!  dmint2 = dmintil(1)**2
!  CALL app_pertfun_bis(1,ddf,eef,nnf)
!  write(*,*)'double integral app pert function (int):',ddf

  write(*,*)'derivative       ',deriv(i)
  write(*,*)'incremental ratio',(ddf-ddf_rem)/incr
  write(*,*)'**********************'
  STOP     

END SUBROUTINE check_deriv_approx


! *********************************************************
SUBROUTINE checkder(deriv,i)
! the order of the variables is: G,Z,g,z
! ----------------- interface ------------------
  REAL(KIND=dkind),DIMENSION(4),INTENT(IN) :: deriv
  INTEGER,INTENT(IN) :: i
! ----------------- end interface --------------
! for dmintil_rms
  INTEGER :: nummin
  REAL(KIND=dkind),DIMENSION(5,nminx) :: ddmintdel2
  REAL(KIND=dkind),DIMENSION(nminx) :: v1min,v2min
  REAL(KIND=dkind),DIMENSION(nminx) :: t1,t2
  REAL(KIND=dkind) :: incr,ddf,ddf_rem,eef
  REAL(KIND=dkind) :: g1,g2,z1,z2
  INTEGER :: nnf
! ----------------------------------------------
  write(*,*)'**********************'
!  CALL app_pertfun(ddf,eef,nnf)
!  write(*,*)'double integral app pert function:',ddf
!  CALL app_pertfun_bis(1,ddf,eef,nnf)
!  write(*,*)'double integral app pert function (int):',ddf
!  ddf_rem=ddf

  ddf_rem = Ah(2,2,1)

  incr=1.d-6
!  i=j  !no. to select the variable for increment 
  IF(i.eq.1)THEN
     write(*,*)'  check for G-derivative  '
     g1 = ky*sqrt(aa)*beta
     g2 = g1 + incr
     elem%coord(2) = sqrt(1-(g2/ky)**2/aa)
     elem%coord(3) = ACOS(g1*ci/g2)*SIGN(1.d0,si)
  ELSEIF(i.eq.2)THEN
     write(*,*)'  check for Z-derivative  '
     g1 = ky*sqrt(aa)*beta
     z1 = G1*ci
     z2 = z1+incr
     elem%coord(3) = ACOS(z2/g1)*SIGN(1.d0,si)
  ELSEIF(i.eq.3)THEN
     write(*,*)'  check for g-derivative  '
     elem%coord(5) = elem%coord(5) + incr
  ELSEIF(i.eq.4)THEN
     write(*,*)' check for z-derivative  '
     elem%coord(4) = elem%coord(4) + incr
  ENDIF
  CALL rot_matrices(elpl,elem)
  CALL dmintil_rms(elpl,elem,nummin,dmintil,DDMINTDEL2=ddmintdel2, &
       & V1MIN=v1min,V2MIN=v2min)
  CALL anomalies_at_min(v1min,v2min,u1sgn,u2sgn,t1,t2)
  CALL int_interval(t1,t2,nummin,pinf,psup,ainf,asup,p1,p2)
  Vh(1,1) = u1sgn(1)-epl*sin(u1sgn(1))
  Vh(2,1) = u2sgn(1)- ea*sin(u2sgn(1))
  CALL coord_tau_Ah(u1sgn(1),u2sgn(1),Ah(:,:,1))

!     CALL app_pertfun(ddf,eef,nnf)
!     write(*,*)'double integral with incremented g:',ddf
!  dmint2 = dmintil(1)**2
!  CALL app_pertfun_bis(1,ddf,eef,nnf)
!  write(*,*)'double integral app pert function (int):',ddf

  ddf = Ah(2,2,1)

  write(*,*)'derivative       ',deriv(i)
  write(*,*)'incremental ratio',(ddf-ddf_rem)/incr
  write(*,*)'**********************'
  STOP     

END SUBROUTINE checkder

! ******************************************************************
! APPROXIMATE RIGHT-HAND SIDE (using polar coordinates)
! INTEGRATION DOMAIN: T = [PINF,PSUP] x [AINF,ASUP]
SUBROUTINE app_pertfun_bis(jh,ddd,eee,nnn)
! ===================== INTERFACE ================================= 
  INTEGER,INTENT(IN) :: jh !index of current min point 
! right hand side,errors,number eval                        
  REAL(KIND=dkind),INTENT(OUT) :: ddd,eee 
  INTEGER,INTENT(OUT) :: nnn
! ======================= END INTERFACE ===========================
  REAL(KIND=dkind) :: vsize,vtau1,detAh!,lbar1,lbar2
! integration extrema                                               
  REAL(KIND=dkind),DIMENSION(6) :: ct,st,t 
! for dqags                                                         
  INTEGER,PARAMETER :: limx=300,limx4=4*limx
  INTEGER :: ier(5),limit,iwork(limx),lenw,last 
! function evaluations                                              
  INTEGER,DIMENSION(5) :: neval
  REAL(KIND=dkind) :: epsabs,epsrel,abserr(5),work(limx4) 
! output of dqags                                                   
  REAL(KIND=dkind),DIMENSION(5) :: rm
! loop indexes                                                      
  INTEGER :: i,j,k,ider,fail_flag
! ==================================================================
! Determinant of matrix A, Delta
  detAh = Ah(1,1,jh)*Ah(2,2,jh)-Ah(1,2,jh)**2
  IF(detAh.le.0.d0)THEN 
     write(*,*)' APP_PERTFUN_BIS, detA=',detAh 
     STOP 
  ENDIF
  Delta = 1.d0/(sqrt(detAh))
  vtau1 = sqrt(Ah(1,1,jh))
! --------------------------------------------------------------------
!rho = rho4 ; sigma = rho3 ; tau = rho1
  rho = 1.d0/(Delta*vtau1)
  sigma = Ah(1,2,jh)/vtau1
  tau = vtau1
! --------------------------------------------------------------------
  lbar1 = Vh(1,jh)
  lbar2 = Vh(2,jh)
! ********************************************************************
!     INTEGRATION EXTREMA COMPUTATION
! ********************************************************************
!     %%%%%%%%%%%%%%%%% t(1) %%%%%%%%%%%%%%%%                           
  t(1) = 0.d0 
!     %%%%%%%%%%%%%%%%% t(2) %%%%%%%%%%%%%%%%                           
  ct(2) = sigma*(pig+p2-lbar2) + tau*(pig+p1-lbar1) 
  st(2) = rho*(pig+p2-lbar2) 
  t(2)  = datan2(st(2),ct(2)) 
  IF (t(2).lt.0.d0) THEN 
     WRITE(*,*)'ERROR in t(2)!' 
  ENDIF
!     %%%%%%%%%%%%%%%%% t(3) %%%%%%%%%%%%%%%%                           
  ct(3) = sigma*(pig+p2-lbar2) - tau*(pig-p1+lbar1) 
  st(3) = rho*(pig+p2-lbar2) 
  t(3)  = datan2(st(3),ct(3)) 
  IF (t(3).lt.0.d0) THEN 
     WRITE(*,*)'ERROR in t(3)!' 
  ENDIF
!     %%%%%%%%%%%%%%%%% t(4) %%%%%%%%%%%%%%%%                           
  ct(4) = -sigma*(pig-p2+lbar2) - tau*(pig-p1+lbar1) 
  st(4) = -rho*(pig-p2+lbar2) 
  t(4)  = datan2(st(4),ct(4)) 
  IF (t(4).gt.0.d0) THEN 
     WRITE(*,*)'ERROR in t(4)!' 
  ENDIF
  t(4) = t(4) + dpig 
!     checking angle between t(3) and t(4)                              
  IF((t(4)-t(3)).gt.pig) THEN 
     WRITE(*,*)'ERROR in DIFFERENCE t(4)-t(3)!' 
  ENDIF
!     %%%%%%%%%%%%%%%%% t(5) %%%%%%%%%%%%%%%%                           
  ct(5) = -sigma*(pig-p2+lbar2) + tau*(pig+p1-lbar1) 
  st(5) = -rho*(pig-p2+lbar2)
  t(5)  = datan2(st(5),ct(5)) 
  IF (t(5).gt.0.d0) THEN 
     WRITE(*,*)'ERROR in t(5)!' 
  ENDIF
  t(5) = t(5) + dpig 
!     checking angle between t(5) and t(1)                              
  IF((t(5)-t(1)).lt.pig) THEN 
     WRITE(*,*)'ERROR in DIFFERENCE t(1)-t(5)!' 
  ENDIF
!     %%%%%%%%%%%%%% t(6) %%%%%%%%%%%%%%%%%%                            
  t(6)= dpig 
!     check                                                             
  DO i=1,6                                                          
!     WRITE(*,*)'t(',i,')=',t(i)                                        
  ENDDO
! ---------- initialization --------
  ddd = 0.d0 
  eee = 0.d0 
  nnn = 0 
! ------ preparing dqags call -----
  epsabs=1.d-8 
  epsrel=1.d-5 
  limit=limx 
  lenw=limx4 
! ********************************************************************
!     APPROXIMATE RIGHT HAND SIDE COMPUTATION                           
! ********************************************************************
  CALL dqags(ffun1,t(1),t(2),epsabs,epsrel,rm(1),abserr(1), &
       & neval(1),ier(1),limit,lenw,last,iwork,work)        

  CALL dqags(ffun2,t(2),t(3),epsabs,epsrel,rm(2),abserr(2), &
       & neval(2),ier(2),limit,lenw,last,iwork,work)        

  CALL dqags(ffun3,t(3),t(4),epsabs,epsrel,rm(3),abserr(3), &
       & neval(3),ier(3),limit,lenw,last,iwork,work)        

  CALL dqags(ffun4,t(4),t(5),epsabs,epsrel,rm(4),abserr(4), &
       & neval(4),ier(4),limit,lenw,last,iwork,work)        

  CALL dqags(ffun5,t(5),t(6),epsabs,epsrel,rm(5),abserr(5), &
       & neval(5),ier(5),limit,lenw,last,iwork,work)        

  DO  k = 1,5 
     ddd = ddd + rm(k)
     eee = eee + abserr(k)*(gmp/dpig**2) 
     nnn = nnn + neval(k) 
  ENDDO
  
  ddd = gmp/dpig**2*Delta*(ddd - dpig*sqrt(dmint2))
!  ddd = (gmp/dpig**2)*ddd
END SUBROUTINE app_pertfun_bis

! ===========================================================       
!  FUNCTIONS defining the right-hand side of the equations
! ===========================================================       
DOUBLE PRECISION FUNCTION ffun1(t) 
  REAL(KIND=dkind),INTENT(IN) :: t
! ------------------ end interface --------------------------------
  REAL(KIND=dkind) :: r1,dr1 
! =================================================================
! INTEGRAND FUNCTION                                                
  r1 = rho*tau*(pig+p1-lbar1)/(rho*cos(t)-sigma*sin(t))
  ffun1 = sqrt(dmint2+r1**2)

  RETURN 
END FUNCTION ffun1

! ******************************************************************
DOUBLE PRECISION FUNCTION ffun2(t) 
  REAL(KIND=dkind),INTENT(IN) :: t
! ------------------ end interface --------------------------------
  REAL(KIND=dkind) :: r2,dr2 
! ================================================================= 
! INTEGRAND FUNCTION                                                
  r2 = rho*(pig+p2-lbar2)/sin(t) 
  ffun2 = sqrt(dmint2+r2**2)

  RETURN 
END FUNCTION ffun2
      
! ******************************************************************
DOUBLE PRECISION FUNCTION ffun3(t) 
  REAL(KIND=dkind),INTENT(IN) :: t
! ------------------ end interface --------------------------------
  REAL(KIND=dkind) :: r3,dr3 
! ================================================================= 
! INTEGRAND FUNCTION                                                
  r3 = -rho*tau*(pig-p1+lbar1)/(rho*cos(t)-sigma*sin(t)) 
  ffun3 = sqrt(dmint2+r3**2)

  RETURN 
END FUNCTION ffun3
                                                                        
! ******************************************************************
DOUBLE PRECISION FUNCTION ffun4(t) 
  REAL(KIND=dkind),INTENT(IN) :: t
! ------------------ end interface -------------------------------- 
  REAL(KIND=dkind) :: r4,dr4 
! ================================================================= 
! INTEGRAND FUNCTION                                                
  r4 = -rho*(pig-p2+lbar2)/sin(t) 
  ffun4 = sqrt(dmint2+r4**2)

  RETURN 
END FUNCTION ffun4
                                                                        
! ******************************************************************
DOUBLE PRECISION FUNCTION ffun5(t) 
  REAL(KIND=dkind),INTENT(IN) :: t
! ------------------ end interface -------------------------------- 
  REAL(KIND=dkind) :: r5,dr5 
! =================================================================
! INTEGRAND FUNCTION                                                
  r5 = tau*rho*(pig+p1-lbar1)/(rho*cos(t)-sigma*sin(t)) 
  ffun5 = sqrt(dmint2+r5**2)

  RETURN 
END FUNCTION ffun5

END MODULE kantorovich

! ============================================================
! Multiply a matrix 2x2 by a vector 2x1
SUBROUTINE prodmv22(y,a,x)
  USE fund_const
  IMPLICIT NONE 
  REAL(KIND=dkind),INTENT(IN) :: a(2,2)
  REAL(KIND=dkind) :: x(2),y(2)
  REAL(KIND=dkind) :: s,z(2)
  INTEGER :: j,l 
  z=x                                                                     
  DO  j=1,2 
     s=0.d0 
     DO  l=1,2 
        s=s+a(j,l)*z(l) 
     ENDDO
     y(j)=s 
   ENDDO
 END SUBROUTINE prodmv22

! ============================================================
! Exchange 2nd and 3rd component of a real vector 4-dim
SUBROUTINE riordina(vec)
  USE fund_const
  IMPLICIT NONE
  REAL(KIND=dkind),DIMENSION(4) :: vec
  REAL(KIND=dkind) :: s

  s = vec(2)
  vec(2) = vec(3)
  vec(3) = s

END SUBROUTINE riordina

! ============================================================
! Exchange 2nd and 3rd component of a integer vector 4-dim
SUBROUTINE riordina_int(vec)
  USE fund_const
  IMPLICIT NONE
  INTEGER,DIMENSION(4) :: vec
  INTEGER :: s

  s = vec(2)
  vec(2) = vec(3)
  vec(3) = s

END SUBROUTINE riordina_int
