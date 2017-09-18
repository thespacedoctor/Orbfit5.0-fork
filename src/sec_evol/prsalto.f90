! *********************************************************         
! jump in the derivatives
! now y=(omega,G,Omega_nodale,Z) al tempo t+h
! nc deve diventare un indice sui pianeti da 2 a 6                            
! *********************************************************         
SUBROUTINE prsalto(y,a,elplan,ddg,ddom,ddomnod,ddzl,nc)
  USE orbit_elements
  USE critical_points
  USE kantorovich, ONLY: coord_tau_Ah,rot_matrices!,elpl,elem
  USE fund_const
  IMPLICIT NONE 
! ===================== INTERFACE =======================
  TYPE(orbit_elem),INTENT(IN) :: elplan ! planet orbital elements
  TYPE(orbit_elem) :: elem,elpl !planet and asteroid elems
  REAL(KIND=dkind),INTENT(IN) :: y(4),a
  INTEGER,INTENT(IN) :: nc
! ------------------------------------------------------------
  REAL(KIND=dkind),DIMENSION(nminx),INTENT(OUT) :: ddg,ddom,ddomnod,ddzl
! ======================= END INTERFACE ============================
! Asteroid elements
  REAL(KIND=dkind) aa,e,beta,ci,si
! Varibles in the argument of dmintil_rms
  INTEGER :: nummin
  REAL(KIND=dkind),DIMENSION(nminx) :: dmintil,ddmintil
  REAL(KIND=dkind),DIMENSION(3,nminx) :: c1min,c2min
  REAL(KIND=dkind),DIMENSION(5,nminx) :: ddmintdel2
  REAL(KIND=dkind),DIMENSION(nminx) :: u1min,u2min,v1min,v2min
  REAL(KIND=dkind),DIMENSION(3,nminx) :: tau_1,tau_2
! Auxiliar variables to compute detA
!  REAL(KIND=dkind),DIMENSION(3) :: tau1,tau2
  REAL(KIND=dkind),DIMENSION(2,2) :: Ah
  REAL(KIND=dkind) :: tau1sctau2,deta,vsize
! Variables to convert ddmintdel2 into derivatives w.r.t. KEP
  TYPE(orbit_elem) :: com2 ! asteroid cometary elements
  REAL(KIND=dkind),DIMENSION(6,6) :: dCOM_dKEP,dKEP_dDEL
  INTEGER :: fail_flag
  REAL(KIND=dkind),DIMENSION(1,5) :: ddmintcom2,ddmintkep2,ddmint_del2 ! dmintil derivatives w.r.t. asteroid COM, KEP, DEL elements 
! Varibles for the jump
  REAL(KIND=dkind),DIMENSION(4) :: ddmin_del2 ! dmin derivatives w.r.t. asteroid DEL elements
  REAL(KIND=dkind),DIMENSION(4) :: ddR ! derivatives pert func jump
! ==================================================================  
! loop indexes                                                      
  INTEGER i,j
! ==================================================================
! planet data                                                      
  INCLUDE 'pldata.h90' 
! ******************************************************************

! Compute cosine of inclination from Z/G for the asteroid
  ci=y(4)/y(2)  ! old ci=zl/(ky*sqrt(a)*beta) 
  IF(ci.lt.1.d0)THEN 
     si=dsqrt(1.d0-ci**2) 
  ELSE 
     si=0.d0 
  ENDIF

  elpl = elplan

! asteroid keplerian orbital elements
  elem=undefined_orbit_elem
  elem%coord(1)=a
  elem%coord(2)=sqrt(1.d0-(y(2)/ky)**2/a)
  elem%coord(3)=acos(ci)
  elem%coord(4)=y(3)
  elem%coord(5)=y(1)
  elem%coord(6)=elpl%coord(6) !not used
  elem%coo='KEP'
  elem%t=elpl%t

  aa=a
  e=dsqrt(1.d0-(y(2)/ky)**2/aa)
  beta=dsqrt(1.d0-e**2)                                         
                                                                        
! ********************************************************************     
!                    COMPUTING JUMPING
! ********************************************************************

! -----------------------------------------------------------------------
! computing dmintil, planet and asteroid position at the local minima, 
! their derivatives w.r.t. COM, the vector tau1, tau2 for all local minima
  CALL dmintil_rms(elpl,elem,nummin,dmintil,C1MIN=c1min,&
       & C2MIN=c2min,DDMINTDEL2=ddmintdel2,V1MIN=v1min,V2MIN=v2min)!, &
!       & TAU_1=tau_1,TAU_2=tau_2)
! -----------------------------------------------------------------------

! Jacobian matrix dCOM_dKEP                                                   
  CALL coo_cha(elem,'COM',com2,fail_flag,dCOM_dKEP)

! Jacobian matrix dKEP_dDEL = d(a,e,I,Omnod,om)/d(g,z,G,Z)
  dKEP_dDEL=0.d0
  dKEP_dDEL(5,1)=1.d0
  dKEP_dDEL(4,2)=1.d0
  dKEP_dDEL(2,3)=-beta/(ky*sqrt(aa)*e)
  dKEP_dDEL(3,3)=ci/(si*ky*sqrt(aa)*beta)
  dKEP_dDEL(3,4)=-1.d0/(ky*sqrt(aa)*beta*si)


  CALL rot_matrices(elpl,elem)
  DO i=1,nummin
     CALL ecc_anom(elpl,v1min(i),u1min(i))
     CALL ecc_anom(elem,v2min(i),u2min(i))

     CALL coord_tau_Ah(u1min(i),u2min(i),Ah)
     deta = Ah(1,1)*Ah(2,2)-Ah(1,2)**2
     write(*,*)'1) detA=',deta

     IF(deta.le.0.d0)THEN 
        write(*,*)' deta=',deta,aa,e,y,nc 
        write(*,*)' V_h=',u1min(i),u2min(i)
        write(*,*)' Ah=',Ah(1,1),Ah(2,2),Ah(1,2)
 !       IF(i.eq.1)THEN  !PROBLEM: if we extract the singularity
                        ! at two minima we have to check also for i=2
           STOP 
 !       ENDIF
     ENDIF

! Derivatives of dmintil w.r.t. KEP
     ddmintcom2(1,1:5)=ddmintdel2(1:5,i)
     ddmintkep2(1,1:5)= MATMUL(ddmintcom2(1,1:5),dCOM_dKEP(1:5,1:5))

! Derivative jump of dmintil w.r. to DEL (= (g,z,G,Z))
     ddmint_del2(1,1:4)= MATMUL(ddmintkep2(1,1:5),dKEP_dDEL(1:5,1:4))          

! Derivative jump of dmin_local w.r. to DEL
     ddmin_del2(1:4)=2*ABS(ddmint_del2(1,1:4))

! Derivative jump of the average perturbative function w.r. to DEL
     DO j=1,4
        ddR(j)=gm((nc+1)/2)*ky**2*dpig/sqrt(deta)*ddmin_del2(j)
     ENDDO
                                                                        
! =========================================================         
!     G DERIVATIVE NORMALIZING JUMP      
! =========================================================         
! Derivate jump of of perturbing function with respect to G
!     ddG(i)=ddR(1)/(dpig**2)
     ddG(i)=ddR(3)/(dpig**2)
	
! =========================================================         
!     omega DERIVATIVE NORMALIZING JUMP
! =========================================================         
! Derivate jump of perturbing function with respect to omega      
!     ddom(i)=ddR(4)/(dpig**2)
     ddom(i)=ddR(1)/(dpig**2)

! =========================================================         
!     Z DERIVATIVE NORMALIZING JUMP                                             
! =========================================================         
! Derivate jump of perturbing function with respect to Z     
!     ddzl(i)=ddR(2)/(dpig**2)
     ddzl(i)=ddR(4)/(dpig**2)
	
! =========================================================         
!     Omega_nodal DERIVATIVE NORMALIZING JUMP
! =========================================================         
! Derivate jump of perturbing function with respect to Omega_nodal     
!     ddomnod(i)=ddR(3)/(dpig**2) 
     ddomnod(i)=ddR(2)/(dpig**2) 

! ==================================================================
! 	NORMALIZING RESULT (average is integral over area)                
! ==================================================================
!    dg=ddg/(dpig**2) 
!    ddom=ddom/(dpig**2) 
!    ddzl=ddzl/(dpig**2)
!    ddomnod=ddomnod/(dpig**2)
  ENDDO

  RETURN 
END SUBROUTINE prsalto
    

SUBROUTINE ecc_anom(elm,effe,u)
  USE orbit_elements
  USE fund_const
  implicit none
  TYPE(orbit_elem),INTENT(IN) :: elm ! Keplerian elements a,e,I,Omega,omega
  REAL(KIND=dkind),INTENT(IN) :: effe ! true anomaly 
  REAL(KIND=dkind),INTENT(OUT) :: u ! eccentric anomaly
! end interface
  REAL(KIND=dkind) :: beta2,rsua,cosu,sinu

  IF(elm%coo.ne.'KEP') THEN
     WRITE(*,*)'prsalto: ERROR! the elements are not Keplerian!'
     STOP
  ENDIF
  
  beta2 = 1.d0-elm%coord(2)**2
  rsua = beta2/(1.d0 + elm%coord(2)*cos(effe))
  cosu = rsua*cos(effe)+ elm%coord(2)
  sinu = rsua*sin(effe)/sqrt(beta2)
  u = datan2(sinu,cosu)

END SUBROUTINE ecc_anom
