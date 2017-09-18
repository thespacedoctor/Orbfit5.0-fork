! ***************************************************************
!SUBROUTINE prdnod(iunref,y,a,dn,nnod)
SUBROUTINE prdnod(y,a,dn,nnod)
 USE fund_const
 USE orbit_elements
 USE planet_orbits
 USE critical_points
 IMPLICIT NONE
! =========INPUT================
! INTEGER,INTENT(IN) :: iunref !unit for mutual elements secular evolution 
! y(1)=omega, y(2)=G=ky*sqrt(a*(1-e**2)), y(3)=Omega, y(4)=zl=G*cos(I)
! a= semimaj axis (au)
 REAL(KIND=dkind),INTENT(IN) :: y(4),a
! =========================== OUTPUT ===============================
! dn= 2 x npl nodal distances
! nnod=planet number times (+1) if ascending, (-1) if descending nod
! indicates the nearest node
! -------------------- planet data  --------------------
! planet data
! INCLUDE 'parbep.h'
! INCLUDE 'masses.h'
 INCLUDE 'pldata.h90'
! local variables
 REAL(KIND=dkind),INTENT(OUT) :: dn(nplax2)
 INTEGER,INTENT(OUT) :: nnod
! ======================= END INTERFACE ============================
 REAL(KIND=dkind) :: e,dnp,dnm,tmp,adn
 INTEGER :: i
 REAL(KIND=dkind) :: apl,epl ! planet elements
 TYPE(orbit_elem) :: elpl ! planet orbital elements
 TYPE(orbit_elem) :: elem ! asteroid orbital elements
 REAL(KIND=dkind) :: mutI,mutom,mutompl ! mutual elements
 INTEGER :: fail_flag
! *****************************************************************
! Compute asteroid eccentricity
 e=sqrt(1.d0-(y(2)/ky)**2/a)

! Asteroid Keplerian Elements
 elem=undefined_orbit_elem
 elem%coo='KEP'
 elem%coord(1)=a
 elem%coord(2)=sqrt(1.d0-(y(2)/ky)**2/a)
 elem%coord(3)=acos(y(4)/y(2))
 elem%coord(4)=y(3)
 elem%coord(5)=y(1)
 elem%coord(6)=0.d0

! initialization
 adn=1.d10
 nnod=1000

! ============ loop on planets ==============
 DO i=inpl,ioupl
    elpl=el_pla(i)
    apl=elpl%coord(1) 
    epl=elpl%coord(2) 

! Compute mutual elements
    CALL mutualrefcha(elpl%coord,elem%coord,mutI,mutom,mutompl)

! WRITE(*,*) 'epl,Ipl',elpl%coord(2:3)
! WRITE(*,*) elem%coord(3:5)
! WRITE(*,*) mutI,mutompl,mutom
!    WRITE(iunref,100) elem%coord(1:2),elpl%coord(1:2),mutom,mutompl
!100 FORMAT(4(f10.7,1x),2(f11.7,1x))
    
! Compute nodal distances
    dnp=a*(1.d0-e**2)/(1.d0+e*cos(mutom))- &
         & apl*(1.d0-epl**2)/(1.d0+epl*cos(mutompl))
    dnm=a*(1.d0-e**2)/(1.d0-e*cos(mutom))- &
         & apl*(1.d0-epl**2)/(1.d0-epl*cos(mutompl))
    dn(2*i-1)=dnp
! write(*,*)'prdnod:dn(2i-1)',dn(2*i-1)

    IF(abs(dnp).lt.adn)THEN
       adn=abs(dnp)
       nnod=i
    ENDIF
    dn(2*i)=dnm
! write(*,*)'prdnod:dn(2i)',dn(2*i)

    IF(abs(dnm).lt.adn)THEN
       adn=abs(dnm)
       nnod=-i
    ENDIF

 ENDDO

END SUBROUTINE prdnod


! ******************************************************************
INTEGER FUNCTION nodcod(nnod)
 IMPLICIT NONE
 INTEGER :: nnod
! ***********************************

 IF(nnod.gt.0)THEN
    nodcod=2*nnod-1
 ELSEIF(nnod.lt.0)THEN
    nodcod=-2*nnod
 ELSE
    nodcod=0
 ENDIF

 RETURN
END FUNCTION nodcod


! ****************************************************************
INTEGER FUNCTION codnod(nc,iasc)
 IMPLICIT NONE
 INTEGER :: nc,iasc
! **********************************
 codnod=(nc+1)/2
 iasc=mod(nc,2)

 RETURN
END FUNCTION codnod


! ******************************************************************
SUBROUTINE prddnod(y,dy,a,dn,ddn,nnod)
  USE fund_const
  USE orbit_elements
  USE planet_orbits
!  USE critical_points
  IMPLICIT NONE
! =========================== INPUT ===============================
! y(1)=omega, y(2)=G=ky*sqrt(a*(1-e**2)), y(3)=Omega, y(4)=zl=G*cos(I)
! dy(i), i=1,4 time derivative of the same
! a= semimajor axis (au)
  REAL(KIND=dkind),INTENT(IN) :: y(4),dy(4),a
! =========================== OUTPUT ==============================
! dn= 2 x npl nodal distances
! ddn= time derivative of the same
! nnod=planet number times (+1) if ascending, (-1) if descending nod
! indicates the nearest node
! --------------- planet data ---------------
! planet data
! INCLUDE 'parbep.h'
! INCLUDE 'masses.h' 
  INCLUDE 'pldata.h90'
  REAL(KIND=dkind),INTENT(OUT) :: dn(2*nplax),ddn(2*nplax)
  INTEGER :: nnod
! local variables
  REAL(KIND=dkind) :: adn
! ====================== END INTERFACE  ============================== 
  TYPE(orbit_elem) :: elpl ! planet orbital elements
  TYPE(orbit_elem) :: elem ! asteroid orbital elements
!  REAL(KIND=dkind) :: mutI,mutom,mutompl ! mutual elements
  REAL(KIND=dkind) :: som,com,si,ci,son,con
  REAL(KIND=dkind) :: somp,comp,sip,cip,sonp,conp
  REAL(KIND=dkind) :: sinmutom,sinmutompl,cosmutom,cosmutompl,cosmutI
  REAL(KIND=dkind),DIMENSION(3) :: N_pl,N_ast,xa,xpl,A_nod,Axa,Axpl
  REAL(KIND=dkind) :: modA,vsize
  REAL(KIND=dkind),DIMENSION(3) :: dNdI,dNdon,dAndI,dAndon,dAdI,dAdon
  REAL(KIND=dkind),DIMENSION(3) :: dxdI,dxdom,dxdon
  REAL(KIND=dkind) :: dmutomdI,dmutomdom,dmutomdon
  REAL(KIND=dkind) :: dcosmutomdI,dcosmutomdom,dcosmutomdon
  REAL(KIND=dkind) :: dmutompldI,dmutompldon
  REAL(KIND=dkind) :: dcosmutompldI,dcosmutompldon
  REAL(KIND=dkind) :: dmutom,dmutompl
  REAL(KIND=dkind) :: apl,epl,e,beta
  REAL(KIND=dkind) :: upeco,upecopl,umeco,umecopl,dnp,dnm,ddnp,ddnm
  REAL(KIND=dkind) :: ddnpde,ddnmde,de,ddnpdmutom,ddnmdmutom,dIdt,depl
  REAL(KIND=dkind) :: ddnpdmutompl,ddnmdmutompl,ddnpdepl,ddnmdepl
  INTEGER :: i,j,k,fail_flag
! Output of mutual_fef, subroutine in critical_points
!  DOUBLE PRECISION,DIMENSION(2) :: dsindcos
!  DOUBLE PRECISION :: dsinmutIdcos
!  DOUBLE PRECISION,DIMENSION(3,2) :: dcosmutI
! Matrix of rotation from inertial to mutual reference frame
!  DOUBLE PRECISION,DIMENSION(3,3) :: rotinmut
! Derivatives of the matrix of rotation from inertial to 
! mutual reference frame
!  DOUBLE PRECISION,DIMENSION(10,3,3) :: drotinmut
! Derivatives of cosom1 and cosom2 (to compute the derivatives
! of the nodal distances)
!  DOUBLE PRECISION,DIMENSION(2,3,2) :: dcosom

! ************************************************************************
! --- Declaretion of asteroid elements ------------------------
  elem%coo='KEP'
  elem%coord(1)=a
  elem%coord(2)=sqrt(1.d0-(y(2)/ky)**2/a)
  elem%coord(3)=acos(y(4)/y(2))
  elem%coord(4)=y(3)
  elem%coord(5)=y(1)
  elem%coord(6)=0.d0

!  a=elem%coord(1) ! asteroid semiaxis
  e=elem%coord(2) ! asteroid eccentricity [e=sqrt(1.d0-(y(2)/ky)**2/a)]
  beta=sqrt(1.d0-e**2) 

! ------------------ Declaretion of variables ----------------------------
  som=sin(elem%coord(5))
  com=cos(elem%coord(5))
  si=sin(elem%coord(3))
  ci=cos(elem%coord(3))
  son=sin(elem%coord(4))
  con=cos(elem%coord(4))
 
! ------------------- Declaretion of vectors  -----------------------------
! Versor pointing to the angular momentum of the asteroid
  N_ast(1)=son*si
  N_ast(2)=-con*si
  N_ast(3)=ci

! Derivative of N_ast w.r.t. I
  dNdI(1)=son*ci
  dNdI(2)=-con*ci
  dNdI(3)=-si

! Derivative of N_ast w.r.t. Omnod
  dNdon(1)=con*si
  dNdon(2)=son*si
  dNdon(3)=0.d0

! Unit vector X(0) pointing to the position of the perihelion 
! of the orbit of the asteroid 
  xa(1)=con*com-son*som*ci
  xa(2)=son*com+con*som*ci
  xa(3)=som*si

! Derivative of X(0) w.r.t. omega
  dxdI(1)=son*som*si
  dxdI(2)=-con*som*si
  dxdI(3)=som*ci

! Derivative of X(0) w.r.t. omega
  dxdom(1)=-con*som-son*com*ci
  dxdom(2)=-son*som+con*com*ci
  dxdom(3)=com*si

! Derivative of X(0) w.r.t. omnod
  dxdon(1)=-son*com-con*som*ci
  dxdon(2)=con*com-son*som*ci
  dxdon(3)=0.d0

! initialization
  adn=1.d10
  nnod=1000

! =========== loop on planets =============
  DO i=inpl,ioupl

     elpl=el_pla(i)

! Compute mutual elements
!     CALL mutual_ref(elpl,elem,mutI,mutompl,mutom,dsindcos, &
!     & dsinmutIdcos,dcosmutI,rotinmut,drotinmut,dcosom)
!     write(*,*) mutI,mutompl,mutom
!     CALL mutualrefcha(elpl%coord,elem%coord,mutI,mutom,mutompl)
!     write(*,*) mutI,mutompl,mutom

! ------------------ Declaretion of variables ----------------------------
     somp=sin(elpl%coord(5))
     comp=cos(elpl%coord(5))
     sip=sin(elpl%coord(3))
     cip=cos(elpl%coord(3))
     sonp=sin(elpl%coord(4))
     conp=cos(elpl%coord(4))

! ------------------- Declaretion of vectors  -----------------------------
! Versor pointing to the angular momentum of the planet
     N_pl(1)=sonp*sip
     N_pl(2)=-conp*sip
     N_pl(3)=cip

! cos(mutI) = DOT_PRODUCT(N_pl,N_ast) 
!     write(*,*) DOT_PRODUCT(N_pl,N_ast)
!     write(*,*) cos(mutI)

! Unit vector X'(0) pointing to the position of the perihelion
! of the orbit of the planet
     xpl(1)=conp*comp-sonp*somp*cip
     xpl(2)=sonp*comp+conp*somp*cip
     xpl(3)=somp*sip

! ==================================================================
! COMPUTE ASCENDING MUTUAL NODE VECTOR                   
     cosmutI=DOT_PRODUCT(N_pl,N_ast) 
! Check that N_pl and N_ast are not parallel                           
     IF ((abs(abs(cosmutI)-1.d0)).lt.1.d-10) THEN 
        WRITE(*,*)'MUTUAL INCLINATION CLOSE TO ZERO!' 
!        read(*,*)
!        mutI=0.d0 
!        mutom=y(1)+y(3) 
!        mutompl=elpl%coord(5)+elpl%coord(4)
! setting derivatives
        cosmutom=com*con-som*son
        sinmutom=som*con+com*son
        dmutomdI=0.d0
        dmutomdom=1.d0
        dmutomdon=1.d0
        cosmutompl=comp*conp-somp*sonp
        sinmutompl=somp*conp+comp*sonp
        dmutompldI=0.d0
        dmutompldon=0.d0
! skip next steps                                                   
        GOTO 18 
     ELSE
! Vector A_nod = N_pl x N_ast pointing to the ascending mutual node
        CALL prvec(N_pl,N_ast,A_nod)
        modA=vsize(A_nod) ! Absulute value of A_nod
     ENDIF
! ==================================================================
 
     cosmutom=DOT_PRODUCT(A_nod(1:3)/modA,xa)
     CALL prvec(A_nod(1:3)/modA,xa,Axa)
     sinmutom=DOT_PRODUCT(Axa,N_ast)
!     write(*,*) sinmutom, sin(mutom)
!     write(*,*) cosmutom,cos(mutom)
! cos(mutom) = DOT_PRODUCT(A_nod/modA,X(0))
!     write(*,*) DOT_PRODUCT(A_nod(1:3)/modA,xa)
!     write(*,*) cos(mutom)

     cosmutompl=DOT_PRODUCT(A_nod(1:3)/modA,xpl)
     CALL prvec(A_nod(1:3)/modA,xpl,Axpl)
     sinmutompl=DOT_PRODUCT(Axpl,N_pl)
!     write(*,*) sinmutompl,sin(mutompl)
!     write(*,*) cosmutompl,cos(mutompl)
! cos(mutompl) = DOT_PRODUCT(A_nod/modA,X'(0))
!     write(*,*) DOT_PRODUCT(A_nod,xpl)/modA
!     write(*,*) cos(mutompl)
!     read(*,*)

! ------------------ Derivative w.r.t. I -------------------------
! Actual is -sin(ommut)*domtdI=<dA_noddI,X(0)>+<A_nod,dX(0)dI>

! Derivative of A_nod w.r.t. I
     CALL prvec(N_pl,dNdI,dAndI)
     
! Derivative oh the versor A_nod/modA w.r.t. I
     dAdI=0.d0
     DO j=1,3
        dAdI(j)=(1.d0/modA)*dAndI(j)- &
             & A_nod(j)*DOT_PRODUCT(A_nod,dAndI)/modA**3
     ENDDO

! Derivative of omega_mutual w.r.t. I
     dcosmutomdI=DOT_PRODUCT(dAdI,xa)+DOT_PRODUCT(A_nod,dxdI)/modA
     dmutomdI=-dcosmutomdI/sinmutom
!     write(*,*)'======================'
!     write(*,*) dcosmutomdI,dcosom(2,1,2)
!     read(*,*)
! Derivative of omega_mutual_prime w.r.t. I
     dcosmutompldI=DOT_PRODUCT(dAdI,xpl)
     dmutompldI=-dcosmutomdI/sinmutompl
!     write(*,*)'======================'
!     write(*,*) dcosmutompldI,dcosom(1,1,2)
!     read(*,*)

! ------------------- Derivative w.r.t. omega -------------------------
! Actual is -sin(ommut)*dommutdI=<A_nod,dX(0)dom>

! Derivative of omega_mutual w.r.t. omega
     dcosmutomdom=DOT_PRODUCT(A_nod,dxdom)/modA
     dmutomdom=-dcosmutomdom/sinmutom
!     write(*,*)'======================'
!     write(*,*) dcosmutomdom,dcosom(2,2,2),dcosom(1,2,2)
!     read(*,*)
! Derivative of omega_mutual_prime w.r.t. omega
! dmutompldom=0.d0

! ------------------- Derivative w.r.t. Omnod --------------------------
! Actual is -sin(ommut)*dommutdI=<dA_noddI,X(0)>+<A_nod,dX(0)dI>

! Derivative of A_nod w.r.t. Omnod
     CALL prvec(N_pl,dNdon,dAndon)

! Derivative oh the versor A_nod/modA w.r.t. Omnod
     dAdon=0.d0
     DO j=1,3
        dAdon(j)=(1.d0/modA)*dAndon(j)- &
             & A_nod(j)*DOT_PRODUCT(A_nod,dAndon)/modA**3
     ENDDO

! Derivative of omega_mutual w.r.t. Omnod
     dcosmutomdon=DOT_PRODUCT(dAdon,xa)+DOT_PRODUCT(A_nod,dxdon)/modA
     dmutomdon=-dcosmutomdon/sinmutom
!     write(*,*)'======================'
!     write(*,*) dcosmutomdon,dcosom(2,3,2)
!     read(*,*)
! Derivative of omega_mutual_prime w.r.t. Omnod
     dcosmutompldon=DOT_PRODUCT(dAdon,xpl)
     dmutompldon=-dcosmutompldon/sinmutompl
!     write(*,*)'======================'
!     write(*,*) 'don', dcosmutompldon,dcosom(1,3,2)
!     read(*,*)

18   CONTINUE
! ========================================================================
! --- Declaretion of planet elements ------------------------
     apl=elpl%coord(1) ! planet semiaxis
     epl=elpl%coord(2) ! planet eccentricity

! ---------- Node distances ---------------------------------------------
     upeco=1.d0+e*cosmutom 
     upecopl=1.d0+epl*cosmutompl
! ascending node distance
     dnp=a*(1.d0-e**2)/upeco-apl*(1.d0-epl**2)/upecopl 

     umeco=1.d0-e*cosmutom
     umecopl=1.d0-epl*cosmutompl
! descending node distance  
     dnm=a*(1.d0-e**2)/umeco-apl*(1.d0-epl**2)/umecopl 

! ---------------- Derivatives w.r.t. e ----------------------------------
! Derivative of the asteroid eccentricity w.r.t. time: dedt=dedG*dGdt
     de=-dy(2)*y(2)/sqrt(ky**2*a*(ky**2*a-y(2)**2)) 

! Ascending node distance derivative w.r.t. e
     ddnpde=-(a/upeco**2)*(2.d0*e*upeco+(1.d0-e**2)*cosmutom)
     ddnpde=ddnpde*de

! Descending node distance derivative w.r.t. e
     ddnmde=-(a/umeco**2)*(2.d0*e*umeco-(1.d0-e**2)*cosmutom)
     ddnmde=ddnmde*de

! ---------------- Derivatives w.r.t. epl  ----------------------------------
! Derivative of the planet eccentricity w.r.t. time: depldt
     depl=0.d0

! Ascending node distance derivative w.r.t. epl
     ddnpdepl=(apl/upeco**2)*(2.d0*epl*upeco-(1.d0-epl**2)*cosmutom)
     ddnpdepl=ddnpdepl*depl

! Descending node distance derivative w.r.t. epl
     ddnmdepl=(apl/umeco**2)*(2.d0*epl*umeco+(1.d0-epl**2)*cosmutom)
     ddnmdepl=ddnmdepl*depl

! ---------------- Derivatives w.r.t. mutom ------------------------------
! Derivative of mutom w.r.t. time, mutom=mutom(I,om,omnod)
     dIdt=1.d0/(ky*sqrt(a)*si*beta)*(dy(2)*ci-dy(4))
     dmutom=dmutomdI*dIdt+dmutomdom*dy(1)+dmutomdon*dy(3)

! Ascending node distance derivative w.r.t. mutom
     ddnpdmutom=a*(1.d0-e**2)*e*sinmutom/upeco**2
     ddnpdmutom=ddnpdmutom*dmutom

! Descending node distance derivative w.r.t. mutom
     ddnmdmutom=-a*(1.d0-e**2)*e*sinmutom/umeco**2
     ddnmdmutom=ddnmdmutom*dmutom

! ---------------- Derivatives w.r.t. mutompl ----------------------------
! Derivative of mutompl w.r.t. time, mutompl=mutompl(I,om,omnod)
     dmutompl=dmutompldI*dIdt+dmutompldon*dy(3)

! Ascending node distance derivative w.r.t. mutompl
     ddnpdmutompl=-apl*(1.d0-epl**2)*epl*sinmutompl/upecopl**2
     ddnpdmutompl=ddnpdmutompl*dmutompl

! Descending node distance derivative w.r.t. mutompl
     ddnmdmutompl=apl*(1.d0-epl**2)*epl*sinmutompl/umecopl**2
     ddnmdmutompl=ddnmdmutompl*dmutompl

! ---------------- Derivative Node Distance ----------------------------
     ddnp=ddnpde+ddnpdmutom+ddnpdmutompl+ddnpdepl
     ddnm=ddnmde+ddnmdmutom+ddnmdmutompl+ddnmdepl
! ************************************************************************

     dn(2*i-1)=dnp
     ddn(2*i-1)=ddnp
     IF(abs(dnp).lt.adn)THEN
        adn=abs(dnp)
        nnod=i
     ENDIF
     
     dn(2*i)=dnm
     ddn(2*i)=ddnm
      IF(abs(dnm).lt.adn)THEN
        adn=abs(dnm)
        nnod=-i
     ENDIF

  ENDDO
  
  RETURN
END SUBROUTINE prddnod
