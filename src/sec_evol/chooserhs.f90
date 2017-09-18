! ***********************************************************       
! Subroutine selecting method for computing right hand side         
! of Hamilton's equations                                           
! ***********************************************************       
SUBROUTINE chooserhs(y,om,omnod,g,zl,aa,ddd,eee,nnn)
  USE fund_const 
  USE orbit_elements
  USE planet_orbits
  USE critical_points
  USE right_hand_side
  USE kantorovich, ONLY : split!,elem,elpl
  IMPLICIT NONE
! ================= INTERFACE ================================ 
  REAL(KIND=dkind),INTENT(IN) :: y(4) !y=[omega,G,Omnod,Z] ! values at previous step
  REAL(KIND=dkind),INTENT(IN) :: om,omnod,g,zl !values at intermediate time steps
  REAL(KIND=dkind),INTENT(IN) :: aa
  REAL(KIND=dkind),DIMENSION(4),INTENT(OUT) :: ddd,eee
  INTEGER,DIMENSION(4),INTENT(OUT) :: nnn
! ================= END INTERFACE ===========================
  TYPE(orbit_elem) :: elem,elpl !planet and asteroid elems
  REAL(KIND=dkind),DIMENSION(4) :: dd,dd1,ee,ee1,dd_tmp
  INTEGER,DIMENSION(4) :: nn,nn1
  REAL(KIND=dkind) :: gmp !GM of planets
  REAL(KIND=dkind) :: ea,apl,epl
  REAL(KIND=dkind) :: eta !decision parameter 
  REAL(KIND=dkind) :: mutI,mutom,mutompl !mutual elems
  INTEGER :: i,j ! loop indexes
  INCLUDE 'pldata.h90'
! nodal distances
  REAL(KIND=dkind) :: dnodp(nplax),dnodm(nplax)
! dmintil
  INTEGER :: nummin  !number of minimal points
  REAL(KIND=dkind),DIMENSION(nminx) :: dmintil !minimal distances
  REAL(KIND=dkind) :: dmin_1,dmin_2 !fst,snd dmintil
! =================================================================
! for secpert
  REAL(KIND=dkind),DIMENSION(6) :: ast_elem
  REAL(KIND=dkind) :: hamil,hamil1,error ! output
  INTEGER :: numb
! for derrhs : test on 2nd-derivatives of R
  REAL(KIND=dkind),DIMENSION(4) :: ddd_ii,eee_ii,derR
  INTEGER,DIMENSION(4) :: nnn_ii
! =================================================================
! *************************************************************     
  ea=sqrt(1.d0-(y(2)/ky)**2/aa) !asteroid eccentricity

! Asteroid Keplerian Elements
  elem = undefined_orbit_elem
  elem%coo = 'KEP'
  elem%coord(1) = aa
  elem%coord(2) = ea
  elem%coord(3) = acos(y(4)/y(2))
  elem%coord(4) = y(3)
  elem%coord(5) = y(1)
  elem%coord(6) = 0.d0

! initialization                                                    
  DO j = 1,4 
     ddd(j) = 0.d0 
     eee(j) = 0.d0 
     nnn(j) = 0
  ENDDO

! loop on number of planets       
  DO 10 i = inpl,ioupl 

     elpl = el_pla(i)

     apl = elpl%coord(1) !planet semiaxis
     epl = elpl%coord(2) !planet ecc
     gmp=gm(i) 

     eta = 1.d0*(sqrt(gmp/39.4769))
!     write(*,*)'eta : ',eta 

! ----------------- NODAL DISTANCE COMPUTATION ------------------
! Compute mutual elements
     CALL mutualrefcha(elpl%coord,elem%coord,mutI,mutom,mutompl)
!     write(*,*)'mutual elements : ',mutI,mutom,mutompl

! Compute absulute value of nodal distances
     dnodp(i)=abs(aa*(1.d0-ea**2)/(1.d0+ea*cos(mutom)) -  &
          &  apl*(1.d0-epl**2)/(1.d0+epl*cos(mutompl))) 
     dnodm(i)=abs(aa*(1.d0-ea**2)/(1.d0-ea*cos(mutom)) -  &
          &  apl*(1.d0-epl**2)/(1.d0-epl*cos(mutompl)))
!     write(*,*)'asc and desc nod distance : ',dnodp(i),dnodm(i),i
! ---------------------------------------------------------------

! ------------------- DMINTIL_RMS CALL---------------------------
     elem%t = elpl%t
     CALL dmintil_rms(elpl,elem,nummin,dmintil)
     dmin_1 = ABS(dmintil(1))
     IF(nummin.eq.1)THEN
        dmin_2 = 1.d10
     ELSEIF(nummin.ge.2)THEN
        dmin_2 = ABS(dmintil(2))
     ENDIF
!    write(*,*)'dmintil values: ',dnodp(i),dnodm(i),i

     IF (dmin_1.gt.eta.AND.dmin_2.gt.eta) THEN 
        CALL rhs2(i,elpl,om,omnod,g,zl,aa,dd,ee,nn)
     ELSE
        CALL split(i,elpl,om,omnod,g,zl,aa,dd,ee,nn)
     ENDIF
     
! check: compare RHS and split outputs
!        write(13,*)'time,planet,nummin',elem%t/365.25+1858.87953d0,i,nummin
!        CALL check_split(dd,dd_tmp)

     DO j = 1,4 
        ddd(j) = ddd(j) + dd(j) 
        eee(j) = eee(j) + ee(j) 
        nnn(j) = nnn(j) + nn(j)
     ENDDO

10 END DO

  RETURN 
END SUBROUTINE chooserhs

SUBROUTINE check_split(dd,dd1)
  USE fund_const
  IMPLICIT NONE
! ------------------------- interface ------------------------------
  REAL(KIND=dkind),DIMENSION(4),INTENT(IN) :: dd,dd1
! ----------------------- end interface -----------------------------
  INTEGER :: k

  DO k=1,4
     write(13,*)dd(k),dd1(k)
     IF (ABS(dd(k)-dd1(k)).gt.1.d-7)THEN
        write(*,*)'chooserhs: RHS2 and SING EXTRACTION too different', &
             & dd1(1),dd(1)
     ENDIF
  ENDDO
 
END SUBROUTINE check_split
