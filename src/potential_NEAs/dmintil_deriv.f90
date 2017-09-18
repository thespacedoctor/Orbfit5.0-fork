! --------------------------------------------------------------------
! SUBROUTINE DMINTIL_DERIV
! written by G.F. Gronchi and C. Tardioli, July 2011
! --------------------------------------------------------------------
! compute minimal distance \tilde{d}_h and their time derivatives
! betwee the orbit of a planet (elea) and the orbit
! of an asteroid (elast).
! the elements of the planets are 
! 1) interpolated from vpla.fil if force_circ=.false., 
! 2) read from pla_dat if force_circ=.true.
! WARNING: note that the time derivatives are computed as Lie derivative of 
! \tilde{d}_h along the averaged vector field, and we assume nonzero
! only the components relative to G,Z,g,z of the asteroid; in particular
! the computation with the planets on non-circular orbits are not exact.

!MODULE derdmin
!  USE fund_const
!  USE orbit_elements
!  IMPLICIT NONE
!  PRIVATE
  
!CONTAINS

  SUBROUTINE dmintil_deriv(ele,elast,nummin,dmintil,deriv)
    USE critical_points
    USE orbit_elements
    USE fund_const
!  USE planet_orbits!, ONLY: placar
    USE right_hand_side
  IMPLICIT NONE
  TYPE(orbit_elem),INTENT(IN) :: elast !asteroid elements
  TYPE(orbit_elem),INTENT(IN) :: ele !Earth elements
  INTEGER,INTENT(OUT) :: nummin
  REAL(KIND=dkind),INTENT(OUT) :: dmintil(nminx) 
  REAL(KIND=dkind),INTENT(OUT) :: deriv(nminx) 
! end interface -----------------------------------------
  TYPE(orbit_elem) :: elem ! Keplerian asteroid elements
  TYPE(orbit_elem) :: elea ! Keplerian Earth elements
  TYPE(orbit_elem),DIMENSION(8) :: elpl ! Keplerian elements of a planet
  REAL(KIND=dkind) :: grad_dmintil(4,nminx) 
  REAL(KIND=dkind) :: averaged_rhs(4) 
  REAL(KIND=dkind),DIMENSION(5,nminx) :: ddmintdel2 !derivs w.r.t COM elems
  REAL(KIND=dkind) :: beta,sinI,cosI
  REAL(KIND=dkind),DIMENSION(5,4):: dKEP_dDEL 
  INTEGER :: fail_flag
  TYPE(orbit_elem) :: com2 ! cometary asteroid elements
  REAL(KIND=dkind),DIMENSION(6,6) :: dCOM_dKEP
  REAL(KIND=dkind),DIMENSION(1,5) :: grad_tmp,vectmp
  INTEGER :: j
  TYPE(orb_uncert) :: unc1,unc2
! for rhs2
  REAL(KIND=dkind) :: G,zl
  !output of rhs2: right hand side,errors,number eval                        
  REAL(KIND=dkind) :: ddd(4),eee(4) 
  INTEGER :: nnn(4) 

  elem=elast
  elea=ele
  IF(elem%coo.ne.'KEP')THEN
!     write(*,*)'converting to KEP'
     CALL coo_cha(elem,'KEP',elem,fail_flag)
     IF(fail_flag.ge.5)THEN
        write(*,*)'fail_flag=',fail_flag,'stopping program'
        STOP
     ENDIF
  ENDIF
  IF(elea%coo.ne.'KEP')THEN
!     write(*,*)'converting to KEP'
     CALL coo_cha(elea,'KEP',elea,fail_flag)
     IF(fail_flag.ge.5)THEN
        write(*,*)'fail_flag=',fail_flag,'stopping program'
        STOP
     ENDIF
  ENDIF

  beta = sqrt(1.d0-elem%coord(2)**2)
  sinI=sin(elem%coord(3))
  cosI=cos(elem%coord(3))

! compute Delaunay's elements
  G=ky*sqrt(elem%coord(1))*beta
  zl=G*cosI

  DO j=1,4
     averaged_rhs(j)=0.d0
  ENDDO

! Call rhs2 to compute the derivative of the right hand side 

  DO j=inpl,ioupl
     elpl(j)=undefined_orbit_elem
     elpl(j)%t=elea%t
!     IF(force_circ)THEN
        ! compute directly Keplerian elements
        CALL pla_dat(ap)
        elpl(j)%coord(1)=ap(j)
        elpl(j)%coord(2:6)=0.d0
        elpl(j)%coo='KEP'
!     ELSE
!        IF(j.eq.3)THEN
!           CALL placar(13,elea%t,elpl(j)%coord,1)
!        ELSE
!           CALL placar(j,elea%t,elpl(j)%coord,1)
!        ENDIF
!        elpl(j)%coo='CAR'
!        CALL planet_elems(j,elpl(j)%t,elpl(j)%coord(1:5))
!        elpl(j)%coo='EQU'
!        CALL coo_cha(elpl(j),'KEP',elpl(j),fail_flag)
!        IF(fail_flag.lt.5)THEN
!                   write(*,*)'j,elpl KEP:',j,elpl(j)%coord
!        ELSE
!           write(*,*)'fail_flag=',fail_flag,'stopping program'
!           STOP
!        ENDIF
!     ENDIF


! to use singularity extraction
!  CALL choose_rhs(elem%coord(5),G,zl,elem%coord(1),ddd,eee,nnn)
! ------------------------------------------------------
     CALL rhs2(j,elpl(j),elem%coord(5),elem%coord(4),G,zl, &
          & elem%coord(1),ddd,eee,nnn)
!  write(*,*)'ddd',ddd
     averaged_rhs(1) = averaged_rhs(1) - ddd(1) !-dR/dG
     averaged_rhs(2) = averaged_rhs(2) - ddd(3) !-dR/dZ
     averaged_rhs(3) = averaged_rhs(3) + ddd(2) ! dR/dg
     averaged_rhs(4) = averaged_rhs(4) + ddd(4) ! dR/dz

 ENDDO

! Jacobian matrix dKEP_dDEL
  dKEP_dDEL=0.d0
  dKEP_dDEL(5,1)=1.d0
  dKEP_dDEL(4,2)=1.d0
  dKEP_dDEL(2,3)=-beta/(ky*sqrt(elem%coord(1))*elem%coord(2))
  dKEP_dDEL(3,3)=cosI/(sinI*ky*sqrt(elem%coord(1))*beta)
  dKEP_dDEL(3,4)=-1.d0/(ky*sqrt(elem%coord(1))*beta*sinI)

! Jacobian matrix dCOM_dKEP
  CALL coo_cha(elem,'COM',com2,fail_flag,dCOM_dKEP)

! computing dmintil and their derivatives w.r.t. COM   
  CALL dmintil_rms(elea,elem,nummin,dmintil,DDMINTDEL2=ddmintdel2)
!  write(*,*)'nummin,dmintil(1:nummin)',nummin,dmintil(1:nummin)
!  write(*,*)'ddmintdel2',ddmintdel2(1:5,1)

  DO j=1,nummin
     vectmp(1,1:5)=ddmintdel2(1:5,j)
     grad_tmp(1,1:5) = MATMUL(vectmp(1,1:5),dCOM_dKEP(1:5,1:5))
     grad_dmintil(1:4,j) = MATMUL(grad_tmp(1,1:5),dKEP_dDEL(1:5,1:4))
     deriv(j) = DOT_PRODUCT(grad_dmintil(1:4,j),averaged_rhs(1:4))
  ENDDO

END SUBROUTINE dmintil_deriv

SUBROUTINE pla_dat(ap)
!  USE planet_masses
  IMPLICIT NONE
  DOUBLE PRECISION ap(9)
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
  RETURN 
END SUBROUTINE pla_dat

