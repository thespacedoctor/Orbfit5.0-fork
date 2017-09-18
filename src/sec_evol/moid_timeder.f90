! ******************************************************************
! computes the time derivative of the secular evolution of the moid
! ******************************************************************
! written by G.F. Gronchi and C. Tardioli, 2011
! ------------------------------------------------------------------
! to switch to circular coplanar orbits for the planets
! set force_circ=.true. in MODULE planet_orbits.f90
PROGRAM moid_timeder
  USE critical_points,ONLY: nminx
  USE fund_const
  USE orbit_elements
  USE output_control
  USE planet_orbits!, ONLY: placar
  USE dyn_param
!  USE attributable
!  USE reference_systems, ONLY: observer_position,pvobs
!  USE station_coordinates, ONLY: statcode
  IMPLICIT NONE
!  INTEGER,PARAMETER :: dkind=KIND(1.d0)
  INTEGER,PARAMETER :: jmax=10000
  INTEGER :: j,npla(jmax)
  INTEGER :: iunmoid
  REAL(KIND=dkind) :: time,Hbar,zl 
  REAL(KIND=dkind) :: t_yr(jmax),omeg(jmax),Omnod(jmax),tcros,d_thres
  REAL(KIND=dkind) :: ecc(jmax),inc(jmax),numpl(jmax)
  TYPE(orbit_elem) :: elem ! asteroid elements
  TYPE(orb_uncert) :: unc_elem
  TYPE(orbit_elem) :: elea ! Earth elements
  REAL(KIND=dkind) :: deriv(nminx) 
  INTEGER :: fail_flag
  INTEGER :: nummin
  DOUBLE PRECISION,DIMENSION(nminx) :: dmintil,phatime  
!  all data on planetary semimajor axis, masses, etc.
  INCLUDE 'pldata.h90'
! ==============OPTIONS============================
  DOUBLE PRECISION :: elemt
  CHARACTER(LEN=60) eledir ! file names 
!  CHARACTER(LEN=60) file,elefil,eledir ! file names 
  CHARACTER(LEN=6) :: progna
  CHARACTER(LEN=80) :: run

  CHARACTER(LEN=9) :: name ! asteroid name 
  LOGICAL ::  eof, err
  CHARACTER*80 :: catnam
  INTEGER :: iunin
  INTEGER :: norb !number of orbits in the catalog
  INTEGER :: nlsloc,lsloc(ndyx)
  INTEGER, PARAMETER :: norbx=100000
  INTEGER :: h
! ===========================================================
! options
  progna='compmo' 
  run='moid_timeder'
  CALL compop

  d_thres=0.02d0
  write(*,*)'thresold=',d_thres
 
  CALL filopn(iunmoid,'moid_timeder.fla','unknown')

  write(*,*)'iunmoid=',iunmoid

! *** reading asteroid orbital data ***
  catnam='./epoch/neodys.ctc'
!  catnam='./catalog.ctm'
  CALL filopn(iunin,catnam,'OLD')
  write(*,*)'iunin=',iunin
  CALL oporbf(catnam,iunin)
  norb=0
  write(*,*)'name  time  MOID   t-deriv'
  DO j=1,norbx
     IF(j/100.gt.(j-1)/100) THEN
 !       WRITE(*,*) j
     ENDIF
! --- asteroid ---
!     CALL read_elems(elem,name,eof,nlsloc,err,lsloc,UNC=unc_elem,FILE='catalog.ctm',UNIT=iunin)
     CALL read_elems(elem,name,eof,nlsloc,err,lsloc,UNC=unc_elem,FILE='neodys.ctc',UNIT=iunin)
     CALL coo_cha(elem,'KEP',elem,fail_flag) 

     IF(.not.unc_elem%succ)THEN 
        WRITE(*,*)'Asteroid ',name,' covariance matrix not found.'
        WRITE(*,*)'skipping computation'
        WRITE(ierrou,*)'Asteroid ',name,' covariance matrix not found.'
        WRITE(ierrou,*)'skipping computation'
        !        STOP 
        CYCLE
     ENDIF
     
     IF(eof) EXIT
     norb=norb+1

! Earth's elements at the same time                                    
     elea=undefined_orbit_elem
     elea%t=elem%t
     IF(force_circ)THEN
        ! compute directly Keplerian elements
        CALL pla_dat(ap)
        elea%coord(1)=ap(3)
        elea%coord(2:6)=0.d0
        elea%coo='KEP'
     ELSE
!        CALL placar(13,elea%t,elea%coord,1) !Earth-Moon barycenter
        CALL placar(3,elea%t,elea%coord,1) !Earth only
        elea%coo='CAR'
        CALL coo_cha(elea,'KEP',elea,fail_flag) 
     ENDIF

     CALL pladat  ! HINT: alcune variabili erano gia' state assegnate

     tcros=0.d0       
     CALL dmintil_deriv(elea,elem,nummin,dmintil,deriv)     
     DO h=1,nummin
        IF(abs(dmintil(h)).le.d_thres.and.elem%h_mag.le.22)THEN
           IF(dmintil(h).lt.0.d0.and.deriv(h).lt.0.d0)THEN
              phatime(h)=(-d_thres-dmintil(h))/deriv(h)
           ELSEIF(dmintil(h).lt.0.d0.and.deriv(h).gt.0.d0)THEN
              phatime(h)=(d_thres-dmintil(h))/deriv(h)
              tcros = -dmintil(h)/deriv(h)
           ELSEIF(dmintil(h).gt.0.d0.and.deriv(h).lt.0.d0)THEN
              phatime(h)=(-d_thres-dmintil(h))/deriv(h)
              tcros = -dmintil(h)/deriv(h)
           ELSEIF(dmintil(h).gt.0.d0.and.deriv(h).gt.0.d0)THEN
              phatime(h)=(d_thres-dmintil(h))/deriv(h)
           ELSE
              WRITE(*,*)'LIMITING CASE: STOP'
              STOP 
           ENDIF
           elemt=elem%t/365.25+1858.87953
           write(iunmoid,101)name,elem%h_mag,elemt,dmintil(h),deriv(h),phatime(h),tcros
           write(*,*)'nummin=',nummin
           write(*,101) name,elem%h_mag,elemt,dmintil(h),deriv(h),phatime(h),tcros
        ENDIF
        
        !        write(iunmoid,*) t_yr(j)+1858.9+time/365.25,dmintil(h),deriv(h)
        !        write(*,*) t_yr(j)+1858.9+time/365.25,dmintil(h),deriv(h)
     ENDDO
  ENDDO
101 FORMAT(a9,2x,f6.2,2x,f10.3,2x,2(f15.8,2x),2(f18.3,2x))
  
  WRITE(*,*)' number of orbits in input ', norb
  CALL clorbf
  
  CALL filclo(iunin,' ')  
!  CALL filclo(iunevol,' ')
  CALL filclo(iunmoid,' ')

CONTAINS
! *****************************************************************  
  SUBROUTINE compop 
    IMPLICIT NONE
    INTEGER            :: le,iunout
    CHARACTER(LEN=100) :: file 
    LOGICAL            :: ireq,found
    CHARACTER(LEN=60)  :: comment 
! =============================  
    CALL initopt(progna,run,'mop')  
! read option for physical model and integration method                 
    CALL rmodel(1)    
    ! Output files: for control                 
    file=run//'.mou' 
    CALL rmsp(file,le) 
    call filopn(iunout,file(1:le),'UNKNOWN') 
! =============WHERE TO FIND ASTEROID ELEMENTS===================    
!    ireq=.true. 
!    comment='asteroid elements directory'
!    CALL input_cha_opt(progna,'eledir',eledir,ireq,found,comment,iunout)
! =============CHECK DERIVATIVES===================
!    ireq=.true.
!    comment='check derivatives'
!    CALL input_log_opt(progna,'chk_der',chk_der,ireq,found,comment,iunout)
  END SUBROUTINE compop

END PROGRAM moid_timeder


