! *********************************************************
! select, from the ASTDYS asteroid list, the ones that could have
! perihelion distance smaller than a given quantity
! and are brighter than a given magnitude
! (reading catalogs .ctm)  
! *********************************************************
! written by G.F. Gronchi & G.B. Valsecchi (October 2010)
! *********************************************************
PROGRAM possible_PNEAs
  USE orbit_elements
  USE fund_const
  USE output_control
  USE dyn_param
  IMPLICIT NONE
  INTEGER, PARAMETER :: norbx=1000000
! ===========for read_elems========================
  INTEGER :: h
  CHARACTER*80 :: catnam
  INTEGER :: iunin,iunpnea
  INTEGER :: norb !number of orbits in the catalog
! ==============OPTIONS============================
  CHARACTER(LEN=6) :: progna
  CHARACTER(LEN=80) :: run
! ===============ASTEROID/COMET====================
  CHARACTER(LEN=9) :: name ! asteroid name 
  LOGICAL ::  eof
  CHARACTER(LEN=60) file,eledir ! file names 
  INTEGER :: le
  INTEGER fail_flag
  TYPE(orbit_elem) :: eq,elkep
  TYPE(orb_uncert) :: unc
! =======================================================
  REAL(KIND=dkind) :: hmax
! for min perihelion distance
  REAL(KIND=dkind) :: max_pd
  REAL(KIND=dkind) :: ecc,emax,qmin,ci,si
  INTEGER :: count,sel
  INTEGER :: j  !loop indexes
  LOGICAL :: err
  INTEGER :: nlsloc
! =====================================================================

! options
  progna='compmo' 
!  run='lostPHAs'
  run='potential_NEAs'
!  CALL compop

  rhs=1
  count=0

  WRITE(*,*)'input max_pd:'
  READ(*,*) max_pd
  WRITE(*,*)'input hmax:'
  READ(*,*) hmax
13 WRITE(*,*)'numb or multiopp? (type 1 or 2):'
  READ(*,*) sel

! Open output files
  CALL rmsp(run,le)

  IF(sel.eq.1)THEN
     catnam='./orbits/allnum.ctm'
     CALL filopn(ierrou,run(1:le)//'.numb.err','unknown') !error file
     CALL filopn(iunpnea,'candidate_PNEAs.numb','unknown')
  ELSEIF(sel.eq.2)THEN
     catnam='./orbits/ufitobs.ctm'
     CALL filopn(ierrou,run(1:le)//'.mopp.err','unknown') !error file
     CALL filopn(iunpnea,'candidate_PNEAs.mopp','unknown')
  ELSE
     write(*,*)'invalid selection!'
     GOTO 13
  ENDIF

! *** reading orbital data ***
  CALL filopn(iunin,catnam,'OLD')
  CALL oporbf(catnam,iunin)
  norb=0
  DO h=1,norbx
     IF(h/10000.gt.(h-1)/10000) THEN
        WRITE(*,*) h
     ENDIF

     CALL read_elems(eq,name,eof,nlsloc,err,UNC=unc,FILE='allnum.ctm',UNIT=iunin)
     IF(eof) EXIT
     norb=norb+1

     CALL coo_cha(eq,'KEP',elkep,fail_flag)
     IF(fail_flag.ge.5)THEN
        WRITE(*,*)'potential_NEAs: error! failed coord change'
        STOP
     ENDIF

! compute max eccentricity and min perihelion distance
     IF(elkep%h_mag.le.hmax)THEN

        si=sin(elkep%coord(3))
        ci=cos(elkep%coord(3))
        ecc = elkep%coord(2)
        emax = sqrt(si**2 + (ecc*ci)**2)
        qmin = elkep%coord(1)*(1.d0-emax)

        IF(qmin.le.max_pd)THEN
!           WRITE(*,*) 'name,qmin,hmag:',name,qmin,elkep%h_mag
           WRITE(iunpnea,100) name,qmin,elkep%h_mag
           count=count+1
        ENDIF
     ENDIF

  ENDDO

100 FORMAT(a9,2x,f10.5,2x,f8.2)

  WRITE(*,*)' number of orbits in input ', norb
  WRITE(*,*)' number of orbits satisfying the criterion', count
  CALL clorbf

33 CONTINUE
  CALL filclo(iunpnea,' ')


END PROGRAM Possible_PNEAs
