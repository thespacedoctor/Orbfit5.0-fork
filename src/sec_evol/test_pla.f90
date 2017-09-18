PROGRAM test_pla
  USE fund_const
  USE orbit_elements
  USE planet_orbits
  IMPLICIT NONE

  CHARACTER*60 :: eledir
  CHARACTER*6 :: progna 
  CHARACTER*80 :: run,elefil
  INTEGER :: iunlst,iunpla,iunpla2,le
  INTEGER,PARAMETER :: maxint=50000 !YR
  INTEGER,PARAMETER :: step=100 !YR
  INTEGER,PARAMETER :: maxint2=50 !YR
  INTEGER,PARAMETER :: step2=30 !MJD
  TYPE(orbit_elem) :: elea,elea2
  INTEGER :: h,k,fail_flag,nevp
  DOUBLE PRECISION :: princ
  CHARACTER*3 :: ext !extension file planets secular elements
                     !'vpla.dat' or 'vpla.fil'
! **********************************************************************
  INCLUDE 'pldata.h90'

  CALL pladat
  progna='propne' 
  run='propneo'
  CALL optpro(progna,run,iunlst,eledir)

!  DO 1 k=1,2
!     IF(k.eq.1)THEN
        ext='fil'
!     ELSE
!        ext='dat'
!     ENDIF
  CALL read_pla(ext)
! Open output files
  CALL filnam('./planets','evol',ext,elefil,le) 
  CALL filopn(iunpla,elefil,'unknown')
  
  CALL filnam('./planets','evol_pla',ext,elefil,le)
  CALL filopn(iunpla2,elefil,'unknown')
 
  DO h=1,nstepx
     IF(tpla(h).le.0)THEN
        GOTO 3
     ENDIF
  ENDDO
3 CONTINUE
  nevp=h-1 !max tpla not zero

! Earth setting
  elea=undefined_orbit_elem
  elea%coord(6)=0.d0
  elea%coo='EQU'
  elea2=undefined_orbit_elem
  elea2%coord(6)=0.d0

  DO  h=0,maxint,step
     elea%t=teph0+(tpla(1)+h)*365.25d0
     IF(tpla(1)+h.gt.tpla(nevp))THEN
!           WRITE(*,*) 'end of planet elemens, epoch ',h
        GOTO 5
     ENDIF
     CALL planet_elems(3,elea%t,elea%coord(1:5))
     WRITE(iunpla,104) FLOAT(h),elea%coord(2:5) !time in yr 
104  FORMAT(f10.3,4(2x,f15.10))
  ENDDO
5 CONTINUE

  DO h=0,INT(maxint2*365.25d0),step2
     elea%t=teph0+tpla(1)*365.25d0+FLOAT(h)
     IF(tpla(1)+h/365.25d0.gt.tpla(nevp))THEN
 !          WRITE(*,*) 'end of planet elemens, epoch ',h
        GOTO 9
     ENDIF
     IF(ext.eq.'fil')THEN
        CALL planet_elems(3,elea%t,elea%coord(1:5))
        elea2%t=elea%t
        CALL placar(13,elea2%t,elea2%coord,1)
        elea2%coo='CAR'
        CALL coo_cha(elea2,'EQU',elea2,fail_flag)
        IF(fail_flag.ge.5)THEN
           write(*,*)'fail_flag=',fail_flag,'stopping program'
           STOP
        ENDIF !time in MJD
        WRITE(iunpla2,108) FLOAT(h)/365.25d0,elea%coord(2:5), &
             & elea2%coord(2:5)
     ELSE
        CALL planet_elems(3,elea%t,elea%coord(1:5))
        WRITE(iunpla2,108) FLOAT(h)/365.25d0,elea%coord(2:5)
     ENDIF
108  FORMAT(f10.3,8(2x,f15.10))
109  FORMAT(f10.3,4(2x,f15.10))
  ENDDO
9 CONTINUE

  CALL filclo(iunpla,' ')
  CALL filclo(iunpla2,' ')
!1 ENDDO

END PROGRAM test_pla


! ===================================================================   
! OPTPRO                                                                
! ===================================================================   
! input options, for the propagator and the specific main program       
! input: progna = program name (6 characters)                           
!        run    = run identifier (80 characters, up to 76 non-blank)    
  SUBROUTINE optpro(progna,run,iunlst,eledir) 
    character*6 progna 
    character*80 run
    integer iunlst 
    character*(*) eledir 
! ==========END INTERFACE============================================   
    CHARACTER*80 neolist
    integer le,iunit 
    character*12 file 
    logical found 
    LOGICAL ireq,fail1,fail 
    CHARACTER*7 prognp 
! =============================                                         
    CALL initopt(progna,run,'nop') 
! ==============================                                        

! read option for physical model and integration method                 
    CALL rmodel(1)    

! ============= OPEN ASTEROID NAME LIST===================              
    fail=.false. 
    ireq=.true. 
    prognp=progna//'.' 
    CALL rmsp(prognp,le) 
    CALL rdncha(prognp,'neolist',neolist,ireq,found,fail1,fail)
    IF(fail)THEN 
       WRITE(*,*)' optpro: list of asteroid names required' 
       STOP 
    ENDIF
    CALL filopn(iunlst,neolist,'OLD') 
! ============= WHERE TO FIND/PUT ASTEROID ELEMENTS===================  
    ireq=.true. 
    CALL rdncha(prognp,'eledir',eledir,ireq,found,fail1,fail)
    IF(fail)THEN 
       WRITE(*,*)' optpro: asteroid elements directory required' 
       STOP 
    ENDIF

  END SUBROUTINE optpro
