!
! Main program for orbit determination and ephemeris generation
!
PROGRAM ephem
  USE file_oper
  USE output_control
  USE astrometric_observations
  USE orbit_elements
  USE fund_const
  USE ofod
  USE rdopt_of
  USE dyn_param
  IMPLICIT NONE


! DRYX: FOUND IN 'parobx.h90'
! MAX NUMBER OF OBSERVATIONS
  INTEGER nobx,nob2x
  PARAMETER (nobx=10000,nob2x=nobx*2)
! INTERFACE WITH blocks etc.
! block variables: max no blocks, max no residuals in block   
! block allocation:
  INTEGER, PARAMETER :: nblx= 1000 ! max number of passages for one asteroid
  DOUBLE PRECISION, PARAMETER ::  dtblock= 30.D0 ! min gap for block termination
  INTEGER, PARAMETER :: nxinbl=80 ! 2*max number of observations in one passage
!INTEGER, PARAMETER :: nxinbl_large=800, nbl_large=50 ! 2*max no obs, number of exceptional cases
  INTEGER, PARAMETER :: nxinbl_large=400, nbl_large=50 ! 2*max no obs, number of exceptional cases
! --------------------------  


! DRYX: FOUND IN parnob.h90
  INTEGER nobjx,nobj1x
  PARAMETER (nobjx=2)
  PARAMETER (nobj1x=nobjx+1)
! --------------------------


! DRYX: FOUND IN 'comlsf.h90'
! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: December 7, 1998
! ---------------------------------------------------------------------
! Options for least squares orbital fit
!
! delcr       -  Convergency control (correction norm)
! iiclsf      -  Initialization check
!
  INTEGER iiclsf
!      DOUBLE PRECISION delcr
  COMMON/cmlsf1/iiclsf
!      COMMON/cmlsf2/delcr
! -------------------------

! DRYX: FOUND IN 'comeph.h90'
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 28, 1997
! ---------------------------------------------------------------------
! Options for ephemeris generation (ORBFIT)
!
!
! iiceph      -  Initialization check
! kepobj      -  List of objects for which ephemeris is done
! nepobj      -  Number of entries in list kepobj
! teph1       -  Starting time for ephemeris (MJD, TDT)
! teph2       -  Ending time for ephemeris (MJD, TDT)
! dteph       -  Ephemeris stepsize (d)
! idsta       -  Observatory code
! ephtsc      -  Time scale
! ephfld      -  Output fields
!
      INTEGER idsta,kepobj(3),nepobj,iiceph
      DOUBLE PRECISION teph1,teph2,dteph
      CHARACTER*10 ephtsc
      CHARACTER*100 ephfld
      COMMON/cmeph1/idsta,kepobj,nepobj,iiceph
      COMMON/cmeph2/teph1,teph2,dteph
      COMMON/cmeph3/ephtsc
      COMMON/cmeph4/ephfld
! -------------------------

! =========================== DECLARATIONS ===========================
! Fixed FORTRAN I/O units
  INTEGER                            :: unirep,uniele,unidif,unieph
! File names
  CHARACTER(LEN=80)                  :: run,file,optfil,eleout

! Object names and input directories
  CHARACTER(LEN=80)                  :: name(nobj1x),nameo(nobj1x),dir(nobj1x)
  CHARACTER(LEN=80)                  :: namof(nobj1x)

! OBSERVATIONS
  INTEGER                            :: nt,n(nobjx)
! new data types
  TYPE(ast_obs),DIMENSION(nobx)      :: obs
  TYPE(ast_wbsr),DIMENSION(nobx)     :: obsw
  CHARACTER(LEN=20)                  :: error_model ! weighing model
  INTEGER                            :: ip1(nobjx),ip2(nobjx)
! Max number of orbital element files (for each object)
  INTEGER                            :: nifx
  PARAMETER (nifx=10)

! ORBITS
  TYPE(orbit_elem),DIMENSION(nobj1x) :: elem
  TYPE(orb_uncert),DIMENSION(nobj1x) :: elem_unc
  INTEGER                            :: nelft,nelf1(nobjx)
  DOUBLE PRECISION                   :: mass(nobj1x)
  LOGICAL                            :: deforb(nobj1x),defcn(nobj1x)
  CHARACTER(LEN=10)                  :: oetype
  CHARACTER(LEN=120)                 :: comele(nobj1x),elft(nifx),elf1(nifx,nobjx)

! Non-gravitational perturbations
  INTEGER                            :: nd
  
! Ephemerides
  DOUBLE PRECISION                   :: oeptim
  INTEGER                            :: lr,unit,nobj,op(10),lf,nfound, i
  LOGICAL                            :: found,opdif,opele,oepset,needrm
  
  INTEGER                            :: fail

  INTEGER lench
  EXTERNAL lench

! verbosity levels for an interactive program, but some a bit less
  verb_pro=10
  verb_clo=10
  verb_dif=10
  verb_mul=5
  verb_rej=5
 
! Fixed FORTRAN I/O units:
! report file
  unirep=1
! output orbital elements
  uniele=2
! log of difcor
  unidif=3
! ephemerides
  unieph=7

! Flags indicating whether output orbital element file and .odc file
! are to be opened or not
  opdif=.true.
  opele=.true.
  
! Initialization fail in coordinate change
  fail=0

  WRITE(*,110)
110 FORMAT(' Run name =')
  READ(*,100) run
100 FORMAT(A)
  lr=lench(run)
  WRITE(*,111) run(1:lr)
111 FORMAT(' Run name = ',A)
  
  file=run(1:lr)//'.olg'
  OPEN(unirep,FILE=file,STATUS='UNKNOWN')
  
! ========================== INITIALIZATION ==========================
  nobj=0
  CALL namini
  CALL libini
! Initializing units
  iicfil=36
  iiclsf=36
  iiceph=36
! ========================= INPUT OF OPTIONS =========================
! Read option files
  CALL filopf(unit,'orbfit.def',found)
  IF(found) THEN
     CALL rdnam(unit)
     CALL filclo(unit,' ')
  ELSE
     WRITE(*,210)
  END IF
210 FORMAT(' WARNING: no default option file (orbfit.def) found')
  optfil=run(1:lr)//'.oop'
  INQUIRE(FILE=optfil,EXIST=found)
  IF(found) THEN
     CALL filopn(unit,optfil,'OLD')
     CALL rdnam(unit)
     CALL filclo(unit,' ')
  ELSE
     lf=lench(optfil)
     WRITE(*,211) optfil(1:lf)
  END IF
211 FORMAT(' WARNING: no run-specific option file (',A,') found')

! Check of keywords
  CALL rdklst('orbfit.key')
  CALL chkkey

! Initialization of JPL ephemerides
  CALL trange

! Input of options
  CALL rmodel(1)
  CALL rdopto(run,op,name,nameo,namof,dir,nobj,elft,nelft,          &
       &            elf1,nelf1,eleout,oeptim,oepset,oetype,               &
       &            nifx,'nifx',error_model)
  CALL ofinip(run)
  IF(nobj.GT.nobjx) STOP '**** ORBFIT: nobj > nobjx ****'
! Check list of perturbing asteroids
  IF(nameo(2).ne.' ') THEN
     CALL selpert2(nameo(1),nameo(2),nfound)
  ELSE
     CALL selpert(nameo(1),found)
  END IF
! Additional options (not always required)
  needrm=.false.
! Options for initial orbit determination
  IF(op(1).GT.0)THEN
     CALL iodini
     needrm=.true.  
  ENDIF
! Options for differential correction
  IF(op(2).GT.0) THEN
     CALL rdoptf
     needrm=.true.
  END IF
! Options for identifications
  IF(op(3).GT.0) THEN
     CALL rdopti
     needrm=.true.
  END IF
! Options for ephemerides
  IF(op(4).GT.0) THEN
     CALL rdopte
     needrm=.true.
  END IF
  IF(oepset) needrm=.true.
! ====================== INPUT OF OBSERVATIONS =======================
  CALL ofiobs(unirep,name,namof,dir,nobj,n,nt,ip1,ip2,obs,obsw,error_model)
! Find if there are radar data
  DO i=1,nobj                                          
     CALL radar_ob(obs(1:n(i))%type,n(i)) 
  END DO
! ==================== INPUT OF ORBITAL ELEMENTS =====================
  IF(op(1).EQ.0)THEN
     CALL ofiorb(unirep,elft,nelft,elf1,nelf1,name,nameo,nobj,elem,elem_unc,deforb,defcn,mass,comele,nifx,nd)
  ELSE
     nd=6+nls
  END IF
! =================== INITIAL ORBIT DETERMINATION ====================
  IF(op(1).GT.0) THEN
     WRITE(*,*) 'ORBFIT: initial orbit determination'
     CALL ofinod(unirep,eleout,op(1),name,namof,deforb,defcn,dir,nobj,obs,obsw,n,nt,&
          & ip1,elem,comele,error_model,nd)
  END IF
! Change coordinates in EQU
  DO i=1,nobj
     CALL coo_cha(elem(i),'EQU',elem(i),fail)
     IF(fail.GT.0) THEN
        WRITE(*,*) ' ORBFIT: wrong coordinate change (in EQU) ', elem(i)
        STOP
     END IF
  END DO
! =========== LEAST SQUARES ORBITAL FIT (SEPARATE ORBITS) ============
  IF(op(2).GT.0) THEN
     WRITE(*,*) 'ORBFIT: least squares orbital fit'
     IF(opdif) THEN
        file=run(1:lr)//'.odc'
        OPEN(unidif,FILE=file,STATUS='UNKNOWN')
        opdif=.false.
     END IF
     CALL ofofit(unirep,uniele,unidif,opele,eleout,op(2),op(5),name,namof,deforb,defcn,elem,elem_unc,  &
          & comele,dir,nobj,obs,obsw,n,nt,ip1,ip2,error_model,oeptim,oepset,oetype,nd)
  END IF

! ======================= ORBIT IDENTIFICATION =======================
  IF(op(3).GT.0) THEN
     WRITE(*,*) 'ORBFIT: orbit identification'
     IF(opdif) THEN
        file=run(1:lr)//'.odc'
        OPEN(unidif,FILE=file,STATUS='UNKNOWN')
        opdif=.false.
     END IF
     CALL ofiden(unirep,uniele,unidif,opele,eleout,op(5),name,namof,deforb,defcn,elem,elem_unc,     &
          & mass,comele,dir,nobj,obs,obsw,n,nt,ip1,ip2,gms,oetype,error_model,nd)
  END IF
      
! ================= PROPAGATION OF ORBITAL ELEMENTS ==================
  IF(oepset) THEN
     IF(.NOT.opele) CLOSE(uniele)
     opele=.true.
     WRITE(*,*) 'ORBFIT: propagation of orbital elements'
     CALL ofprop(unirep,uniele,opele,eleout,name,deforb,defcn,elem,elem_unc,mass,comele,nobj,gms, &
          & oeptim,oepset,oetype,nd)
  END IF

! ============================ EPHEMERIDES ===========================
  IF(op(4).GT.0) THEN
     WRITE(*,*) 'ORBFIT: generating ephemerides'
     file=run(1:lr)//'.oep'
     OPEN(unieph,FILE=file,STATUS='UNKNOWN')
     print *, 'unieph:', unieph
     print *, 'name:', name
     print *, 'deforb - is the orbit defined for each obect (usually T T T):', deforb
     print *, 'defcn - are covariance/normal matrices defined (usually F T T? ... only 1 object so maybe always F?):', defcn
     print *, 'elem:', elem
     print *, 'elem_unc:', elem_unc
     print *, 'mass:', mass
     print *, 'comele:', comele
     print *, 'nobj -- object number:', nobj
     print *, 'nobj1x -- object array size (usually 3):', nobj1x
     CALL ofephe(unieph,name,deforb,defcn,elem,elem_unc,mass,comele,nobj)
  END IF
  
  CALL ofclrf
  
END PROGRAM ephem
