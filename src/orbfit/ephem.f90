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

  INCLUDE 'parlib.h90' 

  character(len=*), parameter :: version = '1.0'
  character(len=32) :: arg
  character(len=200) :: oefilepath
  character(len=8) :: date
  character(len=10) :: time
  character(len=5) :: zone
  logical :: do_time = .false.
  INTEGER iobs
  INTEGER lenld

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
!INTEGER, PARAMETER :: nxinbl_large=800, nbl_large=40 ! 2*max no obs, number of exceptional cases
  INTEGER, PARAMETER :: nxinbl_large=400, nbl_large=40 ! 2*max no obs, number of exceptional cases
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
  CHARACTER(LEN=2048)                  ::  stringInput

! Object names and input directories
  CHARACTER(LEN=80)                  :: name(nobj1x),nameo(nobj1x),dir(nobj1x)
  CHARACTER(LEN=80)                  :: namof(nobj1x)

! OBSERVATIONS
  INTEGER                            :: nt,n(nobjx),ln
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
  CHARACTER(LEN=3)                   :: iobs_str

! Non-gravitational perturbations
  INTEGER                            :: nd
  
! Ephemerides
  DOUBLE PRECISION                   :: oeptim
  INTEGER                            :: lr,unit,nobj,op(10),lf,nfound, i
  LOGICAL                            :: found,opdif,opele,oepset,needrm
  
  INTEGER                            :: fail

! CL-ARUGUMENTS
  CHARACTER(LEN=40)                  :: objname
  CHARACTER(LEN=15)                  :: obscode, mjd


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

! COLLECT THE CL-ARGUMENTS
  do i = 1, command_argument_count()
     call get_command_argument(i, arg)
     select case (arg)
     case ('-h', '--help')
        call print_help()
        stop
     end select
  end do

  if (command_argument_count() /= 3 .AND. command_argument_count() /= 4) then
      call print_help()
      stop
  else
    do i = 1, command_argument_count()

       call get_command_argument(i, arg)

       if (i == 1) then
          obscode = trim(arg)
       else if (i == 2) then
          mjd = trim(arg)
       else if (i == 3) then
          objname = trim(arg)
       end if
    end do
  end if

! CONVERT THE OBSCODE
  call statcode(obscode,iobs)
    
! ========================== INITIALIZATION ==========================
  nobj=0
  CALL namini
  CALL libini
! Initializing units
  iicfil=36
  iiclsf=36
  iiceph=36
! ========================= INPUT OF OPTIONS =========================
! CODE FALLS OVER IF I REMOVE THIS
  CALL filopf(unit,'nofile.def',found)

! ADD COMMAND-LINE ARGUMENTS TO SETTINGS
  stringInput = "ephem.epoch.start = MJD "// mjd //" UTC"
  CALL rdstropt(stringInput) 
  stringInput = "ephem.epoch.end = MJD "// mjd //" UTC"
  CALL rdstropt(stringInput) 

  write (stringInput, "(A17,I4)") "ephem.obscode =  ", iobs
  CALL rdstropt(stringInput) 
  
  stringInput = "object1.name =  " // objname 
  CALL rdstropt(stringInput) 
  stringInput = "object1.inc_name =  " // objname
  CALL rdstropt(stringInput) 

  if(SCAN(objname, "'") /= 0) then
    stringInput = "object1.name =  " // objname 
    CALL rdstropt(stringInput) 
    stringInput = "object1.inc_name =  " // objname
    CALL rdstropt(stringInput) 
  else
    stringInput = "object1.name =  '" // objname // "'"
    CALL rdstropt(stringInput) 
    stringInput = "object1.inc_name =  '" // objname // "'"
    CALL rdstropt(stringInput) 
  end if 

  lenld=lench(dlibd)

  if  (command_argument_count() == 3) then
    stringInput = "object1.inc_files = " // dlibd(1:lenld) // "/astorb.dat[BA2]"
    CALL rdstropt(stringInput) 
  else if (command_argument_count() == 4) then
      do i = 1, command_argument_count()
         call get_command_argument(i, oefilepath)
         if (i == 4) then
            stringInput = "object1.inc_files = " // trim(oefilepath) // "[BA2]"
            CALL rdstropt(stringInput) 
         end if
      end do
  end if
  ! stringInput = "object1.inc_files = " // dlibd(1:lenld) // "/astorb.dat[BA2]"
  ! CALL rdstropt(stringInput) 


! HARDWIRED SETTINGS
  stringInput = "operations.init_orbdet = 0"
  CALL rdstropt(stringInput) 
  stringInput = "operations.diffcor = 0"
  CALL rdstropt(stringInput) 
  stringInput = "operations.ident = 0"
  CALL rdstropt(stringInput) 
  stringInput = "operations.ephem = 1"
  CALL rdstropt(stringInput) 
  stringInput = "object1.obs_fname = dummy"
  CALL rdstropt(stringInput) 
  stringInput = "object1.obs_dir = ."
  CALL rdstropt(stringInput) 
  
  stringInput = "ephem.objects = 1 2 3"
  CALL rdstropt(stringInput) 
  stringInput = "ephem.step = 1.0"
  CALL rdstropt(stringInput) 
  stringInput = "ephem.appmot.format = 'rectangular'"
  CALL rdstropt(stringInput) 
  stringInput = "ephem.appmot.units = arcsec/h"
  CALL rdstropt(stringInput) 
  stringInput = "propag.output_des = false"
  CALL rdstropt(stringInput) 
  stringInput = "propag.ab_mag = true"
  CALL rdstropt(stringInput) 
  stringInput = "propag.iast = 0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.ilun = 1"
  CALL rdstropt(stringInput) 
  stringInput = "propag.imerc = 1"
  CALL rdstropt(stringInput) 
  stringInput = "propag.iplut = 1"
  CALL rdstropt(stringInput) 
  stringInput = "propag.irel  = 1"
  CALL rdstropt(stringInput) 
  stringInput = "propag.filbe = 'CPV' "
  CALL rdstropt(stringInput) 
  stringInput = "propag.iclap = 1"
  CALL rdstropt(stringInput) 
  stringInput = "propag.iaber = 1"
  CALL rdstropt(stringInput) 
  stringInput = "propag.istat = 1"
  CALL rdstropt(stringInput) 
  stringInput = "propag.npoint= 100"
  CALL rdstropt(stringInput) 
  stringInput = "error_model.name='fcct14'"
  CALL rdstropt(stringInput) 
  stringInput = "output.elements='EQU'"

  CALL rdstropt(stringInput) 
  stringInput = "ephem.timescale = UTC"
  CALL rdstropt(stringInput) 
  stringInput = "ephem.fields = mjd,coord,mag,delta,r,elong,phase,glat,appmot"
  CALL rdstropt(stringInput) 
  stringInput = "propag.ngr_opt=.FALSE."
  CALL rdstropt(stringInput) 
  stringInput = "propag.iyark=0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.det_drp=0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.ipa2m=0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.drpa2m=0.d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.A2=0.d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.yark_exp=2.d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.det_outgas=0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.ioutgas=0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.a1ng=0.d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.a2ng=0.d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.a3ng=0.d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.dtdelay=0.d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.iyarpt=0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.yardir='.'"
  CALL rdstropt(stringInput) 
  stringInput = "propag.sep_viol=.false."
  CALL rdstropt(stringInput) 
  stringInput = "propag.eta_sep=0.d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.imet=0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.deltos=1.d-9"
  CALL rdstropt(stringInput) 
  stringInput = "propag.error=1.d-13"
  CALL rdstropt(stringInput) 
  stringInput = "propag.iord=8"
  CALL rdstropt(stringInput) 
  stringInput = "propag.hms=10.d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.hmax_me=6.d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.epms= 1.0d-12"
  CALL rdstropt(stringInput) 
  stringInput = "propag.iork= 8"
  CALL rdstropt(stringInput) 
  stringInput = "propag.eprk= 1.00d-10"
  CALL rdstropt(stringInput) 
  stringInput = "propag.lit1= 10"
  CALL rdstropt(stringInput) 
  stringInput = "propag.lit2= 4"
  CALL rdstropt(stringInput) 
  stringInput = "propag.llev=12"
  CALL rdstropt(stringInput) 
  stringInput = "propag.hev=10.d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.eprk_r=1.0d-10"
  CALL rdstropt(stringInput) 
  stringInput = "propag.lit1_r=10"
  CALL rdstropt(stringInput) 
  stringInput = "propag.lit2_r=4"
  CALL rdstropt(stringInput) 
  stringInput = "propag.lit1_rc=10"
  CALL rdstropt(stringInput) 
  stringInput = "propag.lit2_rc=10"
  CALL rdstropt(stringInput) 
  stringInput = "propag.eprk_c=1.00d-8"
  CALL rdstropt(stringInput) 
  stringInput = "propag.lit1_c=10"
  CALL rdstropt(stringInput) 
  stringInput = "propag.iusci=0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.dmea=0.1d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.dmoon=0.d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.dmjup=0.7d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.dmast=0.02d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.dter=0.1d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.ites=4"
  CALL rdstropt(stringInput) 
  stringInput = "propag.irad=0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.amrat=0.d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.amratsec=0.d0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.itide=0"
  CALL rdstropt(stringInput) 
  stringInput = "propag.ipla=2"
  CALL rdstropt(stringInput) 
  stringInput = "propag.sofa_tdb=.F."
  CALL rdstropt(stringInput) 


! Check of keywords
  CALL rdklst('orbfit.key')
  CALL chkkey

! Initialization of JPL ephemerides
  CALL trange


! Input of options
  CALL rmodel(1)
  CALL rdopto_fork(op,name,nameo,namof,dir,nobj,elft,nelft,          &
       &            elf1,nelf1,eleout,oeptim,oepset,oetype,               &
       &            nifx,'nifx',error_model)

  ! CALL ofinip(run)
  IF(nobj.GT.nobjx) STOP '**** ORBFIT: nobj > nobjx ****'
! Check list of perturbing asteroids
  IF(nameo(2).ne.' ') THEN
     CALL selpert2(nameo(1),nameo(2),nfound)
  ELSE
     CALL selpert(nameo(1),found)
  END IF
! Additional options (not always required)
  needrm=.false.

! Options for ephemerides
  IF(op(4).GT.0) THEN
     CALL rdopte
     needrm=.true.
  END IF

  IF(oepset) needrm=.true.

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
     ! file=run(1:lr)//'.oep'
     ! OPEN(unieph,FILE=file,STATUS='UNKNOWN')
     ! print *, 'unieph:', unieph
     ! print *, 'name:', name
     ! print *, 'deforb - is the orbit defined for each obect (usually T T T):', deforb
     ! print *, 'defcn - are covariance/normal matrices defined (usually F T T? ... only 1 object so maybe always F?):', defcn
     ! print *, 'elem:', elem
     ! print *, 'elem_unc:', elem_unc
     ! print *, 'mass:', mass
     ! print *, 'comele:', comele
     ! print *, 'nobj -- object number:', nobj
     ! print *, 'nobj1x -- object array size (usually 3):', nobj1x
     CALL ofephe_stdout(unieph,objname,obscode,deforb,defcn,elem,elem_unc,mass,comele,nobj)
  END IF
  
  ! CALL ofclrf

contains



  subroutine print_help()
    print '(a)', 'usage:'
    print '(a)', '   ephem <obscode> <mjd> <objectName>'
    print '(a)', ''
    print '(a)', '   <obscode>:        observatory code (use 500 for geocentric)'
    print '(a)', '   <mjd>:            the modified julian date of the ephemeris you wish to generate (UTC)'
    print '(a)', '   <objectName>:     the ID of the asteroid you wish to generate an ephemeris for (MPC number or name)'
    print '(a)', '' 
    print '(a)', 'cmdline options:'
    print '(a)', ''
    print '(a)', '  -h, --help        print usage information and exit'
  end subroutine print_help

  ! DRYX: 20170925
  !  *****************************************************************    
  !  *   *    
  !  *      RDSTROPT        *    
  !  *   *    
  !  *     Reads a namelist from string and stores in common     *    
  !  *   *    
  !  *****************************************************************    
  ! 
  ! INPUT:    IUN       -  Input FORTRAN unit         
  ! 
  SUBROUTINE rdstropt(stringInput) 
    USE option_input
    USE char_str
    IMPLICIT NONE 
    INTEGER iun 
    INCLUDE 'parcmc.h90' 
    CHARACTER(LEN=2048) stringInput
    CHARACTER*(lchx) rec,rec1,rec2,key1,key2,keyt,val1,infile,defcat
    CHARACTER*(lchx) rest 
    CHARACTER*100 inpfl 
    LOGICAL opnd,error 
    INTEGER kr,ldc,lf,lr,ipc,lr2,iuna,ipu,lk,lv,ip1,ip2,lk2,i 
    INTEGER lench 
    EXTERNAL lench 
    IF(iicnam.NE.36) STOP '**** rdnam: internal error (01) ****' 
  ! Name of the input file        
    INQUIRE(iun,OPENED=opnd,NAME=infile) 
    IF(.NOT.opnd) STOP '**** rdnam: internal error (02) ****' 
    lf=lench(infile) 
    CALL chkpdf(lf,lcfx,'lcfx') 
    kr=0 
    defcat=' ' 
    ldc=0 
  1 READ(stringInput,100) rec 
  100 FORMAT(a) 
    kr=kr+1 
    rec1=rec 
    CALL rmsp(rec1,lr) 
    IF(lr.LE.0) GOTO 10 
  ! Compute length excluding comments       
    lr=lench(rec) 
  !**   IF(lr.LT.1) GOTO 1        
    ipc=INDEX(rec(1:lr),comcha) 
    IF(ipc.EQ.1) GOTO 10  
    IF(ipc.EQ.0) THEN 
       rec1=rec(1:lr) 
    ELSE 
       rec1=rec(1:ipc-1) 
       lr=lench(rec1) 
       IF(lr.LT.1) GOTO 10 
    END IF
  ! Processing of "INPUT:" special keyword  
    rec2=rec1 
    CALL norstr(rec2,lr2) 
    IF(lr2.LT.6) GOTO 4 
    IF(rec2(1:6).EQ.'INPUT:') THEN 
       CALL strcnt(rec2(7:lr2),inpfl,rest,error) 
       IF(error) GOTO 21 
       CALL filopl(iuna,inpfl) 
       CALL rdnam1(iuna) 
       CALL filclo(iuna,' ') 
       GOTO 10
    END IF
  4 CONTINUE 
  ! Keyword field       
    ipu=INDEX(rec1(1:lr),'=') 
    IF(ipu.EQ.0) THEN 
       CALL rmsp(rec1,lk) 
       CALL chkfln(lk,lckx,'keyword',kr,infile) 
       key1=rec1(1:lk) 
       val1=' ' 
       GOTO 2 
    END IF
    key1=rec1(1:ipu-1) 
    CALL rmsp(key1,lk) 
    CALL chkfln(lk,lckx,'keyword',kr,infile) 
  ! Value field         
    val1=rec1(ipu+1:) 
    CALL norstr(val1,lv) 
    CALL chkfln(lv,lcvx,'value',kr,infile) 
  2 CONTINUE 
  ! Handling of default category  
    IF(key1(lk:lk).EQ.'.') THEN 
       ldc=lk-1 
       defcat=key1(1:ldc) 
       GOTO 10
    END IF
    IF(key1(1:1).EQ.'.') THEN 
       IF(ldc.LE.0) THEN 
          WRITE(*,101) key1(1:lk),infile(1:lf),kr 
          STOP '**** rdnam: abnormal end ****' 
       END IF
       keyt=defcat(1:ldc)//key1(1:lk) 
       key1=keyt 
       lk=lk+ldc 
    END IF
  101 FORMAT(' ERROR: missing default category declaration'/  &
       &       '        ambiguous keyword "',a,'"'/   &
       &       '        (file "',a,'", record',i4,')')
  ! Check for parentheses         
    ip1=INDEX(key1,'(') 
    IF(ip1.EQ.0) THEN 
       ip2=INDEX(key1,')') 
       IF(ip2.NE.0) THEN 
          WRITE(*,102) key1(1:lk),infile(1:lf),kr 
          STOP '**** rdnam: abnormal end ****' 
       END IF
    ELSE 
       key2=key1(ip1+1:) 
       lk2=lench(key2) 
       IF(lk2.LE.0.OR.key2(lk2:lk2).NE.')') THEN 
          WRITE(*,102) key1(1:lk),infile(1:lf),kr 
          STOP '**** rdnam: abnormal end ****' 
       END IF
    END IF
  102 FORMAT(' ERROR: illegal use of parentheses'/  &
       &       '        in keyword "',a,'"'/&
       &       '        (file "',a,'", record',i4,')')
  ! Look if the key is already present in the namelist
    DO 3 i=1,nne 
       IF(key1(1:lk).EQ.keys(i)) THEN 
          vals(i)=val1 
          namif(i)=infile 
          krecnm(i)=kr 
          kuord(i)=0 
          GOTO 10 
       END IF
  3 END DO
    nne=nne+1 
    CALL chkpdf(nne,nnex,'nnex') 
    i=nne 
    keys(i)=key1 
    vals(i)=val1 
    namif(i)=infile 
    krecnm(i)=kr 
    kuord(i)=0 
    GOTO 10 
  10 CONTINUE 
    RETURN 
  ! Error messages      
  21 CONTINUE 
    WRITE(*,106) infile(1:lf),kr 
  106 FORMAT(' ERROR: illegal "INPUT:" statement'/  &
           &       '        (file "',a,'", record',i4,')')
    STOP '**** rdnam: abnormal end ****' 
  END SUBROUTINE rdstropt


! Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: February 16, 1999
! Hacked by DRYX September 26, 2017
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         R D O P T O - fork                    *
!  *                                                               *
!  *                   Read options for ORBFIT                     *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    NIFX      -  Physical dimension of arrays ELFT,ELF1
!           CNIFX     -  Real name of NIFX parameter
!
! OUTPUT:   OP        -  Selection of execution steps:
!                          1 = initial orbit determination
!                          2 = differential correction
!                          3 = identifications
!                          4 = ephemerides
!                          5 = magnitude determination
!           NAME      -  Object names
!           NAMEO     -  Names to be used for looking in orb. el. files
!           NAMOF     -  Names of observation files (w/o extension)
!           DIR       -  Observations directory
!           NOBJ      -  Number of objects
!           ELFT      -  Input orbital element files (for all objects)
!           NELFT     -  Number of input element files (for all objects)
!           ELF1      -  Input orbital element files (for each object)
!           NELF1     -  Number of input element files (for each object)
!           OEFILE    -  Output orbital element file (possibly blank)
!           OEPTIM    -  Epoch of output elements (MJD, TDT)
!           OEPSET    -  Flag stating that an output epoch is requested
!           OETYPE    -  Type of output elements (CAR/EQU/KEP/EQP)
!           ERROR_MODEL - Error model file name
!
  SUBROUTINE rdopto_fork(op,name,nameo,namof,dir,nobj,elft,nelft,    &
       &                  elf1,nelf1,oefile,oeptim,oepset,oetype,         &
       &                  nifx,cnifx,error_model)
    
    INCLUDE 'sysdep.h90'
    INTEGER,          INTENT(IN)  :: nifx
    CHARACTER(LEN=*), INTENT(IN)  :: cnifx
    CHARACTER*(80),   INTENT(OUT) :: name(2),nameo(2),namof(2) 
    CHARACTER(LEN=*), INTENT(OUT) :: dir(2), elft(nifx), elf1(nifx,2)
    INTEGER,          INTENT(OUT) :: op(5),nobj,nelft,nelf1(2)
    CHARACTER(LEN=*), INTENT(OUT) :: oefile, oetype
    DOUBLE PRECISION, INTENT(OUT) :: oeptim
    LOGICAL,          INTENT(OUT) :: oepset
    CHARACTER(LEN=*), INTENT(OUT) :: error_model 
! ***********************************************************************************
    INTEGER          :: i,lr,mjd,mjde,ln,ld
    DOUBLE PRECISION :: sec,sece
    LOGICAL          :: found,fail1,fail
    CHARACTER        :: cc*100,scale*10
    
    INTEGER lench
    EXTERNAL lench

! Required execution steps (by default, all are selected but ephemeris)
    fail=.FALSE.
    CALL rdnint('operations.','init_orbdet',op(1),.false.,            &
         &            found,fail1,fail)
    CALL rdnint('operations.','diffcor',op(2),.false.,                &
         &            found,fail1,fail)
    CALL rdnint('operations.','ident',op(3),.false.,                  &
         &            found,fail1,fail)
    CALL rdnint('operations.','ephem',op(4),.false.,                  &
         &            found,fail1,fail)
    CALL rdnint('operations.','magfit',op(5),.false.,                 &
         &            found,fail1,fail)
    
    nobj=1
! Observations directory
    IF(.NOT.fail1) THEN
       CALL rdncha('object1.','obs_dir',dir(1),.false.,              &
            &                found,fail1,fail)
    END IF
! Initial condition file
    nelf1(1)=0
    CALL rdmcha('object1.','inc_files',elf1(1,1),nelf1(1),nifx,       &
         &            cnifx,.false.,found,fail1,fail)
    IF(.NOT.found) CALL rdodin(name(1),nelf1(1),elf1(1,1),nifx)
    CALL rdncha('object1.','inc_name',nameo(1),.false.,               &
         &            found,fail1,fail)
    IF(.NOT.found) nameo(1)=name(1)
    CALL rmsp(nameo(1),ln)
    CALL rdncha('object1.','obs_fname',namof(1),.false.,              &
         &            found,fail1,fail)
    IF(.NOT.found) THEN
       namof(1)=name(1)
       CALL rmsp(namof(1),ln)
    END IF
      
! Epoch and type of output elements
    oeptim=0.d0
    CALL rdntim('output.','epoch',cc,mjd,sec,scale,.false.,           &
         &            oepset,fail1,fail)
    IF(oepset) THEN
       CALL cnvtim(mjd,sec,scale,mjde,sece,'TDT')
       oeptim=mjde+sece/86400.d0
    END IF
    
    CALL rdncha('output.','elements',oetype,.false.,                  &
         &            found,fail1,fail)
    
    IF(fail) STOP '**** rdopto_fork: abnormal end ****'
        
! Error model
    CALL rdncha('error_model.','name',error_model,.false.,found,fail1,fail)
    IF(fail) STOP '**** rdopto_fork: abnormal end ****'
    
  END SUBROUTINE rdopto_fork

! ==============================================                        
!  statcode, codestat                                                   
! conversion from/to exadecimal to/from numeric code for observing stati
! answer to mess done by MPC in April 2002                              
! ==============================================                        
  SUBROUTINE statcode(obsstr,iobs) 
  IMPLICIT NONE 
! input                                                                 
  CHARACTER*3 obsstr 
! output                                                                
  INTEGER iobs 
! end interface                                                         
  CHARACTER*1 alfanum 
  INTEGER hundreds, temp 
  LOGICAL isnum 
  READ(obsstr,100)iobs 
100 FORMAT(1x,i2) 
  READ(obsstr,'(A1)')alfanum 
  IF(isnum(alfanum))THEN 
     READ(alfanum,'(I1)')hundreds 
     iobs=iobs+100*hundreds 
  ELSEIF(alfanum.eq.' ')THEN 
     iobs=iobs 
  ELSE 
     temp=ichar(alfanum) - 55 
     IF(temp.gt.35) temp = temp - 6 
     iobs=iobs+100*temp 
  ENDIF 
  RETURN 
  END SUBROUTINE statcode     

  
END PROGRAM ephem
