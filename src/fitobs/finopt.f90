! ===================================================================   
! FINOPT                                                                
! ===================================================================   
! input options, for the propagator and the specific main program       
! input: progna = program name (6 characters)                           
!        run    = run identifier (80 characters, up to 76 non-blank)    
SUBROUTINE finopt(progna,run,astna0,astnap,error_model) 
  USE output_control
  USE tp_trace, ONLY: tpplane
  USE eval_risk, ONLY: givenmass, massgiven
  USE fund_const
  USE astrometric_observations, ONLY: radius, ons_name
  USE offlov_checktp, ONLY: shrinkea
  USE cobweb, ONLY: hmax,sigx,ndir,np,ndir2,np2
  USE dyn_param, ONLY : ngr_opt
  USE multi_store, ONLY : prob_sampl
  USE cla_store, ONLY : dmeacontr
  implicit none 
  character*6,INTENT(IN) :: progna 
  character*80,INTENT(IN) :: run 
  CHARACTER*(*),INTENT(OUT) :: astna0,astnap 
  CHARACTER*(20),INTENT(OUT) :: error_model ! error model file name
! hidden output through bizset: controls for bizarre orbits
  DOUBLE PRECISION ecclim,samin,samax,phmin,ahmax,qmax 
! ==========END INTERFACE============================================   
  integer le
  character*90 file 
  LOGICAL ireq,found
  CHARACTER*60 comment 
  CHARACTER*100 filnam
  INTEGER iunout, iuncovar ! now only local names
! modification to input mass from physical observations
  INCLUDE 'parlib.h90' 
  CHARACTER*200 filmass
  LOGICAL ok
  INTEGER iunmass
! ===================================================================== 
  CALL initopt(progna,run,'fop')                                         
! ===================================================================== 
! Output files: for control and results, for covariance  
! === WARNING === 
! File fou to be opened here because rmodel could write in it
  filnam=run//'.fou' 
  CALL rmsp(filnam,le) 
  CALL filopn(iunout,filnam,'UNKNOWN')
  iun_log=iunout ! later iunout to be abolished
  filnam=run//'.fga' 
  CALL rmsp(filnam,le) 
  CALL filopn(iuncovar,filnam,'UNKNOWN')
  iun_covar=iuncovar
! ===================================================================== 
! read option for physical model and integration method                 
  CALL rmodel(1) 
! ===================================================================== 
! initializations for Gauss method                                      
  CALL iodini 
! ===================================================================== 
! Output files: for errors, close approaches, for propagator parameters 
  filnam=run//'.err' 
  CALL rmsp(filnam,le)
  CALL filopn(ierrou,filnam,'UNKNOWN') 
  numerr=0 
  filnam=run//'.clo' 
  CALL rmsp(filnam,le)
  CALL filopn(iuncla,filnam,'UNKNOWN') 
  numcla=0
  filnam=run//'.pro' 
  CALL rmsp(filnam,le)
  CALL filopn(ipirip,filnam,'UNKNOWN')
! =============================
! asteroid name(s)                                                      
! =============SOME ASTEROID NAME NEEDED===================
  ireq=.false.
  astna0=run
  comment='first arc asteroid name'
  CALL input_cha_opt(progna,'astna0',astna0,ireq,found,comment,iunout)
  CALL rmsp(astna0,le)
  ons_name=.false.
  comment='not a designation'
  CALL input_log_opt(progna,'ons_name',ons_name,ireq,found,comment,iunout)
  ireq=.false. 
  astnap=' '
  comment='second arc asteroid name'
  CALL input_cha_opt(progna,'astnap',astnap,ireq,found,comment,iunout)
  CALL rmsp(astnap,le)
! =============SELECTED ERROR MODEL===================   
  error_model=' '
  ireq=.false.
  comment='error model file'
  CALL input_cha_opt(progna,'error_model',error_model,ireq,found,comment,iunout)
! =============FOR COBWEB ============================
  comment='number of columns/directions'
  CALL input_int_opt(progna,'cob_ndir',ndir,ireq,found,comment,iunout)
  comment='number of rows/points in ray'
  CALL input_int_opt(progna,'cob_np',np,ireq,found,comment,iunout)
  comment='number of columns/directions in the densified grid'
  CALL input_int_opt(progna,'ndir_grid',ndir2,ireq,found,comment,iunout)
  comment='number of rows/points in ray in the densified grid'
  CALL input_int_opt(progna,'np_grid',np2,ireq,found,comment,iunout)
  comment='sigma max in cobweb'
  CALL input_rea_opt(progna,'cob_sigx',sigx,ireq,found,comment,iunout)
  comment='max absolute magnitude'
  CALL input_rea_opt(progna,'cob_hmax',hmax,ireq,found,comment,iunout)
! ======selection of target plane=====================
!  tpplane is set to false in rmodel; true is default in fitobs
  tpplane=.true.
  comment='use b-plane as target plane' 
  CALL input_log_opt(progna,'tpplane',tpplane,ireq,found,comment,iunout)
! ======is mass given in a separate file?=====================
  massgiven=.false.
  comment='mass given in a separate file'
  CALL input_log_opt(progna,'massgiven',massgiven,ireq,found,comment,iunout)
  IF(massgiven)THEN
! get mass known from other sources, if available
     filmass=dlibd//'/'//astna0//'.mass'
     CALL rmsp(filmass,le)
     INQUIRE(file=filmass,exist=ok)
     IF(ok)THEN
        CALL filopn(iunmass,filmass(1:le),'old')
        READ(iunmass,*) givenmass, radius
        massgiven=.true.
        radius=radius/aukm ! convert from km to AU
        CALL filclo(iunmass,' ')
     ELSE
        WRITE(*,*) filmass(1:le),'not found, no mass data'
        massgiven=.false.
        STOP
     ENDIF
  ENDIF
! =============bizarre control=======================   
! ================limit eccentricity============================            
  ireq=.false. 
  ecclim=0.d0  ! leave at the bizset default
  comment='limit eccentricity'
  CALL input_rea_opt(progna,'ecclim',ecclim,ireq,found,comment,iunout)
! ================limit semimajor axis============================            
  ireq=.false. 
  samax=0.d0 !leave at the bizset default
  comment='max semimajor axis'
  CALL input_rea_opt(progna,'samax',samax,ireq,found,comment,iunout)
  ireq=.false. 
  samin=0.d0 !leave at the bizset default
  comment='min semimajor axis'
  CALL input_rea_opt(progna,'samin',samin,ireq,found,comment,iunout)
! ================limit perihelion, aphelion============================
  ireq=.false. 
  ahmax=0.d0 !leave at the bizset default
  comment='max aphelion'
  CALL input_rea_opt(progna,'ahmax',ahmax,ireq,found,comment,iunout)
  ireq=.false. 
  phmin=0.d0 !leave at the bizset default
  comment='min perihelion'
  CALL input_rea_opt(progna,'phmin',phmin,ireq,found,comment,iunout)
  qmax=0.d0 !leave at the bizset default
  comment='max perihelion'
  CALL input_rea_opt(progna,'qmax',qmax,ireq,found,comment,iunout)
! ================control on bizarre orbits================ 
  CALL bizset(ecclim,samin,samax,phmin,ahmax,qmax) 
! shrink of Earth cross section
  ireq=.false. 
  shrinkea=1.d0
  comment='shrink of Earth cross section'
  CALL input_rea_opt(progna,'shrinkea',shrinkea,ireq,found,comment,iunout)
  IF(shrinkea.GT.1.d0)THEN
     WRITE(*,*)'resintpopt: shrinkea must be smaller than 1, shrinkea=',shrinkea
     STOP
  ENDIF
! =================use of IP sampling/sigma_lov sampling===================
  ireq=.false.
  prob_sampl=.false.
  comment='use of constant steps in IP'
  CALL input_log_opt(progna,'prob_sampl',prob_sampl,ireq,found,comment,iunout)
! ================selection of smaller disk (for total ignore of close app)=====
  ireq=.false.  
  comment='max dist. for inclolinctp'
!  dmeacontr=dmea
  CALL input_rea_opt(progna,'dmeacontr',dmeacontr,ireq,found,comment,iunout)
END SUBROUTINE finopt
