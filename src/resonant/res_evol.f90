! ==================================================================
!  PROGRAM RES_EVOL
!  written by S. Marò and G.F. Gronchi, July 2016
! ==================================================================
!  Uses propne.nop for options.
!  Generates proper elements files from eq1 file 
! ==================================================================
PROGRAM res_evol
  USE fund_const
  USE output_control
  USE orbit_elements
  USE rkg_local, ONLY:ismax
  USE dyn_param, ONLY: ndimx
  USE right_hand_side, ONLY: kappap,kappaa,n_res
  IMPLICIT NONE
! ==================================================================
!  asteroid counter, total, maximum
  INTEGER astnum,nastmax,num_obj,loop
  PARAMETER (nastmax=100000)
!  asteroid name                                           
  CHARACTER*11 name,namev(1)
! ======================= state variables ==========================
!  orbital elements, epoch times input data (for rdastb.f)
!  REAL(KIND=dkind) elkep(6)
  TYPE(orbit_elem) :: elequ,elkep
! ========================================================
!  these quantities are not used in this program, but
!  they are needed in calling subroutine rdelem                
  REAL(KIND=dkind) cov0(ndimx,ndimx,1),norm0(ndimx,ndimx,1),mass,massv(1)
  REAL(KIND=dkind) hmag,gmag,enne
  REAL(KIND=dkind) :: hmagv(1),gmagv(1)
! ========================================================
!  all data on planets: here we use only planetary masses
  INCLUDE 'pldata.h90'
! ============= cleverpert call ===============
  REAL(KIND=dkind) g,zl ! input (om,a are given from the database)     
  REAL(KIND=dkind) hamil0,error ! output
  INTEGER numb 
! =============================================                     
  CHARACTER*3 eltype,eltypev(1) 
  CHARACTER*80 comele,comelev(1) 
!  successful input flags, for elements                              
  LOGICAL defcov,defcovv(1),defelem,defelemv(1),neodys_elem 
! ============ propagation =====================
  REAL(KIND=dkind) t1,t1v(1) ! current time (in Julian days)
  REAL(KIND=dkind) elem1(6) ! elements
! ========= input control ======================                    
!  file names                                                        
  CHARACTER*60 file,elefil,eledir,elemfiles(4) 
  INTEGER le,lnam,num_files 
  CHARACTER*6 :: progna 
  CHARACTER*80 :: run
!  logical units                                                     
  INTEGER :: iunlog,iunit,iunlst,iunhand,iunpre 
!  loop index                                                        
  INTEGER :: i
  INTEGER :: fail_flag
  REAL(KIND=dkind) :: deltat
! ========================================================

! ============ input options ==================
  INCLUDE 'neoopt.h90' 
  progna='propne' 
  run='propneo'
  CALL optpro(progna,run,iunlst,eledir)
  OPEN(2,file='prsingeq.opt',status='old')                      
  read(2,*) eprk,epnod,isrk ! integration parameters
  read(2,*)
  read(2,*) deltat ! time interval for secular evolution
  CLOSE(2) 
! =======================================================           
!  Output files: for control and results                             
  call filopn(iunlog,'prop.log','UNKNOWN') 
  call filopn(iunhand,'prop.names.fail','unknown') 
  astnum = 0 

! =========== ENTER MAIN LOOP ===========================           
  DO 1 loop=1,nastmax
! =======================================================           
     astnum = astnum + 1 
! =================READ NEXT NAME =======================           
     write(*,*) 'Name (Ctrl-D to quit)?' 
5    READ(*,100,end=999) name 
!     5       READ(iunlst,100,end=999) name                             
100  FORMAT(a11) 
     CALL rmsp(name,lnam) 
     IF(name(1:1).eq.'!')THEN 
        WRITE(iunhand,*) name,'Commented near line ',astnum 
        GOTO 5 
     ENDIF
!----------insert resonance and resonant planet--------
     write(*,*) 'insert resonance and resonant planet ka,kp,n_res'
     read(*,*) kappaa,kappap,n_res

     namev(1) = name
     defelem=.false. 
     defelemv(1)=.false.
     defcov=.false. 
     defcovv(1)=.false.
!     error file to be opened                                           
     CALL filnam('./err',name,'err',file,le) 
     CALL filopn(ierrou,file(1:le),'unknown') 
     numerr=0 
!     put flag in each log file                                         
     WRITE(iunlog,*)'***************************************' 
     WRITE(iunlog,*)' Processing ',name,' number ',astnum 
     WRITE(iunlog,*)'***************************************' 
     WRITE(*,*)'***************************************' 
     WRITE(*,*)'  Processing ',name,' number ',astnum 
     WRITE(*,*)'***************************************' 
! ================ INPUT ELEMENT SETS ====================          
     CALL filnam(eledir,name,'eq1',elefil,le) 
     INQUIRE(file=elefil(1:le),exist=neodys_elem) 
     IF(neodys_elem)THEN 
        elemfiles(1) = elefil(1:le) 
        num_files = 1 
     ELSE 
        WRITE(ierrou,*)' elements not found for ',name 
        numerr=numerr+1 
        GOTO 2 
     ENDIF
     num_obj=1 
     CALL rdelem(iunlog,namev(1),num_obj,elemfiles(1),num_files,defelemv(1),   &
     & defcovv(1),eltypev(1),t1v(1),elem1(1),cov0(1:ndimx,1:ndimx,1), &
     & norm0(1:ndimx,1:ndimx,1),massv(1),hmagv(1),gmagv(1),comelev(1))  
! most of the variables above have dimension = num_obj; 
! in this case we write 1 explicitely in its place. 
     IF(.not.defelemv(1))THEN 
        WRITE(iunlog,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
        WRITE(iunlog,*)'!!WARNING! WARNING! WARNING! WARNING!!' 
        WRITE(iunlog,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
        WRITE(iunlog,*)'Asteroid ',name,' number ',astnum,          &
 &           ' not found in asteroid catalog.'                      
        WRITE(iunhand,*)name,'not in catalogs.' 
        GOTO 2 
     ENDIF
     eltype = eltypev(1)
     t1=t1v(1)
     mass=massv(1)

! equinoctial elements of the asteroid
     elequ=undefined_orbit_elem
     elequ%coord(1:6)=elem1(1:6)
     elequ%coo=eltype
     elequ%t=t1

!     coordinate change
!     CALL coocha(elem1,eltype,bigg,elkep,'KEP',enne) 
     CALL coo_cha(elequ,'KEP',elkep,fail_flag) 
!     open output files
     CALL filnam('.',name,'res',elefil,le) 
     CALL filopn(iunit,elefil(1:le),'unknown')
     CALL filnam('.',name,'pre',elefil,le) 
     CALL filopn(iunpre,elefil(1:le),'unknown') 
!     PLANET DATA
     CALL pladat 
!     transformation of input elements
!     zl=ky*sqrt(elkep(1)*(1.d0-(elkep(2)**2)))*cos(elkep(3)) 
!     g=ky*sqrt(elkep(1)*(1.d0-elkep(2)**2)) 
!     check
     write(*,*)'propneo:a,e,ai,omnod,omeg:',t1,elkep%coord(1:5)

     write(*,*)'iunit=',iunit,'elefil',elefil(1:le)

! writing mev-file header
     WRITE(iunit,200) t1,elkep%coord(1)
200  FORMAT(f13.6,1x,f18.15)
!     ==========================================================

!     ==========================================================
!     evolution computation: elkep = (a,e,I,Om,om,l)
!     HINT!:angles are passed in radians
!     CALL prsingeq(iunit,iunpre,t1,elkep(1),elkep(2),elkep(3), &
!          & elkep(4),elkep(5),zl)
     CALL prsingeq(iunit,iunpre,elkep,deltat)
!     HINT!:output on file.mev is written by subroutine proutele        
!     ==========================================================        
     CALL filclo(iunit,' ')
     CALL filclo(iunpre,' ') 
!     ==========================================================
!     reopen file.pre to get file .pro, ordered wrt time
     write(*,*)'start sorting!!!' 
     CALL sortpre(name) 
!     =========== CLEAN UP ERR AND CLO FILES ===================
2    CONTINUE
!     close error file                                                  
     IF(numerr.gt.0)THEN 
        CALL filclo(ierrou,' ') 
     ELSE 
        CALL filclo(ierrou,'DELETE') 
     ENDIF
!     =========== END OF MAIN LOOP =============================        
1 ENDDO
  STOP ' Warning, too many objects. Reached end of loop.' 
! =========== EOF EXIT POINT ===================================        
999 CONTINUE 
!     ===== Open file as a completion flag for Perl Script =====        
  CALL filopn(iunit,'complete','unknown') 
  WRITE(iunit,*)'All done!' 
  CALL filclo(iunit,' ') 
!     ==========TERMINATE CLEANLY=========================              
  CALL filclo(iunlog,' ') 
!     CALL filclo(iunlst,' ')                                           
  CALL filclo(iunhand,' ') 
!     ====================================================  
CONTAINS

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

END PROGRAM res_evol
