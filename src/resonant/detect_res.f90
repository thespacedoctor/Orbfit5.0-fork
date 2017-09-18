PROGRAM detect_res
  USE fund_const
  USE output_control
  USE orbit_elements
  USE dyn_param, ONLY: ndimx
  USE right_hand_side, ONLY: kappap,kappaa,n_res,ffd_pert_res,rot_matrix,gmp,sin_flag
  USE planet_orbits
  IMPLICIT NONE

!  asteroid counter, total, maximum
  INTEGER astnum,nastmax,num_obj,loop
  PARAMETER (nastmax=100000)
!  asteroid name                                           
  CHARACTER*11 name,namev(1)
! ======================= state variables ==========================
!  orbital elements, epoch times input data (for rdastb.f)
!  REAL(KIND=dkind) elkep(6)
  TYPE(orbit_elem) :: elequ,elkep
  REAL(KIND=dkind) :: om,omnod,aa,ecc,Gdel,Zdel

! ========================================================
!  these quantities are not used in this program, but
!  they are needed in calling subroutine rdelem                
  REAL(KIND=dkind) :: cov0(ndimx,ndimx,1),norm0(ndimx,ndimx,1),mass,massv(1)
  REAL(KIND=dkind) :: hmag,gmag,enne
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
!  INTEGER :: sin_flag!,res_flag

! for dqags                                                         
  INTEGER :: limx,limx4 
  PARAMETER (limx=500,limx4=4*limx)
  INTEGER :: ier,limit,iwork(limx),lenw,last 
! function evaluations                                              
  INTEGER :: neval,neval_res_c,neval_res_s 
  REAL(KIND=dkind) :: abserr, abserr_res_c,abserr_res_s,work(limx4) 
! output of dqags
  REAL(KIND=dkind) :: rm
  REAL(KIND=dkind) :: np ! mean motion

  REAL(KIND=dkind),PARAMETER :: epsabs = 1.d-9,epsrel=1.d-7
  TYPE(orbit_elem) :: elpl !resonant
  TYPE(orbit_elem) :: elplar

  REAL(KIND=dkind) :: Sres,ares,alpha,beta,Acoe,Bcoe,sigmastar,Vstar,DeltaS
  INTEGER :: j

! ========================================================

! ============ input options ==================
  INCLUDE 'neoopt.h90' 
  progna='propne' 
  run='propneo'
  CALL optpro(progna,run,iunlst,eledir)

!  Output files: for control and results                             
  call filopn(iunlog,'prop.log','UNKNOWN') 
  call filopn(iunhand,'prop.names.fail','unknown') 
  astnum = 0 



! ---- preparing dqags call ---
!  epsabs=1.d-8 
!  epsrel=1.d-5 
  limit=limx 
  lenw=limx4 

  !     PLANET DATA
  CALL pladat 

! loading filtered planetary ephemerides from vpla.fil
  CALL read_pla('dat')

!----------insert resonance and resonant planet--------
  write(*,*) 'insert resonance and resonant planet ka,kp,n_res'
  read(*,*) kappaa,kappap,n_res


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
!     CALL filnam(eledir,name,'eq0',elefil,le) 
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
     CALL filnam('.',name,'sec',elefil,le) 
     CALL filopn(iunit,elefil(1:le),'unknown')
     CALL filnam('.',name,'pre',elefil,le) 
     CALL filopn(iunpre,elefil(1:le),'unknown') 


!     transformation of input elements
!     zl=ky*sqrt(elkep(1)*(1.d0-(elkep(2)**2)))*cos(elkep(3)) 
!     g=ky*sqrt(elkep(1)*(1.d0-elkep(2)**2)) 
!     check
     write(*,*)'propneo:a,e,ai,omnod,omeg:',t1,elkep%coord(1:5)

!     =========== CLEAN UP ERR AND CLO FILES ===================
2    CONTINUE
!     close error file                                                  
     IF(numerr.gt.0)THEN 
        CALL filclo(ierrou,' ') 
     ELSE 
        CALL filclo(ierrou,'DELETE') 
     ENDIF

! ============================================
! compute resonant width


! Planet keplerian elements
     el_pla(n_res)=undefined_orbit_elem
     el_pla(n_res)%t=elkep%t
     el_pla(n_res)%coo='EQU'
     !  write(*,*)'i,el_pla(i)%coord(6) default1',i,el_pla(i)%coord(6)


! interpolate for the i-th planet orbit at asteroid initial time tast0
     CALL planet_elems(n_res,el_pla(n_res)%t,el_pla(n_res)%coord(1:6))
     CALL coo_cha(el_pla(n_res),'KEP',el_pla(n_res),fail_flag)
     IF(fail_flag.ge.5)THEN
        WRITE(*,*)'error in coo_cha! fail_flag=',fail_flag
     ENDIF



! ========================================================== 
! COMPUTATION of A, B, alpha, beta
  om = elkep%coord(5)
  omnod = elkep%coord(4)

!  aa=elkep%coord(1)
  aa = el_pla(n_res)%coord(1)*(abs(kappap)/abs(kappaa))**(-2.0/3.0)
  ares = aa
  write(*,*)'aa = ',aa
  ecc = 0.d0
  DO j=1,50 ! do loop on eccentricity
!  ecc=elkep%coord(2)
     ecc = ecc+0.01
     Gdel=ky*sqrt(aa*(1.d0-ecc**2)) 
     Zdel=Gdel*cos(elkep%coord(3))
     gmp = gm(n_res)   !gmp = (k^2)*m(pla)/m(Sun)

     write(*,*)'gmp=',gmp
! compute rotational matrix and Cartesian coords for planet and asteroid
     CALL rot_matrix(el_pla(n_res),om,omnod,Gdel,Zdel,aa)

     sin_flag=0
     CALL dqags(ffd_pert_res,0.d0,dpig,epsabs,epsrel,Acoe,abserr,neval,ier, &
          &        limit,lenw,last,iwork,work)
     sin_flag=1 
     CALL dqags(ffd_pert_res,0.d0,dpig,epsabs,epsrel,Bcoe,abserr,neval,ier,&
          & limit,lenw,last,iwork,work)

     write(*,*)'Acoe, Bcoe =',Acoe,Bcoe
     
     Sres = ky/kappaa*sqrt(ares)
     
     alpha = 3*ky**4/(kappaa**2*Sres**4)
     
     beta = gmp*2*ky**2/(dpig**2)
     
     sigmastar = atan(Bcoe/Acoe)
     
     Vstar = beta*(Acoe*cos(sigmastar) + Bcoe*sin(sigmastar))

     ! resonance width
     DeltaS = 2.d0*sqrt(abs(Vstar/alpha))

     write(*,*)'alpha, beta =',alpha,beta
     write(*,*)'sigmastar, Vstar =',sigmastar, Vstar
     write(*,*)'DeltaS =',DeltaS
     write(*,*)'Sres =',Sres
     write(*,*)'S =',ky/kappaa*sqrt(elkep%coord(1))
     write(25,*)Sres,DeltaS,ecc

  ENDDO


!     =========== END OF MAIN LOOP =============================        
1 ENDDO
  STOP ' Warning, too many objects. Reached end of loop.' 
! =========== EOF EXIT POINT ===================================        
999 CONTINUE 















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

END PROGRAM detect_res
