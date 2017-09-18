! **********************************************************
! Compute the secular evolution of dmintil
! **********************************************************
! written by C. Tardioli, 2011
! ----------------------------------------------------------
PROGRAM dmint_evol
  USE critical_points
  USE fund_const
  USE orbit_elements
  USE planet_orbits
  USE output_control
  IMPLICIT NONE
! ==================================================================
! -------- - ASTEROID ------------------
  INTEGER :: astnum !asteroid counter
  INTEGER,PARAMETER :: nastmax=100000 !maximum no. of asteroid
  CHARACTER*9 :: name !asteroid name
! --------- OPTIONS ---------------------
  CHARACTER*60 :: elefil,eledir !file names
  INTEGER :: le,lnam,iunlst
  CHARACTER*6 :: progna 
  CHARACTER*80 :: run,moid_name
  CHARACTER*3 :: ext !extension file : ext='dat' or ext='fil'
! ------ INPUT ORBITAL EVOLUTION --------
  INTEGER, PARAMETER :: secmax=100000
  REAL(KIND=dkind),DIMENSION(secmax,6) :: evolast
  REAL(KIND=dkind) :: tevol(secmax),dn(secmax),tast0,tcur
  REAL(KIND=dkind) :: aa !asteroid semimajor axis 
  INTEGER :: nevol
  CHARACTER*3 :: coo
 ! ------- DMINTIL SECULAR EVOLUTION ------
  TYPE(orbit_elem) :: elea,elem
  REAL(KIND=dkind) :: dmintil(nminx),dmintil2(nminx),deriv(nminx)
  REAL(KIND=dkind) :: lambda,step
  INTEGER :: nummin,i !i=index first interval extremum 
  REAL(KIND=dkind),PARAMETER :: tlen=300.d0 ! length of time interval (yr)
  INTEGER,PARAMETER :: nstep=100 ! number of time steps
! ----- LOOP INDEX, CONTROLS AND LOGICAL IUNIT -----------
  INTEGER :: loop,j,h,k,fail_flag,iunmoid
! --------------------------------------------------------
! planet data                                                       
  INCLUDE 'pldata.h90'
! =================================================================
  CALL pladat
! input options
  progna='propne' 
  run='propneo'
  CALL optpro(progna,run,iunlst,eledir,ext)
  astnum=0
100  FORMAT(a9) 
! =========== ENTER MAIN LOOP ===========================           
! ----- Reading asteroid orbital data from file eq1 ------
  DO 1 loop=1,nastmax
     WRITE(*,*) 'Name (Ctrl-D to quit)?' 
5    READ(*,100,end=999) name 
!!     5       READ(iunlst,100,end=999) name                             
     CALL rmsp(name,lnam) 
     IF(name(1:1).eq.'!')THEN 
        WRITE(*,*) name,'Commented near line ',astnum 
        GOTO 5 
     ENDIF
     astnum=astnum+1
     WRITE(*,*)'***************************************' 
     WRITE(*,*)'  Processing ',name,' number ',astnum 
     WRITE(*,*)'***************************************'
! =======================================================
! Planet ephemerides, coo=EQU
! Asteroid evolution : coo=KEP if ext=mev; coo=EQU if ext=fil/dat
! asteroid initial ephemerides epoch = tast0 (MJD)
     write(*,*)'planets on circular orbits',force_circ
     IF(ext.eq.'sec')THEN
        write(*,*)ext
        WRITE(*,*) 'Reading filtered planets ephemerides'
        CALL read_pla('fil')
        WRITE(*,*) 'Reading asteroid secular evolution'
        CALL read_astsec(name,tast0,aa,tevol,evolast,dn,nevol,coo)
        write(*,*)'initial epoch=',tast0
      ELSEIF(ext.eq.'fil'.OR.ext.eq.'dat')THEN
        write(*,*)ext
        WRITE(*,*) 'Reading filtered planets ephemerides'
        CALL read_pla(ext)
        WRITE(*,*) 'Reading filtered asteroid evolution' 
        CALL read_vast(name,ext,tast0,tevol,evolast,nevol,coo) 
        aa=evolast(1,1)
     ENDIF
! ======================================================= 
! Open output files
     CALL filnam('.',name,'dmin',elefil,le)
     write(*,*)'dmint_evol: opening output file ',elefil(1:le)
     CALL filopn(iunmoid,elefil(1:le),'unknown')
     write(iunmoid,103)name,tast0
     write(iunmoid,104)'time','nummin','dmintil(1:nummin)','deriv(1:nummin)'
103  FORMAT(a11,3x,f13.6)
104  FORMAT(6x,a4,2x,a6,3x,a20,10x,a20)
! =======================================================
! Asteroid elements
     elem=undefined_orbit_elem
     elem%coo=coo
     elem%coord(1)=aa
     elem%coord(6)=0.d0
! Earth elements
     elea=undefined_orbit_elem
     elea%coo='EQU'
     elea%coord(6)=0.d0
! =======================================================
! Computing dmintil interpolating tevol time
     step=tlen/float(nstep)
     write(*,*)'step = ',step
     DO 4 h=0,nstep
! Asteroid epoch for interpolation
        elem%t=tast0+(tevol(1)+h*step)*365.25d0
!        elem%t=tast0+tevol(h+1)*365.25d0
        tcur=(elem%t-tast0)/365.25d0
        write(*,'(f10.3)')tcur
        IF(tcur.ge.tevol(nevol))THEN
           WRITE(*,*) 'ENDED sec file',tcur
           GOTO 21
        ENDIF
! Searching for the time interval which contains elem%t
        CALL bin_search_inf(tcur,nevol,tevol,i)
!           write(*,*)'time = ',tcur,' spep=',h
!           write(*,*)'time interval: ',tevol(i),tevol(i+1)
! Interpolating parameter
        lambda = ((elem%t-tast0)/365.25d0-tevol(i))/(tevol(i+1)-tevol(i))
! check lambda
        IF(lambda.lt.0.d0.AND.lambda.ge.-1.d-10)THEN
           lambda=0.d0
        ELSEIF(lambda.gt.1.d0.AND.lambda.le.1.d0+1.d-10)THEN
           lambda=1.d0
        ELSEIF(lambda.lt.-1.d-10.OR.lambda.gt.1.d0+1.d-10)THEN
           WRITE(*,*) 'moidevol: error! lambda = ',lambda
           STOP
        ENDIF
! Interpolated asteroid elements at epoch elem%t
        IF(elem%coo.eq.'KEP')THEN
           elem%coord(2)=(1.d0-lambda)*evolast(i,2)+lambda*evolast(i+1,2)
           DO j=3,5
              CALL ang_interp_var(lambda,evolast(i,j),evolast(i+1,j), &
                   & elem%coord(j))
           ENDDO
        ELSEIF(elem%coo.eq.'EQU')THEN
           elem%coord(2:5)=(1.d0-lambda)*evolast(i,2:5)+ &
                & lambda*evolast(i+1,2:5)
        ELSE
           WRITE(*,*)'Asteroid coordinate not valid',elem%coo
           STOP
        ENDIF
! Earth EQU elements at the same epoch
        elea%t=elem%t
        CALL planet_elems(3,elea%t,elea%coord(1:5))
! Computing dmintil and its derivatives
        IF(ext.eq.'sec')THEN
           CALL dmintil_deriv(elea,elem,nummin,dmintil,deriv)
        ELSE
           CALL dmintil_rms(elea,elem,nummin,dmintil)
        ENDIF
! assign dummy values, to write output file
        IF(nummin.lt.2)THEN
           dmintil(2)=0.d0
           deriv(2)=0.d0
        ENDIF
        IF(ext.eq.'sec')THEN
           WRITE(iunmoid,105) tcur,nummin,dmintil(1:2),deriv(1:2)
        ELSE
           WRITE(iunmoid,106) elem%t/365.25d0,nummin,dmintil(1:2)
        ENDIF
105     FORMAT(f10.2,3x,i1,1x,4(1x,f15.8))
106     FORMAT(f10.3,2x,i1,2(2x,f15.8))
4    ENDDO
! =======================================================
21   CONTINUE
! =========== CLEAN UP ERR AND CLO FILES ================
     CALL filclo(iunmoid,' ')
! =========== END OF MAIN LOOP ==========================        
1 ENDDO
  WRITE(*,*)'insufficient number of iterations, nastmax=',nastmax
  STOP
999 CONTINUE
  WRITE(*,*) 'Number of orbit processed = ',astnum
! =========== TERMINATE CLEANLY ========================= 
!  CALL filclo(iunder,' ')

END PROGRAM dmint_evol

SUBROUTINE bin_search_inf(a,ntot,list,indnum)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: a ! interval to search for in the list
  INTEGER,INTENT(IN) :: ntot! ,i1,i2 ! number of elements in the list
! pointer to the first and last to be considered
  DOUBLE PRECISION,INTENT(IN) :: list(ntot) ! list of names
!  INTEGER,INTENT(IN) :: ordind(ntot) ! ordered addresses of the names
  INTEGER, INTENT(OUT) :: indnum ! address of the max of weak minorants 
                        ! in the list, sorted as indicated by ordind 
! -------------- end interface ----------------------------
  INTEGER :: nlist,fst,nmez,i1,i2
  INTEGER :: i ! loop index
! =========================================================
  nlist=ntot
  fst=1 ! firt extremum of the search interv
  DO i = 1,ntot
     IF(nlist.eq.1) THEN
        IF(fst.eq.ntot)THEN
           IF(a.ge.list(fst))THEN
              indnum=fst
              EXIT
           ELSE
              indnum=fst-1 ! number is just before 
           ENDIF
        ELSEIF(a.ge.list(fst).and.a.lt.list(fst+1))THEN
           indnum=fst
        ELSE
           indnum=fst-1 ! number is just before 
        ENDIF   
        EXIT
     ENDIF
     nmez=nlist/2
     IF(a.lt.list(fst+nmez)) THEN
        fst=fst
        nlist=nmez        
     ELSE
        fst=fst+nmez
        nlist=nlist-nmez        
     ENDIF
  ENDDO
END SUBROUTINE bin_search_inf

! ===================================================================   
! OPTPRO                                                                
! ===================================================================   
! input options, for the propagator and the specific main program       
! input: progna = program name (6 characters)                           
!        run    = run identifier (80 characters, up to 76 non-blank)    
SUBROUTINE optpro(progna,run,iunlst,eledir,ext)
  IMPLICIT NONE 
  character*6 progna 
  character*80 run
  integer iunlst 
  character*(*) eledir
  CHARACTER*3 :: ext !extension file : ext='dat' or ext='fil'
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
! ============= OPEN ASTEROID NAME LIST ==================              
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
! ========= WHERE TO FIND/PUT ASTEROID ELEMENTS ==========
  ireq=.true. 
  CALL rdncha(prognp,'eledir',eledir,ireq,found,fail1,fail)
  IF(fail)THEN 
     WRITE(*,*)' optpro: asteroid elements directory required' 
     STOP 
  ENDIF
! ============ EXTENSION INPUT FILES =================
! sec -> secular evolution; fil -> evolution by filtered data;
! dat -> evolution by pure numerical integration
  WRITE(*,*)'INPUT ext (sec|fil|dat):'
  READ(*,*)ext
!  ireq=.true.
!  CALL rdncha(prognp,'ext',ext,ireq,found,fail1,fail)
!  IF(fail)THEN 
!     WRITE(*,*)' optpro: extension sec/fil/dat required' 
!     STOP 
!  ENDIF
END SUBROUTINE optpro

SUBROUTINE read_vast(nam,ext,tephast0,tevol,evolast,nevol,coo)
  USE fund_const,ONLY:dkind
  USE orbit_elements
  IMPLICIT NONE
  INTEGER,PARAMETER :: secmax=100000
! -------------- interface -----------------------
  CHARACTER*9,INTENT(IN) :: nam
  CHARACTER*3,INTENT(IN) :: ext
  REAL(KIND=dkind),INTENT(OUT) :: tephast0
  REAL(KIND=dkind),DIMENSION(secmax),INTENT(OUT) :: tevol
  REAL(KIND=dkind),DIMENSION(secmax,6),INTENT(OUT) :: evolast
  INTEGER,INTENT(OUT) :: nevol
  CHARACTER*3,INTENT(OUT) :: coo
! ------------ end interface --------------------- 
  CHARACTER*9 :: name
  INTEGER :: iun,status,i,le,lnam
  INTEGER :: ngf(secmax)
! To read the header
  REAL(KIND=dkind) :: t0
  CHARACTER*60 :: commen,elefil
! *************************************************
  coo='EQU'
  name=nam
  CALL rmsp(name,lnam) 
  CALL filnam('./'//name(1:lnam),'vast',ext,elefil,le)
  write(*,*)'read_vast: opening file ',elefil(1:le)
  CALL filopn(iun,elefil(1:le),'old')

! Reading the header
  CALL skip(iun,1) !skipping asterid name
  CALL reaflc(iun,'t0',t0,commen) !reference time and comment 
  tephast0=t0-2400000.5d0
  CALL skip(iun,2)

! Reading evolution
  DO i=1,50000
     READ(iun,*,END=333) tevol(i) ! difference of time (yr) from tephast0
     READ(iun,100,END=333) evolast(i,1:6),ngf(i)
  ENDDO
  write(*,*)'do loop too short!'
  STOP
100  FORMAT(f12.9,4f11.7,f11.8,1x,i9,1p,e12.4)
333 CONTINUE
  nevol=i-1

  CALL filclo(iun,' ')

END SUBROUTINE read_vast


SUBROUTINE read_astsec(nam,tast0,aa,tevol,evolast,dn,nevol,coo)
  USE fund_const
  USE orbit_elements
  IMPLICIT NONE
  INTEGER,PARAMETER :: secmax=100000
! -------------- interface -----------------------
  CHARACTER*9,INTENT(IN) :: nam
  REAL(KIND=dkind),INTENT(OUT) :: tast0,aa
  REAL(KIND=dkind),DIMENSION(secmax),INTENT(OUT) :: tevol,dn
  REAL(KIND=dkind),DIMENSION(secmax,6),INTENT(OUT) :: evolast
  INTEGER,INTENT(OUT) :: nevol
  CHARACTER*3,INTENT(OUT) :: coo
! ------------ end interface --------------------- 
  CHARACTER*9 :: name
  CHARACTER*60 :: elefil
  INTEGER :: iunsec,le,lnam,nit(secmax),nnod(secmax),i
! ************************************************ 
  coo='KEP'
  name=nam
  CALL rmsp(name,lnam) 
!  CALL filnam(name(1:lnam),name(1:lnam),'sec',elefil,le)
  CALL filnam('.',name(1:lnam),'sec',elefil,le)
  write(*,*)'read_astsec: opening file ',elefil(1:le)  
  CALL filopn(iunsec,elefil(1:le),'old')

  READ(iunsec,107) tast0,aa
  DO i=1,secmax
     READ(iunsec,100,END=333) tevol(i),evolast(i,5),evolast(i,4), &
          & evolast(i,2:3),nit(i),dn(i),nnod(i)
     evolast(i,3:5)=evolast(i,3:5)*radeg
     evolast(i,1) = aa
  ENDDO
  WRITE(*,*)'insufficient number of iterations, secmax=',secmax
  STOP
107 FORMAT(f13.6,1x,f18.15)
100 FORMAT(f10.2,1x,f11.5,1x,f11.5,1x,f10.7,1x,f13.7,i3,1x,f13.10,1x,i3)
333 CONTINUE

  nevol=i-1

  CALL filclo(iunsec,' ')

END SUBROUTINE read_astsec
