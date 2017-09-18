
! -------------------------------------------------------------------------
! search for candidate potential NEAs by looking in the plane (omega,e)
! at intersections of the line e=e_* 
! (e_* is the critical value of the eccentricity defined by a(1-e_*)=1.3 AU) 
! with the level curve of the averaged Hamiltonian Hbar.
! A necessary condition for a PNEAs is that these two lines cross each other
! (it is not sufficient because a level curve could be not-connected)
! HINT! the current NEAs could never cross the critical line, because they 
! can stay always inside the NEA region 
! -------------------------------------------------------------------------
! written by G.F. Gronchi, July 2011

PROGRAM candidate_PNEAs
  USE fund_const
  USE orbit_elements
  USE perturbing_function
  USE dyn_param

  IMPLICIT NONE
  INTEGER,PARAMETER :: nastmax=1000000
  INTEGER :: loop,astnum !asteroid counters                               
  CHARACTER*9 :: name !asteroid name
  INTEGER :: lnam,le
  INTEGER :: iunin,iuncan
  CHARACTER*60 eledir,elefil
  LOGICAL :: found,eof
  TYPE(orbit_elem) :: eq,elkep
  TYPE(orb_uncert) :: unc
  INTEGER :: fail_flag
  REAL(KIND=dkind) :: a,ecc,inc,Omnod,omega
  REAL(KIND=dkind) :: estar,G,zl,Hbar
  REAL(KIND=dkind) :: ddd,eee 
  INTEGER :: nnn 
  INTEGER,PARAMETER :: nsamp=180
  INTEGER :: j ! loop index  
  LOGICAL :: err
  INTEGER :: nlsloc
! common in pldata.h90
!  INTEGER nplax,nplax2,npl,inpl,ioupl,ndum
!  parameter (nplax=10)
!  parameter (nplax2=20)
!  DOUBLE PRECISION ap(nplax),gm(nplax)
!  DOUBLE PRECISION ky,bigg 
!  COMMON/pldata/gm,ap,ky,ndum,npl,inpl,ioupl,bigg
!  INCLUDE 'pldata.h90'
! --------------------------------------------------------
  REAL(KIND=dkind) :: omj,Gj,Hj,diff
  LOGICAL :: pnea

! initialization
  eledir='epoch'
  astnum=0
  rhs=1 

  CALL filopn(iuncan,'candidates','unknown')

! ENTER MAIN LOOP
  DO loop=1,nastmax
     astnum = astnum + 1 
! READ NEXT NAME
     write(*,*) 'Name (Ctrl-D to quit)?' 
5    READ(*,100,end=999) name 
100  FORMAT(a9) 
     CALL rmsp(name,lnam) 
     IF(name(1:1).eq.'!')THEN 
        WRITE(*,*) name,'Commented near line ',astnum 
        GOTO 5 
     ENDIF
!     WRITE(*,*)'***************************************' 
!     WRITE(*,*)'  Processing ',name,' number ',astnum 
!     WRITE(*,*)'***************************************' 
! INPUT ORBITAL ELEMENTS
     CALL filnam(eledir,name,'eq0',elefil,le) 
     INQUIRE(file=elefil(1:le),exist=found)
     write(*,*) astnum,elefil 
     IF(.not.found)THEN 
        WRITE(*,*)' elements not found for ',name 
        STOP
     ENDIF
     CALL read_elems(eq,name,eof,nlsloc,err,FILE=elefil)
     !     write(*,'(5(f10.5,2x))') eq%coord(1:5)

     CALL coo_cha(eq,'KEP',elkep,fail_flag)
     IF(fail_flag.ge.5) THEN
        WRITE(*,*)'coo_cha failed'
        STOP
     ENDIF
     a=elkep%coord(1)
     ecc=elkep%coord(2)
     inc=elkep%coord(3)
     Omnod=elkep%coord(4)
     omega=elkep%coord(5)

! check if it is a current NEA
     IF(a*(1-ecc).le.1.3d0)THEN
        WRITE(*,*)'current NEA'
        WRITE(*,*)name,'q=',a*(1-ecc)
        write(iuncan,100)name
        GOTO 111
     ENDIF

     CALL pladat 

! constant Z                                                            
     zl=sqrt(a*(1.d0-(ecc**2)))*cos(inc) 
!     write(*,*)'zl=',zl
     G=ky*sqrt(a*(1.d0-ecc**2)) 
     CALL perturb(omega,G,zl,a,ddd,eee,nnn) 
     Hbar=ddd
     write(*,'(a12,1x,5(f10.5,2x))')'omega,e,Hbar',omega*180.d0/pig,ecc,Hbar

     pnea=.false.
!     omj=0.d0
     estar=1.d0-1.3d0/a ! critical value for NEA
     write(*,*)'estar',estar
     Gj=ky*sqrt(a*(1.d0-estar**2)) 
     DO j=0,nsamp-1

        omj=(j*pig)/(2.d0*nsamp)
!        omj=(j*pig)/nsamp
!        omj=(j*dpig)/nsamp
! Hint: MUST USE THE SAME VALUE OF zl !!!

!        CALL clevpert(omj,Gj,zl,a,ddd,eee,nnn) 
        CALL perturb(omj,Gj,zl,a,ddd,eee,nnn) 
        Hj=ddd
!        write(*,'(a13,1x,5(f10.5,2x))')'omega,Hj,Hbar',omj*180.d0/pig,Hj,Hbar
        IF((j.gt.0).and.(diff*(Hj-Hbar).lt.0.d0))THEN
           write(*,*)'the sign changes:',name
           write(iuncan,100)name
           pnea=.true.
           WRITE(*,'(a12,1x,5(f10.5,2x))')'omega,e,Hbar',omj*180.d0/pig,estar,Hj
           GOTO 111
        END IF
        diff=Hj-Hbar
     ENDDO
111  CONTINUE
  ENDDO
  STOP ' Warning, too many objects. Reached end of loop.' 
! =========== EOF EXIT POINT ===================================        
999 CONTINUE 

  CALL filclo(iuncan,' ')

END PROGRAM candidate_PNEAs
