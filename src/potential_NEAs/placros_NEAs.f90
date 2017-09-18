
! -------------------------------------------------------------------------
! search for number of Earth crossings by looking in the plane (omega,e)
! at intersections of the node crossing lines dnodp = dnodm = 0, where
!       dnodp = a' - a(1-e)/(1 + e*cos(omega)) ,  
!       dnodm = a' - a(1-e)/(1 - e*cos(omega)) , 
! with the level curve of the averaged Hamiltonian Hbar.
! A necessary condition for a PNEAs is that these lines cross each other
! (it is not sufficient because a level curve could be not-connected)
! -------------------------------------------------------------------------
! written by G.F. Gronchi, February 2012

PROGRAM placros_NEAs
  USE fund_const
  USE orbit_elements
  USE perturbing_function
  USE dyn_param

  IMPLICIT NONE
  INTEGER,PARAMETER :: nastmax=1000000
  INTEGER :: loop,astnum !asteroid counters                               
  CHARACTER*9 :: name !asteroid name
  INTEGER :: lnam,le
  INTEGER :: iunin,iuncan,iuncros
  CHARACTER*60 eledir,elefil
  LOGICAL :: found,eof
  TYPE(orbit_elem) :: eq,elkep
  TYPE(orb_uncert) :: unc
  INTEGER :: fail_flag
  REAL(KIND=dkind) :: a,ecc,inc,Omnod,omega
  REAL(KIND=dkind) :: ej,G,zl,Hbar,emax
  REAL(KIND=dkind) :: ddd,eee 
  INTEGER :: nnn 
  INTEGER,PARAMETER :: nsamp=500
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
  REAL(KIND=dkind) :: bigX,bigY,theta,xi,eta,rho,halfrapp
  REAL(KIND=dkind) :: omj,Gj,Hj,diff,ejprev
  LOGICAL :: fst_inside
  INTEGER :: numcros

! initialization
  eledir='epoch'
  astnum=0
  rhs=1 

  CALL filopn(iuncan,'Earth_crossings','unknown')
  CALL filopn(iuncros,'numcross','unknown')

! ENTER MAIN LOOP
  DO loop=1,nastmax
     astnum = astnum + 1 
! READ NEXT NAME
     write(*,*) 'Name (Ctrl-D to quit)?' 
5    READ(*,100,end=999) name 
100  FORMAT(a9) 
101  FORMAT(a9,1x,4(f10.5,2x),2x,a1) 
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

     CALL pladat 

! constant Z                                                            
     zl=sqrt(a*(1.d0-(ecc**2)))*cos(inc) 
!     write(*,*)'zl=',zl
     G=ky*sqrt(a*(1.d0-ecc**2)) 
     CALL perturb(omega,G,zl,a,ddd,eee,nnn) 
     Hbar=ddd
!     write(*,'(a12,1x,5(f10.5,2x))')'omega,e,Hbar',omega*180.d0/pig,ecc,Hbar
!     write(*,*)'starting loop'

     emax = sqrt(1.d0 - zl**2/a) ! radius of Kozai's domain
     halfrapp = ap(3)/(2.d0*a) 
     rho = 1.d0 - halfrapp

     numcros=0 !initialization
! sampling of the ascending node-crossing curve
     fst_inside=.false.
     ejprev=0.d0 !initialization
     DO j=0,nsamp-1
        theta = j*dpig/nsamp
        bigX = rho*cos(theta) 
        bigY = rho*sin(theta) 
        xi = bigX - halfrapp
        eta = bigY
        ej = sqrt(xi**2+eta**2)
        IF(ejprev.gt.emax.and.ej.le.emax)THEN
           fst_inside=.true.
        ELSE
           fst_inside=.false.
        ENDIF
        ejprev = ej !remember ej
! check if inside Kozai's domain
        IF(ej.gt.emax)THEN
!           WRITE(*,*)'outside Kozai domain: ej=',ej,'emax=',emax
           CYCLE
        ENDIF
        omj = ATAN2(eta/ej,xi/ej)
        Gj=ky*sqrt(a*(1.d0-ej**2)) 

! Hint: MUST USE THE SAME VALUE OF zl !!!
        CALL perturb(omj,Gj,zl,a,ddd,eee,nnn) 
        Hj=ddd
!        write(*,'(a13,1x,5(f10.5,2x))')'omega,Hj,Hbar',omj*180.d0/pig,Hj,Hbar
        IF((j.gt.0).and.(.not.fst_inside).and.(diff*(Hj-Hbar).lt.0.d0))THEN
!           write(*,*)'*** ascending node *** The sign changes:',name
           numcros=numcros+1
           write(iuncan,101)name,omj*180.d0/pig,ej,Hj,Hbar,'A'
           WRITE(*,'(a12,1x,5(f10.5,2x))')'(A) omega,e,Hbar',omj*180.d0/pig,ej,Hj
!           GOTO 111
        END IF
        diff=Hj-Hbar
     ENDDO

! sampling of the descending node-crossing curve
     fst_inside=.false.
     ejprev=0.d0 !initialization
     DO j=0,nsamp-1
        theta = j*dpig/nsamp
        bigX = rho*cos(theta) 
        bigY = rho*sin(theta) 
        xi = bigX + halfrapp
        eta = bigY
        ej = sqrt(xi**2+eta**2)
        IF(ejprev.gt.emax.and.ej.le.emax)THEN
           fst_inside=.true.
        ELSE
           fst_inside=.false.
        ENDIF
        ejprev = ej !remember ej
! check if inside Kozai's domain
        IF(ej.gt.emax)THEN
!           WRITE(*,*)'outside Kozai domain: ej=',ej,'emax=',emax
           CYCLE
        ENDIF
        omj = ATAN2(eta/ej,xi/ej)
        Gj=ky*sqrt(a*(1.d0-ej**2)) 

! Hint: MUST USE THE SAME VALUE OF zl !!!
        CALL perturb(omj,Gj,zl,a,ddd,eee,nnn) 
        Hj=ddd
!        write(*,'(a13,1x,5(f10.5,2x))')'omega,Hj,Hbar',omj*180.d0/pig,Hj,Hbar
        IF((j.gt.0).and.(.not.fst_inside).and.(diff*(Hj-Hbar).lt.0.d0))THEN
!           write(*,*)'*** descending node *** The sign changes:',name
           numcros=numcros+1
           write(iuncan,101)name,omj*180.d0/pig,ej,Hj,Hbar,'D'
           WRITE(*,'(a12,1x,5(f10.5,2x))')'(D) omega,e,Hbar',omj*180.d0/pig,ej,Hj
!           GOTO 111
        END IF
        diff=Hj-Hbar
     ENDDO

     IF(numcros.gt.0)THEN
        WRITE(iuncros,102)name,numcros
     ENDIF

111  CONTINUE
! end loop on asteroids
  ENDDO
  STOP ' Warning, too many objects. Reached end of loop.' 
! =========== EOF EXIT POINT ===================================        
999 CONTINUE 

102  FORMAT(a9,2x,i2)

  CALL filclo(iuncan,' ')
  CALL filclo(iuncros,' ')

END PROGRAM placros_NEAs
