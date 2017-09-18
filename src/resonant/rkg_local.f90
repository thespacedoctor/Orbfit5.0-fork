MODULE rkg_local
  USE fund_const
  USE orbit_elements
  USE planet_orbits
  IMPLICIT NONE
  PRIVATE

! =================================================================
!  numerical integrators  parameters                                    
  INTEGER, parameter, public :: ismax=20, itmax=20 
!  coefficients RKG for rkimp
  REAL(KIND=dkind) :: a(ismax,ismax),b(ismax),c(ismax)
!  coefficients interpolator for k array, for kintrp
  REAL(KIND=dkind) :: a1(ismax,ismax)
!  LOGICAL,PARAMETER :: check_der = .true.
  LOGICAL,PARAMETER :: check_der = .false.

  PUBLIC rkgstep,kintrp,check_der

CONTAINS
! =================================================================      
!  SUBROUTINE rkgstep(extr,tast0,t,y,h,yat,aa,isrk,eprk,ck,yi,nit)
  SUBROUTINE rkgstep(extr,tast0,t,y,h,yat,isrk,eprk,ck,yi,nit)
! ===================== INTERFACE =================================
    INTEGER,INTENT(IN) :: extr !   extrapolation parameters 
    REAL(KIND=dkind),INTENT(IN) :: tast0,t ! time (initial, final)
! -----------------------------------------------------------------
    REAL(KIND=dkind),INTENT(IN) :: y(6) ! time, dynamical variables
! y(1)=omega, y(2)=G, y(3)=Omega y(4)=Z y(5)=S, y(6)=sigma
! note that the order of action-angle is not kept
! -----------------------------------------------------------------
    REAL(KIND=dkind),INTENT(IN) :: h! h is the current step
!    REAL(KIND=dkind),INTENT(IN) :: aa ! semimajor axis
    REAL(KIND=dkind),INTENT(OUT) :: yat(6) !output dynamical variables
! parameters for integrator
    INTEGER,INTENT(IN) :: isrk
    REAL(KIND=dkind),INTENT(IN) :: eprk
    REAL(KIND=dkind),INTENT(INOUT) :: ck(ismax,6) !values of RHS of eqs
! intermediate values of dynamical variables
    REAL(KIND=dkind),INTENT(OUT) :: yi(6,ismax) 
    INTEGER,INTENT(INOUT) :: nit
! --------- end interface ----------------------------------------
    INTEGER :: nstep,isfl,lit
    REAL(KIND=dkind) :: hmax! hmax is the maximum step
    REAL(KIND=dkind) :: ck1(ismax,6) 
! ck,ck1 sono l'input e l'output per l'estrapolazione polinomiale
! che mi servira' per la scelta di un passo ridotto
    INTEGER :: i,j ! loop indexes
    INTEGER :: lflag ! flag 
    DATA lflag/0/ 
    SAVE 
! ================================================================== 
! ****************************************************************** 
! Options for the integrator, only for the first time               
    IF(lflag.eq.0)THEN 
       CALL prlegnum(isrk,isfl) 
!       write(*,*) a1(1:4,1:4)
       IF(isfl.ne.0)THEN 
          WRITE(*,*)' coefficienti non trovati per is=', isrk 
          STOP 
       ENDIF
       lit=8 
       nit=0 
       lflag=1 
    ENDIF
! EXTRAPOLATION FROM THE PREVIOUS STEP                              
! ==================================================================        
    IF(extr.eq.1)THEN 
! ================= fixed lenght step ==============================        
       CALL kintrp(ck,ck1,isrk,6)
       DO j=1,isrk 
          DO i=1,6
             ck(j,i)=ck1(j,i) 
          ENDDO
       ENDDO
    ELSEIF(extr.eq.0)THEN 
! variable lenght step (also first step)                            
       DO j=1,isrk 
          DO i=1,6 
             ck(j,i)=0.d0 
          ENDDO
       ENDDO
    ELSEIF(extr.eq.-1)THEN 
! do nothing; use ck as precomputed                                 
    ENDIF
!    write(*,*)ck
! ******************************************************************
! RKG step (RKG is the Runge Kutta Gauss subroutine                 
! that computes the solution of Hamilton's equations                
! ****************************************************************** 
!    write(*,*)t1,t,y,h,a
    CALL prrkg(tast0,t,y,h,yat,isrk,lit,yi,ck,eprk,nit) 

 END SUBROUTINE rkgstep
! ******************************************************************       
!  {\bf  kintrp} ORB8V                                                  
!                                                                       
!  scopo : predizione dei valori di ck (= f(z_j)) per interpolazione da quelli
!          del passo precedente (con polinomio di grado is-1)           
!          serve per avere condizioni iniziali vicine al punto unito    
!          nel procedimento iterativo per risolvere le equazioni        
!          implicite la cui soluzione e' la matrice ck (= z_j).        
!  input :                                                              
!       ck1(ismax,nvar) : ck al passo precedente                        
!       is : numero di stadi del metodo rk                              
!       nvar : numero di variabili nell'equazione di moto               
!  output:                                                              
!       ck(ismax,nvar) : valori interpolati                             
!  osservazione : non va usato nel passo iniziale e se rispetto al      
!                 passo precedente e' cambiato il passo, oppure is.     
!****************                                                       
!   static memory not required                                          
!****************                                                       
  SUBROUTINE kintrp(ck,ck1,is,nvar)
! ===================== INTERFACE =================================
    INTEGER,INTENT(IN) :: nvar ! dimensioni dipendenti da ismax (qui=ismax)
    REAL(KIND=dkind),INTENT(IN) :: ck(ismax,nvar)
    REAL(KIND=dkind),INTENT(OUT) :: ck1(ismax,nvar) 
    INTEGER,INTENT(IN) :: is ! number of intermediate steps
! ======================= END INTERFACE ============================ 
    INTEGER :: j,n,jj ! loops indexes
    REAL(KIND=dkind) :: de ! temporary variables
! ================================================================== 
! ******************************************************************
    DO j=1,is 
       DO  n=1,nvar 
          de=0.d0 
          DO jj=1,is 
             de=de+a1(j,jj)*ck(jj,n)
          ENDDO
          ck1(j,n)=de 
       ENDDO
    ENDDO
  END SUBROUTINE kintrp

! ================================================================== 
!     Runge-Kutta-Gauss, used as a pure single step
! ================================================================== 
!  SUBROUTINE prrkg(tast0,t,y,h,yat,aa,isrk,lit,yi,ck,eprk,nit)
  SUBROUTINE prrkg(tast0,t,y,h,yat,isrk,lit,yi,ck,eprk,nit)
! ===================== INTERFACE ==================================
! INPUT
    REAL(KIND=dkind),INTENT(IN) :: tast0 !tast0 = initial time of asteroid (MJD)
    REAL(KIND=dkind),INTENT(IN) :: t !t is the time for y vector (yr)
    REAL(KIND=dkind),INTENT(IN) :: h !h is the step
!    REAL(KIND=dkind),INTENT(IN) :: aa !asteroid semiaxis
    REAL(KIND=dkind),INTENT(IN) :: y(6) !y(1)=omega,y(2)=G,y(3)=Omega,y(4)=zl,y(5)=S,y(6)=sigma
! OUTPUT
    REAL(KIND=dkind),INTENT(OUT) :: yat(6) ! yat is y vector at time t+h
    INTEGER,INTENT(IN) :: isrk,lit
    INTEGER,INTENT(INOUT) :: nit ! nit is the number of iterations
    REAL(KIND=dkind),INTENT(INOUT) :: yi(6,ismax),ck(ismax,6)
! ======================= END INTERFACE ============================
    REAL(KIND=dkind) :: dery(6),de
    REAL(KIND=dkind) :: eprk 
! controls etc.  
!    REAL(KIND=dkind) :: ep(itmax),t1(ismax),t_int(ismax) 
    REAL(KIND=dkind) :: ep(itmax),t_int(ismax) 
    INTEGER :: i,j,it,jj,n 
! per i loop(forse non servono)
    INTEGER :: s1,s2
    INTEGER :: fail_flag
! ===================== right hand side ============================
    REAL(KIND=dkind) :: sigma,bigelle,om,g,omnod,zl,ddd(6),eee(6) 
    INTEGER :: nnn(6)
! ==================================================================
! dn,dnod                                                           
    REAL(KIND=dkind) :: dn(ismax) 
    INTEGER :: nnod
! ==================================================================
    INCLUDE 'pldata.h90'
! ******************************************************************
!     intermediate times between t and (t+h)
    DO j=1,isrk                                                       
       t_int(j)=t + h*c(j) 
!       write(*,*)'j,t_int',j,t_int(j)
    ENDDO
!     azzeramento                                                       
    DO  it=1,itmax 
       ep(it)=0.d0 
    ENDDO
! Gauss-Seidel for ck's                                             
! begin iterations for ck's                                         
    DO 1 it=1,itmax 
!    main loop                                                         
       DO 11 j=1,isrk 
          DO i=1,6
             de=0.d0 
             DO jj=1,isrk 
                de=de+a(j,jj)*ck(jj,i) 
             ENDDO
             yi(i,j)=de*h+y(i) 
          ENDDO
! right hand side of Hamilton's equations                           
          om=yi(1,j) 
          g=yi(2,j) 
          omnod=yi(3,j)
          zl=yi(4,j)
          bigelle=yi(5,j)
          sigma=yi(6,j)

          
          DO n=inpl,ioupl
!             el_pla(n)=undefined_orbit_elem
!             el_pla(n)%t=t_int(j)
             el_pla(n)%coo='EQU'
             CALL planet_elems(n,tast0+t_int(j)*365.25d0, &
                  & el_pla(n)%coord(1:6))
             CALL coo_cha(el_pla(n),'KEP',el_pla(n),fail_flag)
             IF(fail_flag.ge.5)THEN
                WRITE(*,*)'error in coo_cha! fail_flag=',fail_flag
             ENDIF
          ENDDO
!          write(*,*)'!!!!y prima di chooserhs=',y
     !     write(*,*)'!!!!om,omnod,g,zl,bigelle,sigma',om,omnod,g,zl,bigelle,sigma
          CALL chooserhs(y,om,omnod,g,zl,bigelle,sigma,ddd,eee,nnn)
! note that in the correspondign program in sec_evol/ we changed some signs 
! in the following raws because we gave a different interpretation to 
! the right hand sides

! y(1)=omega, y(2)=G, y(3)=Omega y(4)=Z y(5)=S, y(6)=sigma
          dery(1)=ddd(1) !dH/dG
          dery(2)=ddd(2) !-dH/dg
          dery(3)=ddd(3) !dH/dZ
          dery(4)=ddd(4) !-dH/dz
          dery(5)=ddd(5) !-dH/dsigma
          dery(6)=ddd(6) !dH/dS

!          IF(check_der)THEN
!             write(*,*)'dH/d(G,g,Z,z,sigma,S)'
!             write(*,*) -dery(2),dery(1),-dery(4),dery(3),dery(6),-dery(5)
!             STOP
!          ENDIF

! convergence control                                               
          DO i=1,6
             ep(it)=ep(it)+dabs(dery(i)-ck(j,i)) 
          ENDDO
! updating intermediate values                                      
          DO i=1,6 
             ck(j,i)=dery(i) 
          ENDDO
11     ENDDO
! check if iterazions g-s have been done                            
       ep(it)=ep(it)/isrk
       IF(ep(it).lt.eprk)THEN 
          GOTO 2 
       ELSE 
          IF(it.ge.lit)THEN 
! too many iterations in Gauss-Seidel                               
             WRITE(*,*)' rkg: not convergent in ',it,' iterations' 
             WRITE(*,*)t,(ep(jj),jj=1,it) 
          ENDIF
       ENDIF
1   ENDDO
! too many iterations                                               
    WRITE(*,*)' RKG: NOT CONVERGENT' 
2   CONTINUE 
! ============= computation of new points ==========================
    DO i=1,6
       de=0.d0 
       DO  j=1,isrk 
          de=de+b(j)*ck(j,i) 
       ENDDO
       yat(i)=y(i)+h*de 
    ENDDO
    nit=it-1
  END SUBROUTINE prrkg
! ==================================================================
! ******************************************************************    
!  LEGNUM                                                               
! ******************************************************************       
! reads  Runge--Kutta--gauss coefficients                               
! to be read in file ./lib/rk.coef                                      
!  is=required number of substeps; isfl=0 if found, otherwise           
!  uses closest available                                               
! ==============INTERFACE===========================================
  subroutine prlegnum(is,isfl)
! ===================== INTERFACE ==================================
! INPUT                                            
    INTEGER :: is 
! OUTPUT                                            
    INTEGER :: isfl
! ======================= END INTERFACE ============================
! input unit, current is                                                
    INTEGER :: iun,is1
! loop indexes                                                          
    INTEGER :: i,j,jj 
! skip trick                                                            
    CHARACTER*1 :: cc(ismax)
! ================================================================== 
! ****************************************************************** 
!    static memory not required 
! ****************************************************************** 
! reads RKG coefficients rk                                             
    isfl=-1 
    iun=20 
    open(iun,file='rk.coe',status='old') 
198 read(iun,100,end=199) is1 
100 format(6x,i4) 
!  control on is compatible wtih parameter ismax                        
    if(is1.gt.ismax.or.is1.le.0)goto 199 
    if(is1.eq.is)then 
       read(iun,101)(c(j),j=1,is1) 
101    format(7x,5d24.16) 
       read(iun,102)(b(j),j=1,is1) 
102    format(7x,5d24.16) 
       DO j=1,is1 
          read(iun,104)(a(i,j),i=1,is1) 
       ENDDO
104    format(3x,5d24.16) 
       DO j=1,is1 
          read(iun,106)(a1(i,j),i=1,is1) 
106       format(4x,5d24.16) 
       ENDDO
       isfl=0 
       goto 199 
    else 
       read(iun,201)(cc(j),j=1,is1) 
201    format(7x,5(23x,a1)) 
       read(iun,201)(cc(j),j=1,is1) 
       DO j=1,is1 
          read(iun,204)(cc(jj),jj=1,is1) 
204       format(3x,5(23x,a1)) 
       ENDDO 
       DO j=1,is1 
          read(iun,206)(cc(jj),jj=1,is1) 
206       format(4x,5(23x,a1)) 
       ENDDO 
       isfl=is1 
       goto 198 
    ENDIF
! end read                                                              
199 close(iun) 

  END subroutine prlegnum
  
END MODULE rkg_local
