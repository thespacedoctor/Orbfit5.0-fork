! ******************************************
!  SINGULAR AVERAGED EQUATIONS INTEGRATOR            
! ******************************************
! last modified 11/10/2010 GFG
SUBROUTINE prsingeq(iunit,iunpre,elkep,deltat)
  USE critical_points,ONLY: nminx
  USE fund_const
  USE rkg_local
  USE orbit_elements
  USE planet_orbits
  USE right_hand_side,ONLY: secpert 
  IMPLICIT NONE
! =============================================================
  INTEGER,INTENT(IN) :: iunit,iunpre !output file ids
  TYPE(orbit_elem),INTENT(IN) :: elkep
  REAL(KIND=dkind),INTENT(IN) :: deltat 
! --------- END INTERFACE -------------------------------
  REAL(KIND=dkind) :: tast0 ! asteroid elements initial time
  REAL(KIND=dkind) :: a,e !initial semimajor axis, eccentricity
  INCLUDE 'pldata.h90'
  REAL(KIND=dkind) :: t ! time increment (yr)
  REAL(KIND=dkind) :: y(4),yat(4) ! dynamical variables: 
                                  ! y at tast0+t*365.25
                                  ! yat at tast0 + (t+h)*365.25   
  REAL(KIND=dkind) :: dy(4) ! auxiliary derivatives variable
  REAL(KIND=dkind) :: omfreq,hmax !omega initial frequency, maximum step
  INCLUDE 'neoopt.h90' ! numeric integrator parameters
! epnod (in neoopt.h) is needed for the node closeness check        
  INTEGER nit 
! hmax is the maximum step, h is the current step,
! ck,ck1 are input and output for the polynomial extrapolation
! which will be needed to select a reduced step
  REAL(KIND=dkind) :: h,ck(ismax,4),yi(4,ismax),ck1(ismax,4)
! stopflag is the flag that tells me when to stop
  INTEGER stopfl
! =================================================================
!     extr is a flag that controls EXTRAPOLATION:     
!
!     extr = 1 : extrapolation can be done because we can use the
!     same step of the previous time
!
!     extr = 0 : no extrapolation will be used because we are at first
!     step or we have changed step so that we cannot use previous data
!  
!     extr = -1 : the same coefficients ck (=derivatives of Hbar)
!     are used but the right hand side of equation have been changed 
!     because of the jump
  INTEGER :: extr
! =================================================================
  INTEGER :: issa ! issa=1 asc node cros.; issa=-1 desc. node cross.
! ===================== NODAL DISTANCES ============================
! dn1 is the nodal distance at the first extreme, dn2 at the last one
  REAL(KIND=dkind) :: dn1(nplax2),dn2(nplax2)
! dni= nodal distances at intermediate points; ddni= their derivatives       
  REAL(KIND=dkind) :: dni(nplax2,ismax),ddni(nplax2,ismax)
! dn,ddn are like dni,ddni             
  REAL(KIND=dkind) :: dn(ismax),ddn(ismax)
  REAL(KIND=dkind) :: mutom,mutompl
  INTEGER nnod1,nnod2,nno,nodcod,nc,nnodi(ismax) 
  REAL(KIND=dkind) :: eat,ainc
! ==================================================================
! dt serve per la modifica del passo (regula falsi, estrapolaz.,..) 
! hn e' una variabile ausiliaria, e rappresenta anch'essa il passo  
  REAL(KIND=dkind) :: dt,hn
! ==================== DERIVATIVE JUMP =============================
! These quantities represent the discontinuity jump in the derivatives
! with respect to G and omega
  REAL(KIND=dkind),DIMENSION(nminx) :: ddom,ddg,ddomnod,ddzl
  TYPE(orbit_elem) :: elpl ! planet orbital elements
  INTEGER :: fail_flag
! =============== secpert call =================================
  REAL(KIND=dkind) :: hamil0,error ! output
! =================================================================
  INTEGER :: nst,i,j !loop indexes, units
  REAL(KIND=dkind) :: ddd,eee 
  INTEGER :: nnn 
! ******************************************************************

  dn1(1:nplax2) = 0.d0
  dn2(1:nplax2) = 0.d0

! loading filtered planetary ephemerides from vpla.fil
  CALL read_pla('fil')

! Planet keplerian elements
  DO i = inpl,ioupl
     el_pla(i)=undefined_orbit_elem
     el_pla(i)%t=elkep%t
     el_pla(i)%coo='EQU'
! interpolate for the i-th planet orbit at asteroid initial time tast0
     CALL planet_elems(i,el_pla(i)%t,el_pla(i)%coord(1:5))
     CALL coo_cha(el_pla(i),'KEP',el_pla(i),fail_flag)
     IF(fail_flag.ge.5)THEN
        WRITE(*,*)'error in coo_cha! fail_flag=',fail_flag
     ENDIF
  ENDDO

! initial step computation
  CALL prsingopt(elkep,omfreq,hmax)
  write(*,*)'omfreq = ',omfreq 

! =============== INITIAL CONDITIONS ========================           
! =========== y(1)=omega, y(2)=G, y(3)=Omega y(4)=Z =========           
  a=elkep%coord(1)
  e=elkep%coord(2)
  y(1)=elkep%coord(5)
  y(2)=ky*sqrt(a*(1.d0-e**2)) 
  y(3)=elkep%coord(4)
  y(4)=y(2)*cos(elkep%coord(3))

  tast0=elkep%t ! initial asteroid time (MJD)
  t=0.d0 ! initial time increment (yr)
  h=hmax ! h is the default time step (initially set at the maximum value)
  
  hamil0=0.d0
  !     compute planet elements at new step t+h
  DO i = inpl,ioupl 
     CALL secpert(i,el_pla(i),y(1),y(3),y(2),y(4),a,ddd,eee,nnn)
     hamil0=hamil0+ddd
  ENDDO
  write(*,*)'t,hamil0,zl',t,hamil0,y(4)
  write(23,*)t,hamil0,y(4)

! constant step, first step                                         
  extr=0 
! distance computation and output at first step 
  CALL prdnod(y,a,dn1,nnod1)
  nno=nodcod(nnod1) 
  nit=0 
! ========= writing output on screen =========
  CALL proutele(0,t,y,a,nit,dn1(nno),nnod1)
! ========= writing output on file ============
  CALL proutele(iunit,t,y,a,nit,dn1(nno),nnod1)
!     storage of initial state for stop control                         

  stopfl=0
  CALL prstop(iunpre,tast0,omfreq,t,y,a,stopfl) 
  write(*,*) 'stop_flag = ',stopfl

! *****************************************************************
! *** sostituire PRSTOP con un semplice criterio del tipo
!  IF (t+h.gt.deltat) THEN
!     h = deltat
!     stopfl=1 
!  ENDIF
! *****************************************************************

! INTEGRATION MAIN DO LOOP 
  DO                                                       
     IF (stopfl.lt.0) EXIT ! stop control
3    CONTINUE
!     one step of Runge Kutta Gauss                                     
     CALL rkgstep(extr,tast0,t,y,h,yat,a,isrk,eprk,ck,yi,nit)

     hamil0=0.d0
!     compute planet elements at new step t+h
     DO i = inpl,ioupl 
        el_pla(i) = undefined_orbit_elem
        el_pla(i)%t = tast0 + (t+h)*365.25d0
        el_pla(i)%coo = 'EQU'
!        write(*,*)'tast0+(t+h)*365.25 =',el_pla(i)%t
        CALL planet_elems(i,el_pla(i)%t,el_pla(i)%coord(1:5))
        CALL coo_cha(el_pla(i),'KEP',el_pla(i),fail_flag)
        IF(fail_flag.ge.5)THEN
           WRITE(*,*)'error in coo_cha! fail_flag=',fail_flag
        ENDIF

! y(1)=omega, y(2)=G, y(3)=Omega y(4)=Z 
        WRITE(24,*) 'elpla',el_pla(i)%coord(1:2)
        CALL secpert(i,el_pla(i),yat(1),yat(3),yat(2),yat(4),a,ddd,eee,nnn)
        hamil0=hamil0+ddd
     ENDDO
     write(*,*)'t,hamil0,zl',t+h,hamil0,yat(4)
     write(23,*)t+h,hamil0,yat(4)

!     distance computation and output at new step
     CALL prdnod(yat,a,dn2,nnod2)
     nno=nodcod(nnod2)
     write(*,*) yat,a,dn2(1),nnod2
                                                
!     WRITING OUTPUT ON SCREEN
     CALL proutele(0,t+h,yat,a,nit,dn2(nno),nnod2)

!     distance and its derivatives computation in the intermediate points
     DO j=1,isrk 
        DO i=1,4
           dy(i)=ck(j,i) 
        ENDDO
        CALL prddnod(yi(1:4,j),dy,a,dni(1:2*nplax,j),ddni(1:2*nplax,j),nnodi(j))
     ENDDO
! ======== normal step (before the check) ===========
     hn=h !hn is an auxiliary variable needed to store h

!     ========================================================          
!     CHECK THE POSSIBILITIES THAT CAN HAPPEN AT NODE CROSSING          
!     ========================================================          
     DO nc = 2*inpl-1,2*ioupl 
        
        eat=sqrt(1.d0-yat(2)**2/(ky**2*a)) 
        ainc=acos(yat(4)/(ky*sqrt(a*(1.d0-eat**2)))) 
!     ========                                                          
!     CASE 1                                                            
!     ========                                                          
!     check if node crossing is taking place                            
        IF(abs(dn2(nc)).lt.max(epnod*sin(ainc),5.d-5))THEN 
           WRITE(*,*)'NODE CROSSING',dn2(nc),nc 
           WRITE(*,*)'epnod*sinI=',epnod*sin(ainc) 
                                                                        
!     writing output on file *.pre                                      
           WRITE(iunpre,100)t+h,yat(1)*degrad,yat(3)*degrad,eat,    &
                & ainc*degrad,nnod2                                   
100        FORMAT(f10.2,1x,f11.5,1x,f11.5,1x,f10.7,1x,f13.7,i3) 
                                                          
!     *** modifying array with the JUMP ***
!     compute new RHS coeff as if the jump did not take place
           CALL kintrp(ck,ck1,isrk,4) 
           CALL prsalto(yat,a,el_pla((nc+1)/2),ddg,ddom,ddomnod,ddzl,nc) 
           
           WRITE(*,*)'***************************************' 
           WRITE(*,*)'DERIVATIVE VARIATION AT NODE CROSSING:' 
           WRITE(*,*)'    G derivative jump=',ddg(1)
           WRITE(*,*)'omega derivative jump=',ddom(1)
           WRITE(*,*)'Omega_nod derivative jump=',ddomnod(1)
           WRITE(*,*)'    Z derivative jump=',ddzl(1)
           WRITE(*,*)'***************************************' 
                                                                        
           IF(dn1(nc).gt.0.d0)THEN 
              issa=+1 
           ELSE 
              issa=-1 
           ENDIF
                      
! y(1:4) = (g,G,z,Z)                                                  
           DO j=1,isrk
              ck(j,1)=ck1(j,1)-issa*ddG(1)
              ck(j,2)=ck1(j,2)+issa*ddom(1)
! 
!              ck(j,3)=ck1(j,3)+issa*ddomnod(1)
!              ck(j,4)=ck1(j,4)-issa*ddzl(1)
! changed on 5/7/2016
              ck(j,3)=ck1(j,3)-issa*ddzl(1)
              ck(j,4)=ck1(j,4)+issa*ddomnod(1)
           ENDDO
                                                                        
           WRITE(*,*)'=============================================' 
           WRITE(*,*)'extrapolated derivatives after the JUMP' 
           WRITE(*,*) (ck(1,i),i=1,4) 
           WRITE(*,*)'=============================================' 
                                                                        
!     next step extrapolation will be needed, starting with the jump    
           extr=-1 
                                                                        
!     ========
!     CASE 2       immediately after node crossing                       
!     ========   
        ELSEIF(abs(dn1(nc)).lt.max(epnod*sin(ainc),5.d-5))THEN 
           WRITE(*,*)'AFTER NODE CROSSING',nc 
           WRITE(*,*)'epnod*sinI=',epnod*sin(ainc) 
           WRITE(*,*)'=============================================' 
           WRITE(*,*)'derivatives after the JUMP' 
           WRITE(*,*) (ck(1,i),i=1,4) 
           WRITE(*,*)'=============================================' 
           extr=0 
                                                                        
!     reset step to the original lenght
           hn=hmax

!     ========
!     CASE 3    check if any nodal distance has changed sign
!     ========
        ELSEIF(dn2(nc)*dn1(nc).lt.0.d0)THEN 

!     --- CASE 3.1 -----------------------------------------------------
!     sign change happened between the second last and the last point
           IF(dni(nc,isrk)*dn1(nc).gt.0.d0)THEN 
!     regula falsi to choose next step                                  
              write(*,*)'regula falsi da caso 3.1 a caso 1' 
              dt=h/(dn1(nc)-dn2(nc))*dn2(nc) 
!     HINT! dt e' sempre negativo
              IF(abs(dt).lt.abs(hn))THEN 
                 WRITE(*,*)'controllo: dt prodotto =',dt 
!     modifying step after regula falsi                                 
                 h=h+dt 
                 write(*,*)'la regula falsi produce h=',h 
                 extr=0 
                 GOTO 3 
              ENDIF
!     --- CASE 3.2 -----------------------------------------------------
!     sign change happened in intermediate points so the step is halved 
           ELSE 
!     halving step                                                      
              h=h/2.d0 
              WRITE(*,*)' halving ',h 
              extr=0 
              GOTO 3 
           ENDIF

!     =========
!     CASE 4         no node crossing
!     =========
!     nella situazione piu' comune faccio sempre un controllo se al passo
!     successivo ho dei problemi di node crossing; per fare questo      
!     uso la subroutine prextrdn.f90
        ELSE 
!     write(*,*)'case 4'

!     extrapolation to select next step
           DO j=1,isrk 
              ddn(j)=ddni(nc,j) 
              dn(j)=dni(nc,j) 
           ENDDO
           
!     only if the nodal distance is less than 1AU take
!     into account the extrapolation done with regula falsi
!     and polynomial interpolation with ddn

           IF(abs(dn2(nc)).le.1.d0) THEN 

!     write(*,*)'call prextrdn, nc = ',nc
              CALL prextrdn(dn1(nc),dn2(nc),h,dn,ddn,isrk,dt) 
!     WRITE(*,*)dt
!     dt must have the same sign as h and be small
              IF(dt*h.gt.0.d0.and.abs(dt).lt.abs(hn))THEN 
!     only if the nodal distance is less than 1AU take
!     into account the extrapolation done with
                 IF(abs(dn2(nc)).le.1.d0) THEN 
!     WRITE(*,*)'dt after extrapolation=',dt
!     ====================================
                    hn=dt
!     ====================================
                    extr=0 
                 ENDIF
                                                                        
              ELSE 
!     no node crossing in sight, use constant step
                 extr=1 
              ENDIF
           ENDIF                           
!     end main if cycle
        ENDIF
!     end do loop on nodal distances with all planets
     ENDDO

!     ============== ACCEPTED STEP ================
!     decide if I have to modify stopflag:
!     I notice that output file is open, so
!     in prstop I can say to write in it

!     control
!     write(*,*)'prima di prstop:omeg,omnod',yat(1),yat(3)
     CALL prstop(iunpre,tast0,omfreq,t+h,yat,a,stopfl) 
!     IF(stopfl.lt.0)THEN
!        write(*,*) stopfl
!     ENDIF
!     ======== WRITING OUTPUT ON FILE ========
     CALL proutele(iunit,t+h,yat,a,nit,dn2(nno),nnod2)

!     accepted step
     DO i=1,4 
        y(i)=yat(i) 
     ENDDO
     t=t+h 
     h=hn 
     DO nc = inpl,2*ioupl-inpl+1 
        dn1(nc)=dn2(nc) 
     ENDDO
     nnod2=nnod1 
                                                                        
!     ==== END MAIN DO LOOP ====
  ENDDO
  RETURN 
END SUBROUTINE prsingeq
