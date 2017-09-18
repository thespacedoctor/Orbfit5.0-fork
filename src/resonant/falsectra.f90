! ************************************************************      
! falsectra regula falsi for sector transition                      
! ************************************************************      
SUBROUTINE falsectra(tast0,t1,y1,t2,y2,omfreq,targ,yt,dt) 
  USE fund_const
  USE rkg_local
  IMPLICIT NONE 
! =============================================================     
! INPUT VARIABLES                                                   
! two states with om-target opposite sign                           
  REAL(KIND=dkind) :: tast0,t1,y1(6),t2,y2(6)
  REAL(KIND=dkind) :: omfreq ! frequency for stepsize control
!  REAL(KIND=dkind) :: a ! semimajor axis
  REAL(KIND=dkind) :: targ ! target value of omega 
! =============================================================     
! OUTPUT VARIABLES                                                  
! state at target and time difference with respect to t2            
  REAL(KIND=dkind) :: yt(6),dt 
! =============================================================     
! numeric integrator parameters                                     
  INCLUDE 'neoopt.h90' 
  INTEGER :: nit,extr 
! dt is the current step                                            
  REAL(KIND=dkind) :: ck(ismax,6),yi(6,ismax) 
! ================== local variables ==========================
! regula falsi                                                      
  INTEGER :: it,i 
  REAL(KIND=dkind) :: dom,dom1,dom2,tcur,epsom 
! =============================================================
! *************************************************************
!     write(*,*)t1,y1,t2,y2,omfreq,a,zl,targ,yt,dt                      
!     write(*,*)'target=',targ                                          
!     initialisation                                                    
  epsom=1.d-4 
                                                                        
! ==================== MAIN DO LOOP ===========================
  DO 1 it=1,itmax 
                                                                        
! control on number of regula falsi done                            
     WRITE(*,*)it,' regula falsi' 
                                                                        
! new time by regula falsi                                          
     dt=(t1-t2)*(targ-y2(1))/(y1(1)-y2(1)) 
     write(*,*)'iteration =',it 
     write(*,*)'ho scelto dt =',dt 
     write(*,*)'t1=',t1,' t2=',t2 
                                                                        
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! =============================================================
!     CHECK: t2+dt must be between t1 and t2
! =============================================================
     IF(t2.gt.t1)THEN 
        IF(t2+dt.lt.t1.or.t2+dt.gt.t2)THEN 
           WRITE(*,*)' ahi! t2+dt must be between t1 and t2 ' 
           WRITE(*,*)'t1,t2+dt,t2,dt',t1,t2+dt,t2,dt 
        ENDIF
                                                                        
! (this should not happen)                                          
     ELSEIF(t1.gt.t2)THEN 
        IF(t2+dt.gt.t1.or.t2+dt.lt.t2)THEN 
           WRITE(*,*)' ahi!', t1,t2,dt,t2+dt 
        ENDIF
        
     ELSE 
        write(*,*)' t1=t2 ',t1,t2 
        STOP 
     ENDIF
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
                                                                        
! =============================================================
!                    W A R N I N G   M E S S A G E                          
! =============================================================
! dt should be smaller than pig/(100*omfreq)*1.1                    
     IF(abs(dt).gt.pig/(100*abs(omfreq))*1.1d0)THEN 
        WRITE(*,*)'WARNING: step too long' 
        WRITE(*,*)'dt = ',dt,' GREATER THAN',pig/(100*abs(omfreq)) 
     ENDIF
! =============================================================
                                                                        
! one step with stepsize dt                                         
! no extrapolation (variable step)                                  
     extr=0 
                                                                        
     CALL rkgstep(extr,tast0,t2,y2,dt,yt,isrk,eprk,ck,yi,nit) 
                                                                        
! decide next step                                                  
     dom=yt(1)-targ 
                                                                        
     IF(abs(dom).lt.epsom)THEN 
        RETURN 
     ENDIF
                                                                        
     dom1=y1(1)-targ 
     dom2=y2(1)-targ 
     tcur=t2+dt 
                                                                        
     IF(dom*dom1.lt.0.d0)THEN 
        DO i =1,6
           y2(i)=yt(i) 
        ENDDO
        t2=tcur 
     ELSEIF(dom*dom2.lt.0.d0)THEN 
        DO i =1,6 
           y1(i)=yt(i) 
        ENDDO
        t1=tcur 
     ENDIF
                                                                        
1 ENDDO
! =============================================================   
! if more than itmax iterations have been done                      
  WRITE(*,*)' too many iterations ' 
! =============================================================
  RETURN 
END SUBROUTINE falsectra
