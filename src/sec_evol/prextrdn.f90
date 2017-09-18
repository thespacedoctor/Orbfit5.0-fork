!     subroutine that check if the nodal distance will change           
!     sign if next iteration is done using step h, and decides          
!     which is the more suitable step dt.                               
      SUBROUTINE prextrdn(dn1,dn2,h,dn,ddn,isrk,dt) 
      USE rkg_local
      IMPLICIT NONE 
      INTEGER isrk
      DOUBLE PRECISION,INTENT(IN):: dn1,dn2,h,dn(ismax),ddn(ismax)
      DOUBLE PRECISION,INTENT(OUT):: dt
!     ========= END INTERFACE ============                              
      DOUBLE PRECISION:: dt1 
      DOUBLE PRECISION ddn1(ismax) 
                                                                        
!     regula falsi to select next step                                  
      dt=h/(dn1-dn2)*dn2 
                                                                        
!     control                                                           
!      write(*,*)'*********************************'                    
!      write(*,*)'regula falsi to select next step'                     
!      write(*,*)'dt =',dt                                              
!      write(*,*)'*********************************'                    
                                                                        
!     polynomial extrapolation with the values of the derivatives       
!     of the nodal distance                                             
      CALL kintrp(ddn,ddn1,isrk,1) 
      dt1=-dn2/ddn1(1) 
                                                                        
!      control                                                          
!      write(*,*)'++++++++++++++++++++++++++++++++'                     
!      write(*,*)'extrapolation with ddn gives dt=',dt1                 
!      write(*,*)'++++++++++++++++++++++++++++++++'                     
                                                                        
!     CHECK WHICH IS THE MINIMUM VALUE BETWEEN BOTH CHOICE OF dt        
      IF (abs(dt).le.abs(dt1)) THEN 
         dt = dt 
      ELSEIF (abs(dt).ge.abs(dt1)) THEN 
         dt = dt1 
      ENDIF 
!     control                                                           
!      write(*,*)'choice for dt in prextrdn =',dt                       
                                                                        
    END SUBROUTINE prextrdn
