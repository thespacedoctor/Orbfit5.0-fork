!     ==================================================================
!     subroutine to get angles (passed in radians) in [0,2pi] interval  
      SUBROUTINE rad02pi(lambda,chl) 
      USE fund_const
      IMPLICIT NONE 
!     INPUT                                                             
      DOUBLE PRECISION lambda 
!     OUTPUT                                                            
      DOUBLE PRECISION chl 
!     ====== sin(lambda), cos(lambda) ===========                       
      DOUBLE PRECISION x,y 
                                                                        
      x=sin(lambda) 
      y=cos(lambda) 
!     chl in [-pi,pi]                                                   
      chl=datan2(x,y) 
                                                                        
!     chl in [0,360]                                                    
      IF((chl.ge.0.d0).and.(chl.le.pig)) THEN 
         chl = chl 
      ELSEIF((chl.lt.0.d0).and.(chl.ge.-pig)) THEN 
         chl = dpig + chl 
      ELSE 
         WRITE(*,*)'ERROR!!!!' 
      ENDIF 
                                                                        
!     control                                                           
      if(chl.gt.dpig) then 
         write(*,*)'ERROR ERROR ERROR !!!' 
      endif 
                                                                        
      RETURN 
      END                                           
