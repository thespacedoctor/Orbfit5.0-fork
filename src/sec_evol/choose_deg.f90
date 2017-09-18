!     ===========================================================       
!     subroutine to obtain angles in [-180,180]                         
!     lambda is given in degrees                                        
      SUBROUTINE choose_deg(lambda,chl) 
      USE fund_const
      IMPLICIT NONE 
!     INPUT                                                             
      DOUBLE PRECISION lambda 
!     OUTPUT                                                            
      DOUBLE PRECISION chl 
!     ====== sin(lambda), cos(lambda) ===========                       
      DOUBLE PRECISION x,y 
                                                                        
      x=sin(radeg*lambda) 
      y=cos(radeg*lambda) 
      chl=degrad*atan2(x,y) 
                                                                        
      RETURN 
    END SUBROUTINE choose_deg
