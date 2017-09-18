!     ===========================================================       
!     subroutine to obtain angles in [0,360]                            
SUBROUTINE choosedeg0_2pi(lambda,chl) 
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
!     chl in [-180,180]                                                 
  chl=degrad*atan2(x,y) 
!     chl in [0,360]                                                    
  IF((chl.ge.0).and.(chl.le.180)) THEN 
     chl = chl 
  ELSEIF((chl.lt.0).and.(chl.ge.-180)) THEN 
     chl = 360 + chl 
  ELSE 
     WRITE(*,*)'ERROR!!!!' 
  ENDIF                                                                      
!     check                                                             
  if(chl.gt.360) then 
     write(*,*)'non va bene per niente!!!' 
  endif
                                                                        
END SUBROUTINE choosedeg0_2pi

!     ===========================================================       
!     subroutine to obtain angles in [-180,360]                         
SUBROUTINE choosedegmpi_pi(lambda,chl) 
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
!     chl in [-180,180]                                                 
  chl=degrad*atan2(x,y) 
!     control                                                           
  if((chl.lt.-180).or.(chl.gt.180)) then 
     write(*,*)'ERROR IN CHOOSEDEG-PI_PI!!!' 
  endif
END SUBROUTINE choosedegmpi_pi

! =========================================================================    
! subroutine to obtain adfl(s); that is the ascending/descending node flag.
! This subroutine is required only in the circulation case!!!       
! =========================================================================
SUBROUTINE chooseadfl(circflag,oom,adfl,ascdiscfl) 
  IMPLICIT NONE 
  INTEGER,INTENT(IN) :: circflag,adfl 
  DOUBLE PRECISION, INTENT(IN) :: oom 
  INTEGER,INTENT(OUT) :: ascdiscfl(4)
! end interface 
  DOUBLE PRECISION :: chom 

!     put oom in [-180,180]                                             
  CALL choosedegmpi_pi(oom,chom) 
 
!  ascdicsfl(1) = 1*adfl
                                                                       
  IF (circflag.eq.1) THEN 
     IF(chom.ge.0.and.chom.le.90) THEN 
        ascdiscfl(2)=-1*adfl 
        ascdiscfl(3)=-1*adfl 
        ascdiscfl(4)=1*adfl 
        
     ELSEIF(chom.gt.90.and.chom.le.180) THEN 
        ascdiscfl(2)=1*adfl 
        ascdiscfl(3)=-1*adfl 
        ascdiscfl(4)=-1*adfl 
        
     ELSEIF(chom.ge.-180.and.chom.lt.-90) THEN 
        ascdiscfl(2)=-1*adfl 
        ascdiscfl(3)=-1*adfl 
        ascdiscfl(4)=1*adfl 
        
     ELSEIF(chom.ge.-90.and.chom.lt.0) THEN 
         ascdiscfl(2)=1*adfl 
        ascdiscfl(3)=-1*adfl 
        ascdiscfl(4)=-1*adfl 
        
     ENDIF
                                                                        
  ELSEIF(circflag.eq.-1) THEN 
     IF(chom.ge.0.and.chom.le.90) THEN 
        ascdiscfl(2)=1*adfl 
        ascdiscfl(3)=-1*adfl 
        ascdiscfl(4)=-1*adfl 
        
     ELSEIF(chom.gt.90.and.chom.le.180) THEN 
        ascdiscfl(2)=-1*adfl 
        ascdiscfl(3)=-1*adfl 
        ascdiscfl(4)=1*adfl 
        
     ELSEIF(chom.ge.-180.and.chom.lt.-90) THEN 
        ascdiscfl(2)=1*adfl 
        ascdiscfl(3)=-1*adfl 
        ascdiscfl(4)=-1*adfl 
                                                                        
     ELSEIF(chom.ge.-90.and.chom.lt.0) THEN 
        ascdiscfl(2)=-1*adfl 
        ascdiscfl(3)=-1*adfl 
        ascdiscfl(4)=1*adfl 
                                                                        
     ENDIF
                                                                        
  ENDIF
                                                                        
!     check                                                             
  WRITE(*,*)'IN CHOOSEADFL:',ascdiscfl 
                                                                        
END SUBROUTINE chooseadfl
      
!     ============================================                      
SUBROUTINE choosephi(oom,adfl,phi,chphi) 
  USE fund_const
  IMPLICIT NONE 
!     ======= input  (omega,ascdiscflag, 0<phi<pi/2 ) ========          
!     ======= adfl = 1 (ascending node cros.) ================          
!     ======= adfl = -1 (descending node cros.) ==============          
  DOUBLE PRECISION,INTENT(IN) :: oom,phi 
  INTEGER,INTENT(IN) :: adfl 
!     output ( value of phi chosen)                                     
  DOUBLE PRECISION,INTENT(OUT) :: chphi
! ======= end interface ======================================
 !     ======= true anomaly ===========                                  
  DOUBLE PRECISION :: effe 
!     ========= CASES for phi ===========                               
!     case = 1         0 <= phi < pi/2                                   
!     case = 2      pi/2 <= phi < pi                                     
!     case = 3    3/2*pi <= phi < 2*pi                                   
!     case = 4        pi <= phi < 3/2*pi                                 
  INTEGER case 
!     ========= end interface =====================                     
                                                                        
!     ========== setting effe =====================                     
  IF(adfl.eq.1) THEN 
!     IN RADIANS                                                        
     effe=-radeg*oom 
  ELSEIF(adfl.eq.-1) THEN 
     effe=pig-radeg*oom 
  ENDIF
!     control                                                           
!      write(*,*)'effe=',degrad*effe                                    
                                                                        
!     ============= setting case ==================                     
!     CASE1 (ASCENDING NODE  POST-PERIHELION)   Ux>=0  Uz>0              
  IF(adfl.eq.1.and.sin(effe).ge.0) THEN 
     case=1 
!     CASE2 (DESCENDING NODE  POST-PERIHELION)  Ux>0  Uz<0              
  ELSEIF(adfl.eq.-1.and.sin(effe).gt.0) THEN 
     case=2 
!     CASE3 (ASCENDING NODE  PRE-PERIHELION)    Ux<0  Uz>0              
  ELSEIF(adfl.eq.1.and.sin(effe).lt.0) THEN 
     case=3 
!     CASE4 (DESCENDING NODE  PRE-PERIHELION)   Ux<=0  Uz<0              
  ELSEIF(adfl.eq.-1.and.sin(effe).le.0) THEN 
     case=4 
  ENDIF
!     ============= choosing sector for phi ====================        
  IF(case.eq.1) THEN 
     chphi=phi 
  ELSEIF(case.eq.2) THEN 
     chphi=180-phi 
  ELSEIF(case.eq.3) THEN 
     chphi=360-phi 
  ELSEIF(case.eq.4) THEN 
     chphi=180+phi 
  ENDIF
!     control                                                           
!      write(*,*)'case=',case,'  chphi=',chphi                          
                                                                        
END SUBROUTINE choosephi
      
