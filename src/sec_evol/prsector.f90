!     =========================================================         
      SUBROUTINE prsector(f,nsec,saf,iz,trans,n) 
      USE fund_const
      IMPLICIT NONE 
!     ========== INPUT ============                                     
!     n = number of iterations; nsec = number of sectors chosen         
      INTEGER n,nsec 
      DOUBLE PRECISION f(n) 
!     saf defines the transition band                                   
      DOUBLE PRECISION saf 
!     ========== OUTPUT ===========                                     
!     trans =.false. (not in transition band)                           
      LOGICAL trans(n) 
!     sector found                                                      
!     (in this case the possible choice is in a subdivision like:       
!     ...,0,1,2,3,0,1,2,3,0,1,2,... because nsec is 4)                  
      INTEGER iz(n) 
!     =========================================================         
!     ssiz = sector size                                                
      DOUBLE PRECISION ssiz,fs,df 
!     loop index                                                        
      INTEGER i  
!     ================== end interface ========================         
                                                                        
      ssiz=dpig/nsec 
      DO 1 i=1,n 
         fs=f(i)+saf 
         IF(fs.gt.0.d0)THEN 
            iz(i)=fs/ssiz 
            df=dmod(fs,ssiz) 
            IF(df.lt.2.d0*saf) THEN 
               trans(i)=.true. 
            ELSE 
               trans(i)=.false. 
            ENDIF 
         ELSE 
            iz(i)=fs/ssiz-1 
            df=dmod(-fs,ssiz) 
            IF(df.gt.ssiz-2*saf) THEN 
               trans(i)=.true. 
            ELSE 
               trans(i)=.false. 
            ENDIF 
         ENDIF 
    1 END DO 
      RETURN 
    END SUBROUTINE prsector
                                                                        
!     ======================================================            
!     =============== function fntarpt ======================           
!     ======================================================            
      DOUBLE PRECISION FUNCTION fntarpt(before,after) 
      IMPLICIT NONE 
!     PREVIOUS SECTOR, PRESENT SECTOR                                   
      INTEGER before,after 
!     end interface                                                     
                                                                        
      IF (after-before.eq.-1.or.after-before.eq.3) THEN 
         fntarpt=before 
      ELSEIF (after-before.eq.1.or.after-before.eq.-3) THEN 
         fntarpt=after 
      ENDIF 
                                                                        
      RETURN 
    END FUNCTION fntarpt
