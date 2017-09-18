!     programma di sorting                                              
      SUBROUTINE sortpre(name) 
      IMPLICIT NONE 
!     asteroid name                                                     
      CHARACTER*9  name 
!     name lenght                                                       
      INTEGER lnam 
!     input file data                                                   
      DOUBLE PRECISION time(100),omeg(100),Omnod(100),ecc(100),inc(100) 
      INTEGER nnod(100) 
!     output file data                                                  
      DOUBLE PRECISION ordtime(100),ordomeg(100),ordOmnod(100) 
      DOUBLE PRECISION ordecc(100),ordinc(100) 
      INTEGER ordnnod(100) 
!     loop indexes                                                      
      INTEGER i,j,k,l 
!     max numb of iterations                                            
      INTEGER jmax 
!     =============== END INTERFACE =========================           
                                                                        
!     compute name lenght                                               
      CALL rmsp(name,lnam) 
!     write(*,*)name(1:lnam)//'.pre'                                    
!     input file                                                        
      OPEN(1,file=name(1:lnam)//'.pre',status='old') 
                                                                        
!     initialisation                                                    
      READ(1,*)time(1),omeg(1),Omnod(1),ecc(1),inc(1),nnod(1) 
                                                                        
      ordtime(1)=time(1) 
      ordomeg(1)=omeg(1) 
      ordOmnod(1)=Omnod(1) 
      ordecc(1)=ecc(1) 
      ordinc(1)=inc(1) 
      ordnnod(1)=nnod(1) 
                                                                        
      DO j=2,1000 
         READ(1,*,end=34)time(j),omeg(j),Omnod(j),ecc(j),inc(j),nnod(j) 
                                                                        
         DO i=1,j-1 
            IF(time(j).lt.ordtime(j-i)) THEN 
                                                                        
!     change order                                                      
               ordtime(j-i+1)=ordtime(j-i) 
               ordomeg(j-i+1)=omeg(j-i) 
               ordOmnod(j-i+1)=Omnod(j-i) 
               ordecc(j-i+1)=ecc(j-i) 
               ordinc(j-i+1)=inc(j-i) 
               ordnnod(j-i+1)=nnod(j-i) 
!                                                                       
               ordtime(j-i)=time(j) 
               ordomeg(j-i)=omeg(j) 
               ordOmnod(j-i)=Omnod(j) 
               ordecc(j-i)=ecc(j) 
               ordinc(j-i)=inc(j) 
               ordnnod(j-i)=nnod(j) 
!                                                                       
            ELSEIF(time(j).ge.ordtime(j-i)) THEN 
                                                                        
               ordtime(j-i+1)=time(j) 
               ordomeg(j-i+1)=omeg(j) 
               ordOmnod(j-i+1)=Omnod(j) 
               ordecc(j-i+1)=ecc(j) 
               ordinc(j-i+1)=inc(j) 
               ordnnod(j-i+1)=nnod(j) 
!                                                                       
               GOTO 33 
            ENDIF 
         ENDDO 
   33    CONTINUE 
                                                                        
      ENDDO 
                                                                        
   34 CONTINUE 
      jmax=j-1 
                                                                        
!     output file                                                       
      OPEN(2,file=name(1:lnam)//'.pro',status='unknown') 
                                                                        
!     writing data                                                      
      DO k=1,jmax 
         WRITE(2,100)ordtime(k),ordomeg(k),ordOmnod(k),ordecc(k),       &
     &        ordinc(k),ordnnod(k)                                      
  100    FORMAT(f10.2,1x,f11.5,1x,f11.5,1x,f10.7,1x,                    &
     &        f13.7,i3)                                                 
                                                                        
      ENDDO 
                                                                        
      CLOSE(2) 
      CLOSE(1) 
                                                                        
      RETURN 
    END SUBROUTINE sortpre
