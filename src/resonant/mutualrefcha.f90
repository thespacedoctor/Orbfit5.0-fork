
!     ******************************************************************
!     ***************** CHANGING REFERENCE SYSTEM **********************
!     ******************************************************************
!     ********** written by GIOVANNI F. GRONCHI (2001) *****************
!     ******* Department of Mathematics, UNIVERSITY of PISA ************
!     last modified May 2003 (GFG)
!     ==================================================================
      SUBROUTINE mutualrefcha(ekpl,elkep,mutI,mutom,mutompl)
      USE fund_const
      IMPLICIT NONE 
!     HINT!: all angles are given in radians                           
      DOUBLE PRECISION,INTENT(IN) :: elkep(6)! asteroid keplerian elements 
      DOUBLE PRECISION,INTENT(IN) :: ekpl(6)! Earth keplerian elements  
      DOUBLE PRECISION,INTENT(OUT) :: mutI,mutom,mutompl! mutual variables 
!     =========== end interface ========================================
      DOUBLE PRECISION :: a,e,Inc,om,Omnod 
      DOUBLE PRECISION :: apl,epl,Ipl,ompl,Omnodpl 
!     asteroid and planet coordinates for u = upl = 0                   
      DOUBLE PRECISION :: x0(3),xpl0(3)
!     asteroid and planet coordinates for u = upl = pi/2           
      DOUBLE PRECISION :: xpm(3),xplpm(3) 
      DOUBLE PRECISION :: Ne(3),Nast(3) ! normal vectors
      DOUBLE PRECISION :: ASCMN(3)! ascending mutual node vector 
      DOUBLE PRECISION :: Ivec(3),sinI,cosI 
      DOUBLE PRECISION :: omvec(3),sinom,cosom 
      DOUBLE PRECISION :: omplvec(3),sinompl,cosompl 
      DOUBLE PRECISION :: test ! test variable
      DOUBLE PRECISION :: prscal ! included function  
      EXTERNAL prscal 
!     auxiliary
      DOUBLE PRECISION :: x0tmp(3),xpl0tmp(3),xpmtmp(3),xplpmtmp(3) 
      DOUBLE PRECISION :: Netmp(3),Nasttmp(3),ASCMNtmp(3) 
!     ==================================================================


      apl = ekpl(1) 
      epl = ekpl(2) 
      Ipl = ekpl(3) 
      Omnodpl = ekpl(4) 
      ompl = ekpl(5) 
                                                                        
      a = elkep(1) 
      e = elkep(2) 
      Inc = elkep(3) 
      Omnod = elkep(4) 
      om = elkep(5) 
                                                                        
      CALL compxyz(a,e,Inc,om,Omnod,0.d0,x0tmp) 
      CALL normalize(x0tmp,x0) 
                                                                        
      CALL compxyz(apl,epl,Ipl,ompl,Omnodpl,0.d0,xpl0tmp) 
      CALL normalize(xpl0tmp,xpl0) 
                                                                        
      CALL compxyz(a,e,Inc,om,Omnod,pig/2,xpmtmp) 
      CALL normalize(xpmtmp,xpm) 
                                                                        
      CALL compxyz(apl,epl,Ipl,ompl,Omnodpl,pig/2,xplpmtmp) 
      CALL normalize(xplpmtmp,xplpm) 
                                                                        
!     ==================================================================
!     COMPUTE NORMAL VECTOR TO EARTH ORBIT: Ne                          
      CALL prvec(xpl0,xplpm,Netmp) 
      CALL normalize(Netmp,Ne) 
!                                                                       
!     COMPUTE NORMAL VECTOR TO ASTEROID ORBIT: Nast                     
      CALL prvec(x0,xpm,Nasttmp) 
      CALL normalize(Nasttmp,Nast) 
!     ==================================================================
                                                                        
!     ==================================================================
!     COMPUTE ASCENDING MUTUAL NODE VECTOR ASCMN                        
      cosI =  prscal(Ne,Nast) 
!     check that Ne and Nast are not parallel                           
      IF ((abs(abs(cosI)-1.d0)).lt.1.d-10) THEN 
         WRITE(*,*)'MUTUAL INCLINATION CLOSE TO ZERO!' 
         mutI = 0.d0 
         mutom = om + Omnod 
         mutompl = ompl + Omnodpl 
!     skip next steps                                                   
         GOTO 18 
      ELSE 
         CALL prvec(Ne,Nast,Ivec) 
         CALL normalize(Ivec,ASCMN) 
      ENDIF 
!     ==================================================================
                                                                        
!     ==================================================================
!     COMPUTE ANGLE BETWEEN N1,N2                                       
      CALL lenght(Ivec,sinI) 
      CALL ptriple(Ne,Nast,ASCMN,test) 
      IF (test.gt.0.d0) THEN 
         sinI = sinI 
      ELSEIF (test.lt.0.d0) THEN 
         sinI = -sinI 
      ELSE 
         WRITE(*,*)' !!!!!!!!!!!!!!!!!' 
         WRITE(*,*)' ERROR!!! ERROR!!!' 
         WRITE(*,*)'sinI = ',sinI 
         WRITE(*,*)' !!!!!!!!!!!!!!!!!' 
      ENDIF 
      mutI = datan2(sinI,cosI) 
!     ==================================================================
                                                                        
!     ==================================================================
!     COMPUTE ANGLE BETWEEN ASCMN,xpl0                                  
      cosompl =  prscal(ASCMN,xpl0) 
      CALL prvec(ASCMN,xpl0,omplvec) 
      CALL lenght(omplvec,sinompl) 
      CALL ptriple(ASCMN,xpl0,Ne,test) 
      IF (test.gt.0.d0) THEN 
         sinompl = sinompl 
      ELSEIF (test.lt.0.d0) THEN 
         sinompl = -sinompl 
      ELSE 
         WRITE(*,*)' !!!!!!!!!!!!!!!!!' 
         WRITE(*,*)' ERROR!!! ERROR!!!' 
         WRITE(*,*)'sinompl = ',sinompl 
         WRITE(*,*)' !!!!!!!!!!!!!!!!!' 
      ENDIF 
      mutompl = datan2(sinompl,cosompl) 
!     ==================================================================
                                                                        
!     ==================================================================
!     COMPUTE ANGLE BETWEEN ASCMN,x0                                    
      cosom =  prscal(ASCMN,x0) 
      CALL prvec(ASCMN,x0,omvec) 
      CALL lenght(omvec,sinom) 
      CALL ptriple(ASCMN,x0,Nast,test) 
      IF (test.gt.0.d0) THEN 
         sinom = sinom 
      ELSEIF (test.lt.0.d0) THEN 
         sinom = -sinom 
      ELSE 
         WRITE(*,*)' !!!!!!!!!!!!!!!!!' 
         WRITE(*,*)' ERROR!!! ERROR!!!' 
         WRITE(*,*)'sinom = ',sinom 
         WRITE(*,*)' !!!!!!!!!!!!!!!!!' 
      ENDIF 
      mutom = datan2(sinom,cosom) 
!     ==================================================================
                                                                        
!     ==================================================================
   18 CONTINUE 
!     ==================================================================
                                                                        
      RETURN 
    END SUBROUTINE mutualrefcha
                                                                        

!     ******************************************************************
!     compute position vector (x,y,z) in cartesian coordinates
      SUBROUTINE compxyz(a,e,i,om,Omnod,u,x) 
        IMPLICIT NONE 
        DOUBLE PRECISION,INTENT(IN) :: a,e,i,om,Omnod,u 
        DOUBLE PRECISION,INTENT(OUT) :: x(3) 

!     ======= end interface ============================================
        DOUBLE PRECISION :: beta 
        DOUBLE PRECISION :: x0(3),x1(3),x2(3) 
!     ==================================================================
                                                                        
        beta = dsqrt(1-e**2) 
                                                                        
!     ================== FIRST STEP ===================                 
        x0(1) = a*(cos(u)-e) 
        x0(2) = a*beta*sin(u) 
        x0(3) = 0 
!     =================================================                 
                                                                        
!     ================== SECOND STEP ==================                 
        x1(1) = cos(om)*x0(1) - sin(om)*x0(2) 
        x1(2) = sin(om)*x0(1) + cos(om)*x0(2) 
        x1(3) = x0(3) 
!     =================================================                 
                                                                        
!     ================== THIRD STEP ===================                 
        x2(1) = x1(1) 
        x2(2) = cos(i)*x1(2) - sin(i)*x1(3) 
        x2(3) = sin(i)*x1(2) + cos(i)*x1(3) 
!     =================================================                 
                                                                        
!     ================== FOURTH STEP ==================                 
        x(1) = cos(Omnod)*x2(1) - sin(Omnod)*x2(2) 
        x(2) = sin(Omnod)*x2(1) + cos(Omnod)*x2(2) 
        x(3) = x2(3) 
!     =================================================                 
                                                                        
        RETURN 
      END SUBROUTINE compxyz
                                                                        
!     ******************************************************************
      SUBROUTINE lenght(x,len) 
        IMPLICIT NONE 
        DOUBLE PRECISION,INTENT(IN) :: x(3) 
        DOUBLE PRECISION,INTENT(OUT) :: len 
!     ===== end interface ===============

        len = dsqrt(x(1)**2 + x(2)**2 + x(3)**2) 
        
        RETURN 
      END SUBROUTINE lenght
                                         
                               
!     ******************************************************************
!     get normalized vector x/|x|
      SUBROUTINE normalize(x,xnorm) 
        IMPLICIT NONE 
        DOUBLE PRECISION,INTENT(IN) :: x(3)
        DOUBLE PRECISION,INTENT(OUT) :: xnorm(3)
!     ===== end interface ===============
        DOUBLE PRECISION :: normx 
        
        CALL lenght(x,normx) 
        xnorm(1) = x(1)/normx 
        xnorm(2) = x(2)/normx 
        xnorm(3) = x(3)/normx 
        
        RETURN 
      END SUBROUTINE normalize
                                                                        

!     ******************************************************************
!     compute triple product  <(x X y),z>
      SUBROUTINE ptriple(x,y,z,prtr) 
        IMPLICIT NONE 
        DOUBLE PRECISION,INTENT(IN) :: x(3),y(3),z(3)
        DOUBLE PRECISION,INTENT(OUT) :: prtr 
!     ==== end interface =================
        DOUBLE PRECISION :: prscal 
        EXTERNAL prscal 
        DOUBLE PRECISION :: prvtmp(3) ! auxiliary
        
        CALL prvec(x,y,prvtmp) 
        prtr = prscal(prvtmp,z) 

        RETURN 
      END SUBROUTINE ptriple
