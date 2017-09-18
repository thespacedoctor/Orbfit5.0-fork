! =========================================================             
! perturbing function                                                   
MODULE perturbing_function
  USE fund_const
  IMPLICIT NONE
  PRIVATE
  REAL(KIND=dkind) :: aa,apl,gmp,ee
  REAL(KIND=dkind) :: co,so,beta,ci,si 
  REAL(KIND=dkind) :: cu,su,x,y,z 
  REAL(KIND=dkind) :: usgn,susgn,cusgn,usgn1,susgn1,cusgn1 
!  COMMON/elems/co,so,ci,si,ee,beta,aa,apl,gmp,usgn,usgn1 

! ap = semimajor axis of the planet, ky = Gauss constant 
! gm = (masse pianeti)/(massa Sole)
  INTEGER :: npl,inpl,ioupl,ndum
  INTEGER,PARAMETER :: nplax=10,nplax2=20 
  REAL(KIND=dkind) :: ap(nplax),gm(nplax)
  REAL(KIND=dkind) :: ky,bigg 
!  COMMON/pldata/gm,ap,ky,ndum,npl,inpl,ioupl,bigg

!   REAL(KIND=dkind) :: cu,su,x,y,z 
!  COMMON/elele/x,y,z,cu,su 

  PUBLIC :: perturb,pladat
  PUBLIC :: gm,ap,ky,ndum,npl,inpl,ioupl,bigg
  PUBLIC :: aa,apl,gmp,ee,co,so,beta,ci,si 
  PUBLIC :: cu,su,x,y,z 
  PUBLIC :: usgn,susgn,cusgn,usgn1,susgn1,cusgn1 
CONTAINS
  
  SUBROUTINE perturb(om,g,zl,a,ddd,eee,nnn) 
!    USE fund_const
!    IMPLICIT NONE 
! output: right hand side,errors,number eval                            
    DOUBLE PRECISION,INTENT(OUT) :: ddd,eee 
    INTEGER,INTENT(OUT) :: nnn 
! input: asteroid elements                                              
    DOUBLE PRECISION,INTENT(IN) :: om,g,zl,a 
! ====end interface========================================             
! planet data                                                           
!      INCLUDE 'pldata.h90' 
      INTEGER n 
! ========================================================              
! per dqagp                                                             
      REAL(KIND=dkind),EXTERNAL :: ff
      INTEGER limx,limx4 
      PARAMETER (limx=300,limx4=4*limx) 
      INTEGER ier,leniw,iwork(limx),lenw,last 
! function evaluations                                                  
      INTEGER neval,npts2 
      DOUBLE PRECISION points(4) 
      DOUBLE PRECISION epsabs,epsrel,abserr,work(limx4) 
! output of dqagp                                                       
      DOUBLE PRECISION rm 
! =======================================================               
! operations done only once                                             
      ee=sqrt(1.d0-(g/ky)**2/a) 
      beta=sqrt(1.d0-ee**2) 
      ci=zl/(sqrt(a)*beta) 
      if(ci.lt.1.d0)then 
           si=sqrt(1.d0-ci**2) 
      else 
           si=0.d0 
      endif 
      co=cos(om) 
      so=sin(om) 
!      write(*,*)'a,zl,omega,e',a,zl,om*180.d0/pig,ee

! set accumulators to zero                                              
      ddd=0.d0 
      eee=0.d0 
      nnn=0 
                                                                        
!     control                                                           
!      write(*,*)'inpl,ioupl',inpl,ioupl                                
! loop on number of planets                                             
      DO 10 n=inpl,ioupl 
! dati per la quadratura corrente                                       
         aa=a 
         apl=ap(n) 
         gmp=gm(n) 
                                                                        
! ========= punti singolari usgn,usgn1 ==============                   
! definisco susgn                                                       
      susgn = -(so*beta/(1.d0+ee*co)) 
! cusgn                                                                 
      cusgn =(co + ee)/(1+ee*co) 
! usgn                                                                  
      usgn = atan2(susgn,cusgn) 
! ===================================================                   
                                                                        
! definisco susgn1                                                      
      susgn1 = (so*beta/(1.d0-ee*co)) 
! cusgn1                                                                
      cusgn1 =(ee - co)/(1-ee*co) 
! usgn1                                                                 
      usgn1 = atan2(susgn1,cusgn1) 
! ====================================================                  
                                                                        
! preparativi chiamata dqagp                                            
         epsabs=1.d-8 
         epsrel=1.d-5 
         leniw=limx 
         lenw=limx4 
         npts2=4 
         points=0.d0
         points(1)=usgn 
         points(2)=usgn1 
! perturbing function                                                   
!         CALL dqags(ff,0.d0,dpig,epsabs,epsrel,rm,abserr,neval,ier,    &
!     &          leniw,lenw,last,iwork,work)                            
         CALL dqagp(ff,-pig,pig,npts2,points,epsabs,epsrel,rm,abserr,   &
     & neval,ier,leniw,lenw,last,iwork,work)                            
         ddd=ddd+rm/(dpig**2) 
         eee=eee+abserr/(dpig**2) 
         nnn=nnn+neval 
   10 END DO 

    END SUBROUTINE perturb

    SUBROUTINE pladat 
!      IMPLICIT NONE 
! planet data                                                           
!      INCLUDE 'pldata.h90' 
      DOUBLE PRECISION gk,gk2,gjyr,gjyr2 
      INTEGER i 
! ======== PIANETI presi in considerazione ==============               
      inpl=2 
      ioupl=6 
      npl=ioupl-inpl+1 
! ============ planet semimajor axis ==================                 
      ap(1)= .3870992058d0 
      ap(2)= .7233274811d0 
      ap(3)= 1.0000036214d0 
      ap(4)= 1.5235973464d0  
      ap(5)= 5.2024107723d0  
      ap(6)= 9.5575876779d0  
      ap(7)= 19.3008879212d0  
      ap(8)= 30.2722024706d0  
      ap(9)= 39.7533710065d0  
! ============ masse dei pianeti =======================                
      gm(1)= 0.00000016601d0  
      gm(2)= 0.00000244781d0  
      gm(3)= 0.0000030404d0  
      gm(4)= 0.00000032272d0  
      gm(5)= 0.00095479d0  
      gm(6)= 0.00028589d0  
      gm(7)= 0.000043662d0  
      gm(8)= 0.000051514d0  
      gm(9)= 0.0000000073964d0 
! =========== conversioni di unita' =====================               
!  Gauss constant                                                       
      gk=0.01720209895d0 
      gk2=gk*gk 
!  conversion to internal units: 1AU, 1JYR=365.25 d(Julian year)        
      gjyr=365.25d0 
      gjyr2=gjyr*gjyr 
      bigg=gk2*gjyr2 
      ky=gk*gjyr 
      DO i =1,9 
        gm(i)=gm(i)*bigg 
      ENDDO 
      RETURN 
    END SUBROUTINE pladat

  END MODULE perturbing_function


! ===========================================================           
! subroutine che calcola la media unidimensionale della media unid.     
    DOUBLE PRECISION FUNCTION ff(u) 
      USE fund_const
      USE perturbing_function
      IMPLICIT NONE 
! anomalia eccentrica dell'asteroide                                    
      REAL(KIND=dkind),INTENT(IN) :: u 

! per dqags                                                             
      INTEGER limx,limx4 
      PARAMETER (limx=300,limx4=4*limx) 
      INTEGER neval,ier,leniw,iwork(limx),lenw,last 
      INTEGER npts2 
      DOUBLE PRECISION eta 
      DOUBLE PRECISION epsabs,epsrel,abserr,work(limx4),resul
      DOUBLE PRECISION points(3) 
      REAL(KIND=dkind),EXTERNAL :: f
! calcoli dipendenti solo da u                                          
      cu=cos(u) 
      su=sin(u) 
      x=aa*((cu-ee)*co-beta*su*so) 
      y=aa*((cu-ee)*so+beta*su*co)*ci 
      z=aa*((cu-ee)*so+beta*su*co)*si 
                                                                        
! preparativi per dqags,dqagp                                           
      epsabs=1.d-8 
      epsrel=1.d-5 
      leniw=limx 
      lenw=limx4 
      npts2=3 
      points=0.d0
      points(1)=0.d0 
                                                                        
! integrale                                                             
! parametro di confronto                                                
      eta = 1.d-7 
      IF (abs(usgn - u).le.eta)THEN 
!      call dqagsc(f,0.d0,dpig,epsabs,epsrel,resul,abserr,neval,ier, &   
!     &   leniw,lenw,last,iwork,work)                                   
      call dqagpc(f,-pig,pig,npts2,points,epsabs,epsrel,resul,         &
     & abserr,neval,ier,leniw,lenw,last,iwork,work)                     
        ELSEIF (abs(usgn - u).ge.eta)THEN 
      call dqagsc(f,-pig,pig,epsabs,epsrel,resul,abserr,neval,ier,     &
     &   leniw,lenw,last,iwork,work)                                    
      ENDIF 
! risultato                                                             
      ff=resul
    END FUNCTION ff
! ===========================================================           
! subroutine che calcola la media unidimensionale                       
    DOUBLE PRECISION FUNCTION f(up) 
      USE fund_const
      USE perturbing_function
      IMPLICIT NONE 
! anomalia (eccentrica) del pianeta                                     
      REAL(KIND=dkind),INTENT(IN) :: up 
! end interface
      REAL(KIND=dkind) :: d2,r,xst,yst 
!                                                                       
      xst=apl*cos(up) 
      yst=apl*sin(up) 
                                                                        
!  FUNZIONE INTEGRANDA                                                  
      d2=((x-xst)**2)+((y-yst)**2)+(z**2) 
      f=(gmp/sqrt(d2))*(1.d0-ee*cu) 
                                                                        
      RETURN 
    END FUNCTION f

