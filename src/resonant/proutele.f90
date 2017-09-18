! *************************************************************     
!     OUTPUT WRITING SUBROUTINE                                         
! *************************************************************
!SUBROUTINE proutele(iun,t,y,a,nit,dn,nnod)
SUBROUTINE proutele(iun,t,y,nit,dn,nnod)
  USE fund_const
  IMPLICIT NONE
! ===================== INTERFACE =============================
! y(1)=omega; y(2)=G; y(3)=Omega; y(4)=zl, y(5)=S, y(6)=sigma                                    
  REAL(KIND=dkind),INTENT(IN) :: t,y(6),dn 
  INTEGER,INTENT(IN) :: iun,nit,nnod
! ======================= END INTERFACE =======================
! to check the Hamiltonian                                          
  REAL(KIND=dkind) :: Hbar,error1 
! =============================================================    
  INTEGER numb1 
  REAL(KIND=dkind) :: a,e,beta,ci,si,ainc 
! =============================================================             
  INCLUDE 'pldata.h90' 
! ************************************************************* 
  a= (y(5)/ky)**2
  e = dsqrt(1.d0-(y(2)/ky)**2/a) 
  beta  = dsqrt(1.d0-e**2) 
  
  ci=y(4)/(ky*dsqrt(a)*beta) 
  if(ci.lt.1.d0)then 
     si=dsqrt(1.d0-ci**2) 
  else 
     si=0.d0 
  endif
  
  ainc=datan2(si,ci)*degrad 
                                                                        
! ==================== WRITING DATA (time in yrs) ==================
  IF(iun.gt.0)THEN 
! t,omega,Omnod,a,e,I,sigma,nit,dn,nnod
     WRITE(iun,100)t,y(1)*degrad,y(3)*degrad,a,e,ainc,y(6),                &
          &        nit,dn,nnod
! ============================================================= 
!     controllo del valore dell'Hamiltoniana                            
!     CALL prclevpert(y(1),y(2),zl,a,Hbar,error1,numb1)                 
!     WRITE(iun,*)'proutele Hbar,omega,G=',Hbar,degrad*y(1),y(2)        
!     WRITE(*,*)'proutele Hbar,omega,G',Hbar,degrad*y(1),y(2)           
!     WRITE(*,*)'                                       '               
!     end control 
! =============================================================  

  ELSEIF(iun.eq.0)THEN 
     WRITE(*,100)t,y(1)*degrad,y(3)*degrad,a,e,ainc,y(6),                  &
          &        nit,dn,nnod
  ENDIF
  
100 FORMAT(f10.2,1x,f11.5,1x,f11.5,1x,f10.7,1x,f10.7,1x,f13.7,1x,f13.7,i3,1x,           &
         &     f13.10,1x,i3)

  RETURN 
END SUBROUTINE proutele
