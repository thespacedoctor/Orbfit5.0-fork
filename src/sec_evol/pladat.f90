      SUBROUTINE pladat 
      IMPLICIT NONE 
! planet data                                                           
      INCLUDE 'pldata.h90' 
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
