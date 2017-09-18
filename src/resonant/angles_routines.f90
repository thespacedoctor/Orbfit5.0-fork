
SUBROUTINE ang_interp(lambda,alpha,beta,gamma)
!  USE fund_const,ONLY: dkind
  IMPLICIT NONE
  INTEGER,PARAMETER :: dkind = KIND(1.D0)
  REAL(KIND=dkind),INTENT(IN) :: lambda,alpha,beta
  REAL(KIND=dkind),INTENT(OUT) :: gamma
! end interface
  REAL(KIND=dkind) :: alpha1,beta1,bma,pridif,princ

  alpha1 = princ(alpha)
  beta1 = princ(beta)

  bma=pridif(beta1,alpha1)
  gamma = lambda*bma+alpha1
  gamma = princ(gamma)
END SUBROUTINE ang_interp

! Assuming that the smaller angle between e^(i*alpha) and e^(i*beta) 
! gives the orientation, compute gamma = alpha+lambda*(beta-alpha)
SUBROUTINE ang_interp_var(lambda,alpha,beta,gamma)
!  USE fund_const,ONLY: dkind
  IMPLICIT NONE
  INTEGER,PARAMETER :: dkind = KIND(1.D0)
  REAL(KIND=dkind),INTENT(IN) :: lambda,alpha,beta
  REAL(KIND=dkind),INTENT(OUT) :: gamma
! end interface
  REAL(KIND=dkind) :: ca,sa,cb,sb
  REAL(KIND=dkind) :: bma,cbma,sbma
  REAL(KIND=dkind) :: cg,sg

! check lambda
  IF(lambda.lt.-1.d-10.OR.lambda.gt.1.d0+1.d-10)THEN
     write(*,*)'angles_routines: error! lambda=',lambda
     STOP
  ENDIF
  ca=cos(alpha);sa=sin(alpha)
  cb=cos(beta);sb=sin(beta)
! bma=beta-alpha
  cbma=cb*ca+sb*sa
  sbma=sb*ca-sa*cb
  bma=ATAN2(sbma,cbma) 
  gamma = alpha + lambda*bma
  cg=cos(gamma);sg=sin(gamma)
  gamma=ATAN2(sg,cg)

END SUBROUTINE ang_interp_var
