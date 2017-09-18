MODULE secular_evolution
  USE fund_const
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sec_evol,undefined_sec_evol


! ---------------------------------------
  TYPE sec_evol
     REAL(KIND=dkind) :: t
     REAL(KIND=dkind) :: a,om,nod,ecc,inc
     INTEGER :: cf   ! crossing flag, read from pro files  
                     ! cf=0 when crossings omega=k*pi/2
  END TYPE sec_evol
! ---------------------------------------
! undefined secular evolution
  TYPE(sec_evol), PARAMETER :: undefined_sec_evol = SEC_EVOL( &
&  0.d0,           & ! time (MJD)
&  0.d0,           & ! a    (AU)
&  0.d0,           & ! omega (deg)
&  0.d0,           & ! Omnod (deg)
&  0.d0,           & ! eccentricity
&  0.d0,           & ! inclination (deg)
&  55              & ! default crossing flag
&     )  ! 


END MODULE secular_evolution
