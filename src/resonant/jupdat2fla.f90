! Legge vpla.fil e lo converte in una tabella leggibile da matlab
! in cui si mettono solo i dati di Jupiter
PROGRAM jupdat2fla
  USE fund_const
  USE orbit_elements
  IMPLICIT NONE
  INTEGER,PARAMETER :: secmax=100000
  REAL(KIND=dkind) :: tevol(secmax),evoljup(secmax,6)
  TYPE(orbit_elem) :: elem
  INTEGER :: iundat,iunfla,i,fail_flag 
  INTEGER :: ngf(secmax),le
! To read the header
  REAL(KIND=dkind) :: t0,tephpla0 
  CHARACTER*60 :: commen,elefil
! *************************************************

  rhs=1

  CALL filopn(iundat,'vpla.dat','old')
  CALL filnam('.','vjupdat','fla',elefil,le)  
  CALL filopn(iunfla,elefil(1:le),'unknown')

! Reading the header
  CALL skip(iundat,1) !skipping planet names
  CALL reaflc(iundat,'t0',t0,commen) !reference time and comment 
  tephpla0=t0-2400000.5d0
! file header
  WRITE(iunfla,101) tephpla0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
  CALL skip(iundat,5)

  elem=undefined_orbit_elem
  
! Reading evolution from file vpla.dat and writing it in file vplafil.fla
  DO i=1,secmax
     READ(iundat,*,END=333) tevol(i) ! difference of time (yr) from tephpla0
     CALL skip(iundat,3)
     READ(iundat,100) evoljup(i,1:6),ngf(i)
     CALL skip(iundat,1)
100  FORMAT(f12.9,4f11.7,f11.8,1x,i9,1p,e12.4)
     elem%coo='EQU'
     elem%coord=evoljup(i,1:6)
     elem%t=tevol(i)
     CALL coo_cha(elem,'KEP',elem,fail_flag)
     IF(fail_flag.lt.5)THEN
     ELSE
        write(*,*)'fail_flag=',fail_flag,'stopping program'
        STOP
     ENDIF
     WRITE(iunfla,101) tevol(i),elem%coord(5)*degrad,elem%coord(4)*degrad, &
          & elem%coord(2),elem%coord(3)*degrad,elem%coord(1), &
          & elem%coord(6)*degrad
103  FORMAT(f10.2,2x,f12.9,2x,4(f11.7,2x))
101  FORMAT(f10.2,2x,f11.5,2x,f11.5,2x,f10.7,2x,f13.7,2x,f13.7,2x,f11.5)
  ENDDO
  write(*,*)'do loop too short!'
  STOP
333 CONTINUE

  CALL filclo(iundat,' ')
  CALL filclo(iunfla,' ')

END PROGRAM jupdat2fla
