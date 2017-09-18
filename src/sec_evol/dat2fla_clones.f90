! ***********************************************************************
! Read vast_*.dat, computes the mean values of the evolutions
! and converts the results into KEP elements in a flat format for matlab
! ***********************************************************************
PROGRAM dat2fla_clones
  USE fund_const
  USE orbit_elements
  IMPLICIT NONE
! ----------------------------------------------------------------
  INTEGER,PARAMETER :: secmax=100000,nclox=100
  REAL(KIND=dkind) :: tevol(secmax,nclox),evolast(secmax,6,nclox)
  INTEGER :: ngf(secmax,nclox),nevol,nclones
  REAL(KIND=dkind) :: t0,tephast0 
  TYPE(orbit_elem) :: elem
  TYPE(orbit_elem) :: el_aux1,el_aux2 ! auxiliary
  REAL(KIND=dkind) :: mean_el(secmax,5),std_el(secmax,5)

  CHARACTER*60 :: commen,elefil,equele
  CHARACTER*3 :: numj
  INTEGER :: iundat,iunfla,iunequ,le,i,j,h,fail_flag
! ----------------------------------------------------------------
  rhs=1
  write(*,*)'Insert no. of clones'
  READ(*,*)nclones
  IF(nclones.le.0)THEN
     write(*,*)'no. of clones should be >0'
     STOP
  ENDIF
! output file
  CALL filnam('.','vastdat','fla',elefil,le)  
  CALL filopn(iunfla,elefil(1:le),'unknown')

  CALL filnam('.','vastEQUdat','fla',equele,le)  
  CALL filopn(iunequ,equele(1:le),'unknown')

! Read input file 'vast_*.dat'
  DO j=1,nclones
     WRITE (numj,'(I3)') j
     numj = adjustl(numj)
     CALL filnam('.','vast_'//trim(numj),'dat',elefil,le)  
     WRITE(*,*) elefil(1:le)
     CALL filopn(iundat,elefil(1:le),'old')
     CALL skip(iundat,1) !skipping asterid name
     IF(j.eq.1)THEN
        CALL reaflc(iundat,'t0',t0,commen) !reference time and comment 
        tephast0=t0-2400000.5d0
        CALL skip(iundat,2)
     ELSE
        CALL skip(iundat,3) !initial time is equal to all input files
     ENDIF
     DO i=1,secmax
        READ(iundat,*,END=333) tevol(i,j) !evolution time (yr)
        READ(iundat,100) evolast(i,1:6,j),ngf(i,j)
     ENDDO
     write(*,*)'do loop too short!'
     STOP
333  CONTINUE
     CALL filclo(iundat,' ')
  ENDDO
  nevol=i-1
  write(*,*)'no. evol ',nevol
  write(*,*)'final time ',tevol(nevol,1)
100 FORMAT(f12.9,4f11.7,f11.8,1x,i9,1p,e12.4)
! output header
  WRITE(iunfla,101) tephast0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
  WRITE(iunequ,102) tephast0,0.,0.,0.,0.,0.,0.,0.,0.
! Computing the arithmetic mean and standard deviation 
! in EQU elements, converting into KEP 
! and writing output in file 'vastdat.fla'
  mean_el = 0.d0
  std_el = 0.d0
  elem=undefined_orbit_elem
  el_aux1=undefined_orbit_elem
  el_aux2=undefined_orbit_elem
  DO i=1,nevol
     elem%coo='EQU'
     elem%t=tevol(i,1)
     elem%coord(6)=evolast(i,6,1)
     elem%coord(1:5)=0.d0
     DO j=1,nclones
        elem%coord(1:5)=elem%coord(1:5)+evolast(i,1:5,j)
     ENDDO
     elem%coord(1:5)=elem%coord(1:5)/nclones
     mean_el(i,1:5) =  elem%coord(1:5)
     CALL coo_cha(elem,'KEP',elem,fail_flag)
     IF(fail_flag.lt.5)THEN
     ELSE
        write(*,*)'fail_flag=',fail_flag,'stopping program'
        STOP
     ENDIF
!     WRITE(iunfla,101) elem%t,elem%coord(5)*degrad,elem%coord(4)*degrad, &
!          & elem%coord(2),elem%coord(3)*degrad
!  ENDDO
!101 FORMAT(f10.2,2x,f11.5,2x,f11.5,2x,f10.7,2x,f13.7)

! computing standard deviation
!     DO i=1,nevol
     DO h=1,5
        DO j=1,nclones
           std_el(i,h) = std_el(i,h) + (mean_el(i,h)-evolast(i,h,j))**2
        ENDDO
        std_el(i,h) = sqrt(std_el(i,h)/nclones)
     ENDDO
     el_aux1%coo='EQU'
     el_aux1%t=tevol(i,1)
     el_aux1%coord(6)=evolast(i,6,1)
     el_aux1%coord(1:5)=mean_el(i,1:5) + 3.d0*std_el(i,1:5)
     CALL coo_cha(el_aux1,'KEP',el_aux1,fail_flag)
     IF(fail_flag.lt.5)THEN
     ELSE
        write(*,*)'fail_flag=',fail_flag,'stopping program'
        STOP
     ENDIF
     
     el_aux2%coo='EQU'
     el_aux2%t=tevol(i,1)
     el_aux2%coord(6)=evolast(i,6,1)
     el_aux2%coord(1:5)=mean_el(i,1:5) - 3.d0*std_el(i,1:5)
     CALL coo_cha(el_aux2,'KEP',el_aux2,fail_flag)
     IF(fail_flag.lt.5)THEN
     ELSE
        write(*,*)'fail_flag=',fail_flag,'stopping program'
        STOP
     ENDIF
     
     WRITE(iunfla,101) elem%t,elem%coord(5)*degrad,elem%coord(4)*degrad, &
          & elem%coord(2),elem%coord(3)*degrad, &
          & el_aux1%coord(5)*degrad,el_aux1%coord(4)*degrad, &
          & el_aux1%coord(2),el_aux1%coord(3)*degrad, &
          & el_aux2%coord(5)*degrad,el_aux2%coord(4)*degrad, &
          & el_aux2%coord(2),el_aux2%coord(3)*degrad
     WRITE(iunequ,102) elem%t,mean_el(i,2:5),std_el(i,2:5)
     
  ENDDO
101 FORMAT(f10.2,3(2x,es11.5,2x,es11.5,2x,es11.5,2x,es11.5))
102 FORMAT(f10.2,2(2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5))

  CALL filclo(iunfla,' ')
  CALL filclo(iunequ,' ')
  
END PROGRAM dat2fla_clones
