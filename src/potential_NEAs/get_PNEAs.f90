PROGRAM get_PNEAs
  IMPLICIT NONE
  INTEGER,PARAMETER :: dkind=KIND(1.d0) 
  INTEGER,PARAMETER :: nastx=1000000
  CHARACTER*9 name
  REAL(KIND=dkind) :: a,emin,emax,max_pd,qmin,qmax,lf
  REAL(KIND=dkind) :: aux5,aux6,aux7,aux8,aux9,aux10
  INTEGER :: h

  max_pd=1.3d0
  OPEN(10,file='neos.info',status='old')
  OPEN(9,file='emaxval',status='unknown')
  DO h=1,nastx
     READ(10,105,END=33) name,a,emin,emax,aux5,aux6,aux7,aux8,aux9,aux10,lf
     qmin=a*(1.d0-emax)
     qmax=a*(1.d0-emin)
     IF(qmin.le.max_pd)THEN
           WRITE(*,106) name,qmin,qmax,lf
           WRITE(9,107) name,a*(1.d0+emax),a*(1.d0-emax)
     ENDIF
  ENDDO
  WRITE(*,*)'error! nastx too small:',nastx
33 CONTINUE

!     name a,emin,emax,incmin,incmax,g,s,lf,ncros(i)(i=2,6),eo,i0            
105 FORMAT(a9,1x,f7.4,1x,f6.4,1x,f6.4,1x,f7.3,1x,f7.3,1x,f7.3,    &
         &        1x,f7.3,1x,f10.3,1x,f10.3,1x,f10.3,1x,i2,1x,i2,1x,i2,1x,  &
         &        i2,1x,i2,3x,f6.4,1x,f7.3)                                 

106 FORMAT(1x,a9,2(3x,f10.5),2x,f10.3)
107 FORMAT(1x,a9,2(3x,f10.5))

  CLOSE(9)
  CLOSE(10)

END PROGRAM get_PNEAs
