PROGRAM add_H
  IMPLICIT NONE
  CHARACTER*9 :: name,candname
  REAL(KIND=8) :: qmin,qmax,qmin_est,h_mag
  CHARACTER*8 :: type !nf = numbered, found in neolist
                      !nn = numbered, not in neolist
                      !mf = multiopp, found in neolist
                      !mn = multiopp, not in neolist
  INTEGER,PARAMETER :: nmax=100000
  INTEGER :: h,k
  REAL(KIND=8) :: lf
  REAL(KIND=8),PARAMETER :: eps=1.d-6

  OPEN(3,file='PNEAs',status='unknown')

! NUMBERED FOUND IN NEOLIST
  OPEN(1,file='found_in_neolist.numb',status='old')
  type = 'numb cur'
  DO h=1,nmax
     READ(1,*,END=33) name,qmin,qmax,lf
     OPEN(2,file='candidate_PNEAs.numb',status='old')
     DO k=1,nmax 
        READ(2,*,END=34) candname,qmin_est,h_mag
        IF(candname.eq.name)THEN
           IF(ABS(lf).lt.eps)THEN
              write(3,106)name,qmin,qmax,h_mag,type,'C'
           ELSE
              write(3,106)name,qmin,qmax,h_mag,type,'L'
           ENDIF
        ENDIF
     ENDDO
     write(*,*)'error!'
     STOP
34   CONTINUE
     CLOSE(2)     
  ENDDO
  write(*,*)'error!'
  STOP
33 CONTINUE
  CLOSE(1)

! NUMBERED NOT IN NEOLIST
  OPEN(1,file='not_in_neolist.numb',status='old')
  type = 'numb pot'
  DO h=1,nmax
     READ(1,*,END=35) name,qmin,qmax,lf
     OPEN(2,file='candidate_PNEAs.numb',status='old')
     DO k=1,nmax 
        READ(2,*,END=36) candname,qmin_est,h_mag
        IF(candname.eq.name)THEN
           IF(ABS(lf).lt.eps)THEN
              write(3,106)name,qmin,qmax,h_mag,type,'C'
           ELSE
              write(3,106)name,qmin,qmax,h_mag,type,'L'
           ENDIF
        ENDIF
     ENDDO
     write(*,*)'error!'
     STOP
36   CONTINUE
     CLOSE(2)     
  ENDDO
  write(*,*)'error!'
  STOP
35 CONTINUE
  CLOSE(1)

! MULTIOPP FOUND IN NEOLIST
  OPEN(1,file='found_in_neolist.mopp',status='old')
  type = 'mopp cur'
  DO h=1,nmax
     READ(1,*,END=37) name,qmin,qmax,lf
     OPEN(2,file='candidate_PNEAs.mopp',status='old')
     DO k=1,nmax 
        READ(2,*,END=38) candname,qmin_est,h_mag
        IF(candname.eq.name)THEN
           IF(ABS(lf).lt.eps)THEN
              write(3,106)name,qmin,qmax,h_mag,type,'C'
           ELSE
              write(3,106)name,qmin,qmax,h_mag,type,'L'
           ENDIF
        ENDIF
     ENDDO
     write(*,*)'error!'
     STOP
38   CONTINUE
     CLOSE(2)     
  ENDDO
  write(*,*)'error!'
  STOP
37 CONTINUE
  CLOSE(1)

! MULTIOPP NOT IN NEOLIST
  OPEN(1,file='not_in_neolist.mopp',status='old')
  type = 'mopp pot'
  DO h=1,nmax
     READ(1,*,END=39) name,qmin,qmax,lf
     OPEN(2,file='candidate_PNEAs.mopp',status='old')
     DO k=1,nmax 
        READ(2,*,END=40) candname,qmin_est,h_mag
        IF(candname.eq.name)THEN
           IF(ABS(lf).lt.eps)THEN
              write(3,106)name,qmin,qmax,h_mag,type,'C'
           ELSE
              write(3,106)name,qmin,qmax,h_mag,type,'L'
           ENDIF
        ENDIF
     ENDDO
     write(*,*)'error!'
     STOP
40   CONTINUE
     CLOSE(2)     
  ENDDO
  write(*,*)'error!'
  STOP
39 CONTINUE
  CLOSE(1)


  CLOSE(3)
106 FORMAT(1x,a9,2(3x,f10.5),2x,f8.2,3x,a8,3x,a1)


END PROGRAM add_H
