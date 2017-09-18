! *********************************************************
! SKIP THIS PART
21 DO i = inpl,ioupl
     el_pla(i)=undefined_orbit_elem
     el_pla(i)%t=t1
     CALL placar(i,el_pla(i)%t,el_pla(i)%coord,1,force_circ)
     IF(force_circ)THEN
        el_pla(i)%coo='KEP'
         GOTO 15
     ENDIF
     el_pla(i)%coo='CAR'
     CALL coo_cha(el_pla(i),'KEP',el_pla(i),fail_flag)
     IF(fail_flag.lt.5)THEN
!        write(*,*)'time,i,elpl KEP:',el_pla(i)%t,i,el_pla(i)%coord
     ELSE
        write(*,*)'fail_flag=',fail_flag,'stopping program'
        STOP
     ENDIF
15   CONTINUE
  ENDDO
! *********************************************************

!        CALL placar(i,el_pla(i)%t,el_pla(i)%coord,1,force_circ)
!        IF(force_circ)THEN
!        el_pla(i)%coo='KEP'
!         GOTO 16
!     ENDIF
!        el_pla(i)%coo='CAR' 
!        CALL coo_cha(el_pla(i),'KEP',el_pla(i),fail_flag)
!        IF(fail_flag.lt.5)THEN
!           ! write(*,*)'j,elpl KEP:',j,elpl%coord
!        ELSE
!           write(*,*)'fail_flag=',fail_flag,'stopping program'
!           STOP
!        ENDIF
! 16 CONTINUE
