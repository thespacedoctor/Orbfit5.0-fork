PROGRAM mindist
  USE fund_const
  USE output_control
  IMPLICIT NONE
  INTEGER :: iunlst,iunout,iudmin,iunmars,iunsta
  CHARACTER(LEN=9) :: name
  CHARACTER(LEN=60) :: categ
  INTEGER,PARAMETER :: nlst=10000
  INTEGER,PARAMETER :: lenx=10000
  INTEGER :: lnam,lcat
  INTEGER :: j,h,nt
  REAL(KIND=dkind),DIMENSION(8,lenx) :: dist
  REAL(KIND=dkind) :: aa
  REAL(KIND=dkind),DIMENSION(8) :: distances
  REAL(KIND=dkind),DIMENSION(7) :: d_nomars
  REAL(KIND=dkind) :: mer,ven,ear,mar,jup,sat,ura,net
  REAL(KIND=dkind),DIMENSION(lenx) :: mind,mind_nomars
  LOGICAL :: found

  WRITE(*,*)'input category:'
  READ(*,*)categ
  CALL rmsp(categ,lcat)
  CALL filopn(iunlst,categ(1:lcat)//'.list','old')
  CALL filopn(iunout,'mindist','unknown')
  CALL filopn(iunmars,'mindist_no_mars','unknown')
  CALL filopn(ierrou,'mindist.err','unknown')
  CALL filopn(iunsta,'stable.list','unknown')

  DO j=1,nlst
     mind=1.d10 !initialization
     READ(iunlst,('(a9)'),end=33) name
!     WRITE(*,*)'name=',name
     CALL rmsp(name,lnam)
     INQUIRE(file=name(1:lnam)//'.dmin',exist=found)
     IF(.not.found)THEN
        WRITE(ierrou,*)'the file',name(1:lnam)//'.dmin',' does not exist'
        CYCLE
     ENDIF
     CALL filopn(iudmin,name(1:lnam)//'.dmin','old')
     DO h=1,lenx
        READ(iudmin,*,END=34) nt,aa,dist(1,h),dist(2,h),dist(3,h),dist(4,h),dist(5,h),dist(6,h),dist(7,h),dist(8,h)
        distances(1:8)=abs(dist(1:8,h))
        mind(h) = minval(distances)
        d_nomars(1:3)=distances(1:3)
        d_nomars(4:7)=distances(5:8)
        mind_nomars(h) = minval(d_nomars)
!        write(*,*)'dist', distances
!        write(*,*)'mind(h)',mind(h)
!        read(*,*)
     ENDDO
34   CONTINUE
     CALL filclo(iudmin,' ')
        
     mer=minval(abs(dist(1,1:h-1)))
     ven=minval(abs(dist(2,1:h-1)))
     ear=minval(abs(dist(3,1:h-1)))
     mar=minval(abs(dist(4,1:h-1)))
     jup=minval(abs(dist(5,1:h-1)))
     sat=minval(abs(dist(6,1:h-1)))
     ura=minval(abs(dist(7,1:h-1)))
     net=minval(abs(dist(8,1:h-1)))
     WRITE(iunout,100)name,mer,ven,ear,mar,jup,sat,ura,net, &
          & minval(mind(1:h-1)),minval(mind_nomars(1:h-1))
     WRITE(*,100)name,mer,ven,ear,mar,jup,sat,ura,net, &
          & minval(mind(1:h-1)),minval(mind_nomars(1:h-1))
!     read(*,*)
     IF(jup.ge.1.d0.AND.ear.ge.0.1d0.AND. &
          & minval(mind_nomars(1:h-1)).ge.0.1d0)THEN
        WRITE(iunsta,110)name,aa,mer,ven,ear,mar,jup,sat,ura,net, &
             & minval(mind(1:h-1)),minval(mind_nomars(1:h-1))
     ENDIF
  ENDDO
  WRITE(*,*)'to small nlst:',nlst
33 CONTINUE
100 FORMAT (a9,1x,10(f10.5,1x)  )
110 FORMAT (a9,1x,11(f10.5,1x)  )

CALL filclo(ierrou,' ')
CALL filclo(iunout,' ')
CALL filclo(iunmars,' ')
CALL filclo(iunlst,' ')
CALL filclo(iunsta,' ')

END PROGRAM mindist
