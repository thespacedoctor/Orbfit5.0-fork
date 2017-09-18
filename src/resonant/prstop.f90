! ******************************************************************
! **** Subroutine deciding when informations are enough to stop ****
! ******************************************************************
! om,omzero are in radians                                          
SUBROUTINE prstop(iunpre,tast0,omfreq,t,y,stopfl) 
  USE fund_const
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: iunpre 
  REAL(KIND=dkind),INTENT(IN) :: tast0,omfreq,t,y(6)
  INTEGER,INTENT(INOUT) :: stopfl 
!     stopfl=0  is used in the first step                               
!     stopfl=1  sector change has not been yet                          
!     stopfl=2  after changing the first sector                         
!     (after the transition band)                                       
!     stopfl=-1 stop                                                    
! ====================== END INTERFACE =============================
  INCLUDE 'pldata.h90' 
  INTEGER nsx 
  PARAMETER (nsx=1000) 
  REAL(KIND=dkind) :: om(nsx),e(nsx),a(nsx),omzero,ezero,azero
  REAL(KIND=dkind) :: omnod(nsx),incl(nsx),sigma(nsx),inclzero,omnodzero,sigmazero
  REAL(KIND=dkind) :: nintomz,chomz 
  REAL(KIND=dkind) :: y1(6),y2(6),targ,yt(6),dt 
  REAL(KIND=dkind) :: tzero,t1,t2,tim(nsx) 
  INTEGER nst,sect(nsx) 
  LOGICAL trans(nsx) 
  INTEGER i,j,prevsec,jump1,jump2,postsec,finsec,nprev 
  REAL(KIND=dkind) :: dist,ee 
! eccentricity and inclination at k*pig/2                           
  REAL(KIND=dkind) :: e1,omeg1,e2,omeg2,omnod1,omnod2,incl1,incl2,a1,a2,sigma1,sigma2 
! comparison parameters                                             
  REAL(KIND=dkind) :: epsilon,radiuse,radiusom,saf  
  REAL(KIND=dkind) :: fntarpt 
! need to save local variable nst for next calls.                   
  SAVE 
! ******************************************************************
! OSS:sarebbe meglio mettere un epsilon che dip dalla evel(velocita' di e)
  epsilon=0.1 
  saf=epsilon*(pig/100) 
! creare un radius diverso per e                                    
  radiusom=pig/100 
  radiuse=0.0001 
! control                                                           
! write(*,*)' saf=',saf,' radius=',radiusom, radiuse                
! ================================================================= 
                                                                        
! ==================== FIRST STEP =================================
                                                                        
  IF(stopfl.eq.0)THEN 
! initialisation
     azero=(y(5)/ky)**2
     omzero=y(1) 
     ezero=sqrt(1.d0-y(2)**2/(ky**2*azero))
     omnodzero=y(3)
     inclzero=dacos(y(4)/y(2))
     sigmazero=y(6)
     tzero=t 
     nst=0
     prevsec=0 
     postsec=0 
     finsec=0 
     jump1=0 
     jump2=0 

     nintomz=nint(omzero*degrad) 
! check                                                             
     write(*,*)'omzero,nint(omzero)',omzero*degrad,nintomz
     write(*,*)'saf=',saf*degrad 
                                                                        
! compute nint(omzero*degrad) between -180 and 180 degrees          
     CALL choose_deg(nintomz,chomz) 
! control                                                           
     write(*,*)'chomz=',chomz 
     IF(abs(chomz+180).le.saf*degrad)THEN 
! stopfl is kept zero                                               
        WRITE(*,*)'INITIAL CONDITION IN TRANSITION BAND' 
     ELSEIF(abs(chomz+90).le.saf*degrad)THEN 
! stopfl is kept zero                                               
        WRITE(*,*)'INITIAL CONDITION IN TRANSITION BAND' 
     ELSEIF(abs(chomz).le.saf*degrad)THEN 
! stopfl is kept zero                                               
        WRITE(*,*)'INITIAL CONDITION IN TRANSITION BAND' 
     ELSEIF(abs(chomz-90).le.saf*degrad)THEN 
! stopfl is kept zero                                               
        WRITE(*,*)'INITIAL CONDITION IN TRANSITION BAND' 
     ELSEIF(abs(chomz-180).le.saf*degrad)THEN 
! stopfl is kept zero                                               
        WRITE(*,*)'INITIAL CONDITION IN TRANSITION BAND' 
     ELSE 
! set stopfl to 1                                                   
        stopfl=1 
     ENDIF
! control                                                           
     WRITE(*,*)'stopfl0=',stopfl 
! ================================================================= 
  ELSEIF(stopfl.eq.1)THEN 
! normal operations, before the first sector change                 
! check for closed curve in e,omega plane                           
     ee=sqrt(1.d0-y(2)**2/(ky**2*azero)) 
     dist=max((abs(y(1)-omzero))/radiusom,(abs(ee-ezero))/radiuse)

! ======== ASYMMETRIC LIBRATION DETECTED! stop next time =========
! (asymmetric librations can occur only without sector changes)     
     IF(dist.le.1.d-5)THEN 
! stopfl=-1                                                         
        WRITE(*,*)'%%%%% ASYMMETRIC LIBRATION DETECTED! %%%%%' 
        WRITE(*,*)'time=',t,'eccentricity=',ee,'omega=',y(1) 
     ENDIF
! ================================================================
  ELSEIF(stopfl.eq.2)THEN 
! normal operations, after the first sector change                  
  ELSEIF(stopfl.eq.-1) THEN 
! asymmetric libration!                                             
! ================================================================
  ENDIF
                                                                        
! control                                                           
! WRITE(*,*)' after first if stopfl=',stopfl                        
                                                                        
! ====================== STORAGE OF DATA =========================     
  nst=nst+1 
! control                                                           
  write(*,*)'nst=',nst 
  a(nst)=(y(5)/ky)**2 
  om(nst)=y(1) 
  e(nst)=sqrt(1.d0-y(2)**2/(ky**2*a(nst))) 
  omnod(nst)=y(3)
  incl(nst)=dacos(y(4)/y(2))
  sigma(nst)=y(6)
  tim(nst)=t 
! ================================================================
                                                                        
  CALL prsector(om,4,saf,sect,trans,nst) 
  write(*,*)'nst,sect,trans,om ',nst,sect(nst),trans(nst),om(nst) 
! check for sector transition: not in the transition band           
  IF(.not.trans(nst))THEN 
     IF(stopfl.eq.1)THEN 
! no transition yet: compare with beginning                         
        IF(sect(nst).ne.sect(1))THEN 
           write(*,*)' stopfl=1 ',(sect(i),i=1,nst) 
! transition occurred; find the last non transition band in         
! the previous sector                                               
           DO j=nst-1,1,-1 
! control                                                           
! write(*,*)'contiamo j=',j                                         
                                                                        
              IF(.not.trans(j))THEN 
                 prevsec=sect(j) 
                 postsec=sect(nst) 
                 jump1=postsec-prevsec 
                 WRITE(*,*)'prev,post,jump1 ',prevsec,postsec,jump1 
                 nprev=j 
! transition recorded: regula falsi                                 
! provvisorio                                                       
! transition recorded: regula falsi                                 
                 t1=tim(nprev) 
                 y1(1)=om(nprev) 
                 y1(2)=ky*sqrt(a(nprev)*(1.d0-e(nprev)**2)) 
                 y1(3)=omnod(nprev)
                 y1(4)=y1(2)*cos(incl(nprev))
                 y1(5)=ky*sqrt(a(nprev))
                 y1(6)=sigma(nprev)
                 t2=tim(nst) 
                 y2(1)=om(nst) 
                 y2(2)=ky*sqrt(a(nst)*(1.d0-e(nst)**2)) 
                 y2(3)=omnod(nst)
                 y2(4)=y2(2)*cos(incl(nst))
                 y2(5)=ky*sqrt(a(nst))
                 y2(6)=sigma(nst)
! targ=postsec*pig/2.d0                                             
! write(*,*)'prevsec,postsec',prevsec,postsec                       
! targ=prevsec*pig/2.d0                                             
                 targ=(fntarpt(prevsec,postsec))*pig/2.d0 
                 write(*,*)'target=',targ*degrad 
                 write(*,*)'t1,t2 prima di falsectra:',t1,t2 
                 CALL falsectra(tast0,t1,y1,t2,y2,omfreq,targ,yt,dt)
                 write(*,*)'t1,t2 dopo di falsectra:',t1,t2 
! output                                        
                 a1=(yt(5)/ky)**2
                 e1=sqrt(1.d0-yt(2)**2/(ky**2*a1)) 
                 omeg1=yt(1) 
                 omnod1=yt(3) 
                 incl1=dacos(yt(4)/yt(2))
                 sigma1=yt(6)
! control                                                           
! write(*,*)'prstop incl1:',incl1                                   
                                                                        
! scrivi l'output t2+dt, yt su file elementi propri                 
! control                                                           
! write(*,*)'&&&&&&&&&&&&&&&&&&&&&&&&&&&&'                          
! write(*,*)'CONTROLLO prima di scrivere in .pre!!!!'               
! write(*,*)'dt=',dt,'  t1=',t1,'  t2=',t2                          
! write(*,*)'&&&&&&&&&&&&&&&&&&&&&&&&&&&&' 
               
                 WRITE(iunpre,100)t2+dt,omeg1*degrad,      &
                      &   omnod1*degrad,e1,incl1*degrad,0
! control                                                           
! WRITE(*,*)'ehi!',t2+dt,omeg1*degrad,                              
!    *                    omnod1*degrad,e1,incl1*degrad,0              
                                                                        
                 stopfl=2 
                                                                        
100              FORMAT(f10.2,1x,f11.5,1x,f11.5,1x,f10.7,1x,f13.7,i3)
                 EXIT 
              ENDIF
           ENDDO
!     2             CONTINUE 
        ELSE 
! noting happened                                                   
        ENDIF
                                                                        
! take into account also the possibility of asymmetric libration    
                                                                        
     ELSEIF(stopfl.eq.2)THEN 
! already one transtion, look for the second                        
        IF(sect(nst).ne.postsec)THEN 
           write(*,*)' stopfl=2 ',(sect(i),i=1,nst),postsec 
! transition occurred; find the last non transition band            
! in the previous sector                                            
           DO j=nst-1,1,-1 
              IF(.not.trans(j))THEN 
                 finsec=sect(nst) 
                 jump2=finsec-postsec 
! WRITE(*,*)'post,fin,jump2 ',postsec,finsec,jump2                  
                 nprev=j 
! transition recorded: regula falsi                                 
                 t1=tim(nprev) 
                 y1(1)=om(nprev) 
                 y1(2)=ky*sqrt(a(nprev)*(1.d0-e(nprev)**2)) 
                 y1(3)=omnod(nprev)
                 y1(4)=y1(2)*cos(incl(nprev))
                 y1(5)=ky*sqrt(a(nprev))
                 y1(6)=sigma(nprev)
                 t2=tim(nst) 
                 y2(1)=om(nst) 
                 y2(2)=ky*sqrt(a(nst)*(1.d0-e(nst)**2)) 
                 y2(3)=omnod(nst)
                 y2(4)=y2(2)*cos(incl(nst))
                 y2(5)=ky*sqrt(a(nst))
                 y2(6)=sigma(nst)
! targ=finsec*pig/2.d0                                              
! write(*,*)'postsec,finsec',postsec,finsec                         
                 targ=(fntarpt(postsec,finsec))*pig/2.d0 
                 write(*,*)'target=',targ 
                 CALL falsectra(tast0,t1,y1,t2,y2,omfreq,targ,yt,dt)
! output                                                 
                 a2=(yt(5)/ky)**2
                 e2=sqrt(1.d0-yt(2)**2/(ky**2*a2)) 
                 omeg2=yt(1) 
                 omnod2=yt(3) 
                 incl2=dacos(yt(4)/yt(2))
                 sigma2=yt(6)
! control                                                           
! write(*,*)'prstop incl2:',incl2                                   
                                                                        
! scrivi l'output t2+dt, yt su file elementi propri                 
                 WRITE(iunpre,100)t2+dt,omeg2*degrad,    &
                      &  omnod2*degrad,e2,incl2*degrad,0               
! control                                                           
! WRITE(*,*)'ehi!',t2+dt,omeg2*degrad,                              
!  *      omnod2*degrad,e2,incl2*degrad,0              
                                                                        
                 stopfl=-1 
                 EXIT 
              ENDIF
           ENDDO
!     2             CONTINUE                                            
        ELSE 
! noting happened                                                   
        ENDIF
                                                                        
         ELSE 
! if first step, nothing to do                                      
         ENDIF
      ELSE 
! during transition phase, better to wait                           
      ENDIF
! debugging messages                                                
! WRITE(*,*)' after second if stopfl=',stopfl                       
      IF(stopfl.ne.-1)RETURN 
! write(*,*)' jump1,jump2,prev,post,fin '                           
! write(*,*) jump1,jump2, prevsec,postsec,finsec
! ================================================================
! ================ caso della CIRCOLAZIONE =======================
! ================================================================
      IF(finsec.ne.prevsec)THEN 
! qui dovrei innescare la regula falsi                              
         write(*,*)'circulation' 
         IF(jump1+jump2.eq.2)THEN 
            write(*,*)' direct' 
         ELSEIF(jump1+jump2.eq.-2)THEN 
            write(*,*)' retrograde ' 
         ELSE 
            write(*,*)' ERROR! ', jump1,jump2, prevsec,postsec,finsec 
         ENDIF

! ================================================================
! ================ SYMMETRIC LIBRATION  ==========================
! ================================================================      
      ELSEIF(finsec.eq.prevsec)THEN 
         write(*,*)'symmetric libration' 
         IF(jump1+jump2.ne.0)THEN 
            write(*,*)' ERROR! ', jump1,jump2, prevsec,postsec,finsec 
         ENDIF

! ================================================================
! ================ ASYMMETRIC LIBRATION  =========================
! ================================================================      
! fermarsi un attimo prima della giusta fine                        
! write(*,*)'librazione asimmetrica'                                
! ===============================================================
      ELSE 
! control                                                           
         write(*,*)'qui non ci dovremmo capitare mai!!' 
                                                                        
      ENDIF
                                                                        
      RETURN 
    END SUBROUTINE prstop
