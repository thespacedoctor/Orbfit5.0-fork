clear;close all;

figure(1);
hold on;

for h=1:7;
   h
   filnam=strcat('dmin_',num2str(h))
   filfla=strcat(filnam,'.fla')
   S=load(filfla);
   [s1 s2]=size(S)
   t0=S(1,1);
   t=S(2:s1,1)+t0/365.25+1858.87953;
   dmin1=S(2:s1,3);
   deriv1=S(2:s1,5);
   dmin2=S(2:s1,4);
   deriv2=S(2:s1,6);
   plot(t,dmin1,'.r');
   %   plot(t,dmin2,'r');
end;

clear;
load VAs;
t=VAs(1,2)/365.25+1858.87953;
x(:,1:3)=VAs(:,3:5);
dmin=VAs(:,9);
dmind=VAs(:,10);
[s1 s2]=size(dmin);
xp=x(:,1);
yp=x(:,2);
zp=x(:,3);

figure(1)
%hold on;
tt=-20:0.2:250;
lt=length(tt);

for j=1:s1;
        dminline= dmin(j) + dmind(j)*tt;
        if j==1
            plot(t+tt,dminline,'b'); 
            tc(j) = -dmin(j)/dmind(j) + t;
            plot(tc(j),0,'*');
            plot(t,dmin(j),'o');
        elseif j==2
            plot(t+tt,dminline,'b'); 
            tc(j) = -dmin(j)/dmind(j) + t;
            plot(tc(j),0,'*');
            plot(t,dmin(j),'o');
        elseif j==3
            plot(t+tt,dminline,'b'); 
            tc(j) = -dmin(j)/dmind(j) + t;
            plot(tc(j),0,'*');
            plot(t,dmin(j),'o');
       elseif j==4
            plot(t+tt,dminline,'k','LineWidth',2); 

            tc(j) = -dmin(j)/dmind(j) + t;
            plot(tc(j),0,'*');
            plot(t,dmin(j),'o');
        elseif j==5
            plot(t+tt,dminline,'b'); 
            tc(j) = -dmin(j)/dmind(j) + t;
            plot(tc(j),0,'*');
            plot(t,dmin(j),'o');
        elseif j==6
            plot(t+tt,dminline,'b'); 
            tc(j) = -dmin(j)/dmind(j) + t;
            plot(tc(j),0,'*');
            plot(t,dmin(j),'o');
        elseif j==7
            plot(t+tt,dminline,'b'); 
            tc(j) = -dmin(j)/dmind(j) + t;
            plot(tc(j),0,'*');
            plot(t,dmin(j),'o');
        end
        lab=num2str(j);
        text(t+tt(lt),dminline(lt),lab);
        text(tc(j),0,lab);

end;
t0=min(t+tt(abs(dminline)<0.05));
t1=max(t+tt(abs(dminline)<0.05));
axis([t0 t1 -0.05 0.05]);
min(tc)
max(tc)
delta=max(tc)-min(tc)