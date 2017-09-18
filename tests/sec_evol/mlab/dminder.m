clear;close all;
load VAs;
t=VAs(1,2)/365.25+1858.87953;
x(:,1:3)=VAs(:,3:5);
dmin=VAs(:,9);
dmind=VAs(:,10);
dmin_2=VAs(:,12);
dmind_2=VAs(:,13);
[s1 s2]=size(dmin);
xp=x(:,1);
yp=x(:,2);
zp=x(:,3);

figure(1)
hold on;
tt=-500:1:500;
lt=length(tt);

for j=1:s1;
    dminline= dmin(j) + dmind(j)*tt;
    plot(t+tt,dminline,'k:'); 
    tc(j) = -dmin(j)/dmind(j) + t;
    %    plot(tc(j),0,'*');
    plot(t,dmin(j),'*b');
    lab=num2str(j);
    %    text(t+tt(lt),dminline(lt),lab);
    %text(tc(j),0,lab);
end;
t0=min(t+tt(abs(dminline)<0.05));
t1=max(t+tt(abs(dminline)<0.05));
%axis([t0 t1 -0.05 0.05]);
tcmin=min(tc);
tcmax=max(tc);
delta=tcmax-tcmin
line([tcmin tcmax],[0 0],'LineWidth',2,'Color','Red')


% second minimum
for j=1:s1;
    dminline_2= dmin_2(j) + dmind_2(j)*tt;
    plot(t+tt,dminline_2,'k:'); 
    tc(j) = -dmin_2(j)/dmind_2(j) + t;
    %    plot(tc(j),0,'*');
    plot(t,dmin_2(j),'*b');
    lab=num2str(j);
    %    text(t+tt(lt),dminline(lt),lab);
    %text(tc(j),0,lab);
end;
t0=min(t+tt(abs(dminline_2)<0.05));
t1=max(t+tt(abs(dminline_2)<0.05));
%axis([t0 t1 -0.05 0.05]);
tcmin=min(tc);
tcmax=max(tc);
delta=tcmax-tcmin
line([tcmin tcmax],[0 0],'LineWidth',2,'Color','Red')