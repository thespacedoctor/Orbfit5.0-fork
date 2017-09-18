clear;close all;

load VAs;
[s1 s2]=size(VAs);
t=VAs(2,2)/365.25+1858.87953;
x(:,1:3)=VAs(2:s1,3:5);
dmin=VAs(2:s1,9);
dmind=VAs(2:s1,10);
xp=x(:,1);
yp=x(:,2);
zp=x(:,3);
nom=VAs(1,1); % label of nominal orbit
n1=VAs(2,1);  % label of first VA
n2=VAs(s1,1); % label of last VA
clear VAs;

figure(1);
hold on;

for h=n1:n2;
    h
    filnam=strcat('dmin_',num2str(h))
    filfla=strcat(filnam,'.fla')
    S=load(filfla);
    [s1 s2]=size(S)
    t0=S(1,1);
    tevol=S(2:s1,1)+t0/365.25+1858.87953;
    dmin1=S(2:s1,3);
    deriv1=S(2:s1,5);
    dmin2=S(2:s1,4);
    deriv2=S(2:s1,6);
    if h==nom
        plot(tevol,dmin1,'b','LineWidth',2);
    else
        plot(tevol,dmin1,'g');
    end;
    %   plot(tevol,dmin2,'r');
end;


figure(1)
hold on;
tt=0:0.2:250;
lt=length(tt);

for j=1:n2-n1+1;
    dminline= dmin(j) + dmind(j)*tt;
    if j==nom-n1+1;
    plot(t+tt,dminline,'b','LineWidth',2,'LineStyle','--'); 
    else
        %    plot(t+tt,dminline,'k:'); 
    end
    tc(j) = -dmin(j)/dmind(j) + t;
    %    plot(tc(j),0,'*');
    plot(t,dmin(j),'+k');
    lab=num2str(j);
    %    text(t+tt(lt),dminline(lt),lab);
    %text(tc(j),0,lab);
end;
t0=min(t+tt(abs(dminline)<0.05));
t1=max(t+tt(abs(dminline)<0.05));
%axis([t0 t1 -0.05 0.05]);
axis([1970 2200 -0.003 0.03]);
tcmin=min(tc);
tcmax=max(tc);
delta=tcmax-tcmin
line([tcmin tcmax],[0 0],'LineWidth',3,'Color','Red')

xlabel('time (yr)')
ylabel('orbit distance (AU)')
%title('crossing time estimate for 1979 XB')
%text(tcmin-2,-0.001,'t1');
%text(tcmax+2,-0.001,'t2');
text(tcmin-1,0,'(');
text(tcmax-1,0,')');

% dmin zero line
%line([1980 2200], [0 0],'Color','Black')

tau1=2084.5;
tau2=2168.3;
text(tau1-1,0,'[');
text(tau2,0,']');
eps=0.00015
line([tau1 tau2], [eps eps],'Color',[0.2,0.5,0.5],'LineWidth',2)
%text(tau1-1,-0.001,'tau1');
%text(tau2+1,-0.001,'tau2');

print -depsc confronto_1979XB.eps;

