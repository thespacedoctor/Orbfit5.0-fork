clear;close all;
load VAs;
[s1 s2]=size(VAs);
t=VAs(2,2)/365.25+1858.87953;
x(:,1:3)=VAs(2:s1,3:5);
dmin=VAs(2:s1,9);
dmind=VAs(2:s1,10);
%xp=x(:,1);
%yp=x(:,2);
%zp=x(:,3);
nom=VAs(1,1); % label of nominal orbit
n1=VAs(2,1);  % label of first VA
n2=VAs(s1,1); % label of last VA
clear VAs;

figure(1)
hold on;
tt=-5000:1:5000;
lt=length(tt);

for j=1:n2-n1+1;
    dminline= dmin(j) + dmind(j)*tt;
    if j==nom-n1+1;
        j
    plot(t+tt,dminline,'b','LineWidth',1.5); 
    else
    plot(t+tt,dminline,'k:'); 
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
line([tcmin tcmax],[0 0],'LineWidth',2,'Color','Red')

xlabel('time (yr)')
ylabel('orbit distance (AU)')
%title('crossing time estimate for 1979 XB')
text(tcmin-3,-0.001,'t_1');
text(tcmax+3,-0.001,'t_2');
