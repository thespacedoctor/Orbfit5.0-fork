
rx=max(rI);rm=min(rI);nnr=size(rI);nr=nnr(2)
rIp=rI;dr=(rx-rm)/(nr-1);rI=rx:-dr:rm;
V4=[-.5 0 .5]
V2=[-2 -1 0 1 2]
V1=[-1 0 1]
V5=[-5 -4 -3 -2 -1 0 1 2 3 4 5]
figure(1);
hold off
d10=gq+sq-g6-s6;
contour(ra,rI,d10,V4,'LineWidth',2,'Color',[.6 0 0])
pause
hold on
gtext('z1')
hold on
d13=gq-2*g6+g5;
contour(ra,rI,d13,V4,'LineWidth',2,'Color',[.6 0 0])
pause
gtext('g-2*g6+g5')
d46=sq-s6-g5+g6;
contour(ra,rI,d46,V4','LineWidth',2,'Color',[.6 0 0])
pause
gtext('s-s6-g5+g6')
d2=gq-g6;
contour(ra,rI,d2,V2,'LineWidth',2,'Color',[.6 0 0])
pause
gtext('g-g6')
d1=gq-g5;
contour(ra,rI,d1,V2,'LineWidth',2,'Color',[.6 0 0])
pause
gtext('g-g5')
d4=sq-s6;
contour(ra,rI,d4,V2,'LineWidth',2,'Color',[.6 0 0])
pause
gtext('s-s6')
z2=(gq-g6)*2+(sq-s6);
contour(ra,rI,z2,V4,'LineWidth',2,'Color',[.6 0 0])
pause
gtext('z2')
z3=(gq-g6)*3+(sq-s6);
contour(ra,rI,z3,V4,'LineWidth',2,'Color',[.6 0 0])
pause
gtext('z3')
d17=(gq-sq-g5+s6);
contour(ra,rI,d17,V4,'LineWidth',2,'Color',[.6 0 0])
pause
gtext('g-s-g5+s6')
pause

title('Secular resonances and proper elements vers. 9.0, e=0.075')
xlabel('Proper a (AU)')
ylabel('Proper inclination (DEG)')


load numb.syn;
N=numb(:,1);a=numb(:,3);e=numb(:,4);sI=numb(:,5);
g=numb(:,7);s=numb(:,8);
LCE=numb(:,9);H=numb(:,2);I=asin(sI)*180/pi;
clear numb;
load numb.sig;
sa=numb(:,2);se=numb(:,3);ssI=numb(:,4);
clear numb;
h=H<15&e>0.05&e<0.1;
f=LCE>50&h;
fs=(se>0.005|ssI>0.005)&h;
%plot(a(h),I(h),'.');
save secres.mat



plot(a(fs),I(fs),'.k');
%plot(a(f),I(f),'.r');


figure(2);
hold off
plot(g(fs),s(fs),'.k')
hold on;

%plot(g(f),s(f),'.r')

