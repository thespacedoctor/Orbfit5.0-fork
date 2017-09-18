clear;close all;

load mevfile;
[S1 S2]=size(mevfile);
t0=mevfile(1,1);
a=mevfile(1,2);
time=mevfile(2:S1,1);
omega=mevfile(2:S1,2);
Omnod=mevfile(2:S1,3);
ecc=mevfile(2:S1,4);
incl=mevfile(2:S1,5);

nstep=S1-1;


% ================================
figure(1);
mov = avifile('secevol.avi');

u = 0:0.1:2*pi+0.1;
nn = size(u);
nu = nn(2);
up = u;

for i=1:3;
pause(1);

%nstep;
 e=ecc(i);
 om=omega(i);
 Omn=Omnod(i);
 inc=incl(i);

%%%%%%%%%%% ASTEROID ORBIT %%%%%%%%%%%
 x1 = a*(cos(u)-e);
 y1 = a*sqrt(1-e^2)*sin(u);
 z1 = 0;

 x3 = x1*cos(om)-y1*sin(om);
 y3 = (x1*sin(om)+y1*cos(om))*cos(inc);
 z3 = (x1*sin(om)+y1*cos(om))*sin(inc);

 x = x3*cos(Omn)-y3*sin(Omn);
 y = x3*sin(Omn)+y3*cos(Omn);
 z = z3;


%%%%%%%%%%% EARTH ORBIT %%%%%%%%%%%%%%%%
ap=1.00020368;
ep=0.01644844;
incp=0.00043074;
Omnp=154.50192116;
omp=309.71597223;

xp1 = ap*(cos(up)-ep);
yp1 = ap*sqrt(1-ep^2)*sin(up);
zp1 = 0;

xp3 = xp1*cos(omp)-yp1*sin(omp);
yp3 = (xp1*sin(omp)+yp1*cos(omp))*cos(incp);
zp3 = (xp1*sin(omp)+yp1*cos(omp))*sin(incp);

xp = xp3*cos(Omnp)-yp3*sin(Omnp);
yp = xp3*sin(Omnp)+yp3*cos(Omnp);
zp = zp3;

% axis equal;
 set(gca,'DataAspectRatio',[1 1 1],'Color','white');


hold off;
h=plot3(x,y,z);
set(h,'LineWidth',2);

hold on;
plot3(xp,yp,zp);
 
%camlight left;
%lighting phong;
Frame = getframe(gca);
mov = addframe(mov,Frame);

%clear x1; clear y1; clear z1; clear x3; clear y3; clear z3; clear x; clear y; clear z;
%clear xp1,yp1,zp1,xp3,yp3,zp3,xp,yp,zp;

end;
%mov = close(mov);
