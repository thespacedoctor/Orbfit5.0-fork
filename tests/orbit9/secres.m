clear
datam
ngrid=data(1)
ra=data(2):data(3):(data(2)+data(3)*(ngrid-1));
rI=data(4):data(5):(data(4)+data(5)*(ngrid-1));
re=data(6):data(7):(data(6)+data(7)*(ngrid-1));
w=data(8)
om=data(9)
load gs.fla
load forced.fla
for i=1:ngrid;
   for j=1:ngrid;
     gq(i,j)=gs(1-i+j*ngrid,1);
     sq(i,j)=gs(1-i+j*ngrid,2);
     xi(i,j)=forced(1-i+j*ngrid,1);
     et(i,j)=forced(1-i+j*ngrid,2);
   end;
end;
clear gs i j
g5=4.25749319
g6=28.24552984
g7=3.08675577
g8=0.67255084
s6=-26.34496354
s7=-2.99266093
s8=-0.69251386

%n=data(1);i=0:n-1;
%ra=data(2)+i*data(3);
%rI=data(4)+i*data(5);
