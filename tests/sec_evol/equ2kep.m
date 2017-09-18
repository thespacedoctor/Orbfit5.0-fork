clear;
load equfile;
a=equfile(:,1);
h=equfile(:,2);
k=equfile(:,3);      
p=equfile(:,4);
q=equfile(:,5);
clear equfile;
[s1 s2]=size(a);

for i=1:s1;
ecc(i)=sqrt(h(i)^2+k(i)^2);
tanImez(i)=sqrt(p(i)^2+q(i)^2); 
inc(i)=2*atan(tanImez(i));
omtil(i)=atan2(h(i)/ecc(i),k(i)/ecc(i));
Omnod(i)=atan2(p(i)/tanImez(i),q(i)/tanImez(i));
omega(i)=omtil(i)-Omnod(i);
end;
omdeg=omega*180/pi + 360;