clear;close;

load elemevol;
t=elemevol(:,1);
inc=elemevol(:,2);
omnod=elemevol(:,3);
omega=elemevol(:,4);

%load incl
%t=incl(:,1);
%inc=incl(:,2);
figure(3)
plot(t,inc)

%clear;
%load node
%t=node(:,1);
%omnod=node(:,2);
figure(4)
plot(t,omnod)

%clear;
%load peric
%t=peric(:,1);
%omega=peric(:,2);
figure(5)
plot(t,omega)
