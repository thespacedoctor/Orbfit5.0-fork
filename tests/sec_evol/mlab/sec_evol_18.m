clear;close all;
load astevol;

load vast1;
vastfil=vast1;
for i=2:18
fil1=strcat('vast',int2str(i));
fil=load(fil1);
vastfil = vastfil + fil;
end;
vastfil=vastfil/18;

t=astevol(:,1);

torb9=vastfil(:,1);
len=length(torb9);

step=10;

tfil=torb9(1:step:len);
omega=unwrap(vastfil(:,2));
omegafil=omega(1:step:len);
Omnod=unwrap(vastfil(:,3));
Omnodfil=Omnod(1:step:len);
ecc=vastfil(:,4);
eccfil=ecc(1:step:len);
inc=unwrap(vastfil(:,5));
incfil=inc(1:step:len);


figure(1);
subplot(2,2,1)
plot(t,astevol(:,2));
hold on;
plot(tfil,omegafil,'r');
xlabel('time');ylabel('omega');

subplot(2,2,2)
plot(t,astevol(:,3));
hold on;
plot(tfil,Omnodfil,'r');
xlabel('time');ylabel('Omnod');

subplot(2,2,3)
plot(t,astevol(:,4));
hold on;
plot(tfil,eccfil,'r');
xlabel('time');ylabel('ecc');

subplot(2,2,4)
plot(t,astevol(:,5));
hold on;
plot(tfil,incfil,'r');
xlabel('time');ylabel('incl');


figure(2);
plot(astevol(:,2),astevol(:,4));
hold on;
plot(omegafil,eccfil,'r');
xlabel('omega');ylabel('ecc');

clear astevol;clear vastfil;
