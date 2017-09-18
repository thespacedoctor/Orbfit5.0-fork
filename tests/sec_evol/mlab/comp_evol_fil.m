% compare secular evolution with pure numerical integration
clear;
close all;
load vastfil.fla;
[x1 y1]=size(vastfil);
tfil0=vastfil(1,1);
decim=input('decimation factor (integer)');
t=1858.87953 + tfil0/365.25 + vastfil(2:decim:x1,1);
omega=vastfil(2:decim:x1,2);
Omnod=vastfil(2:decim:x1,3);
ecc=vastfil(2:decim:x1,4);
inc=vastfil(2:decim:x1,5);

load secfile.fla;
[x2 y2]=size(secfile);
tast0=secfile(1,1);
ts=1858.87953 + tast0/365.25 + secfile(2:x2,1);
oms=secfile(2:x2,2);
onods=secfile(2:x2,3);
es=secfile(2:x2,4);
is=secfile(2:x2,5);

tmax=8000;

figure(1);

subplot(221)
plot(t(t<tmax-11),ecc(t<tmax-11),'b');
hold on;
plot(ts(ts<tmax+11),es(ts<tmax+11),'r');
xlabel('tempo')
ylabel('ecc');
%axis([2000 8000 0.222 0.227]);

subplot(223)
plot(t(t<tmax-20),inc(t<tmax-20),'b');
hold on;
plot(ts(ts<tmax+11),is(ts<tmax+11),'r');
xlabel('tempo');
ylabel('inc');
%axis([2010 8000 10.2 10.9]);

subplot(222)
plot(t(t<tmax-11),omega(t<tmax-11),'b');
hold on;
plot(ts(ts<tmax+11),oms(ts<tmax+11),'r');
xlabel('tempo');
ylabel('omega');
%axis([2010 8000 180 250]);

subplot(224)
plot(t(t<tmax-20),Omnod(t<tmax-20),'b');
hold on;
plot(ts(ts<tmax+11),onods(ts<tmax+11),'r');
xlabel('tempo');
ylabel('Omega nodale');
%axis([2010 8000 265 305]);

%suplabel('Asteroide (433) Eros','t');

%print -depsc evol_c.eps
%print -deps evol_bw.eps

figure(2);
plot(omega(t<tmax-11),ecc(t<tmax-11),'b');
hold on;
plot(oms(ts<tmax+11),es(ts<tmax+11),'r');
xlabel('argomento del pericentro (\omega)');
ylabel('eccentricity');
%axis([180 250 0.222 0.227]);
%title('Asteroide (433) Eros')

%print -depsc evol_omecc_c.eps
%print -deps evol_omecc_bw.eps

