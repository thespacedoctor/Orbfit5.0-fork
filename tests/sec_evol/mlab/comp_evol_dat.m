% compare secular evolution with pure numerical integration
clear;
close all;
load vastdat.fla;
[x1 y1]=size(vastdat);
tdat0=vastdat(1,1);

%t=1858.87953 + tdat0/365.25 + vastdat(2:decim:x1,1);
%omega=vastdat(2:decim:x1,2);
%Omnod=vastdat(2:decim:x1,3);
t=1858.87953 + tdat0/365.25 + vastdat(2:x1,1);
omega=vastdat(2:x1,2);
Omnod=vastdat(2:x1,3);

%for j=1:length(Omnod);
%    if Omnod(j) > 180;
%        Omnod(j) = Omnod(j)-360;
%    end;
%end;

%ecc=vastdat(2:decim:x1,4);
%inc=vastdat(2:decim:x1,5);
ecc=vastdat(2:x1,4);
inc=vastdat(2:x1,5);

decim=input('decimation factor (integer)');

load secfile.fla;
[x2 y2]=size(secfile);
tast0=secfile(1,1);
%ts=1858.87953 + tast0/365.25 + secfile(2:x2,1);
%oms=secfile(2:x2,2);
%onods=secfile(2:x2,3);
%es=secfile(2:x2,4);
%is=secfile(2:x2,5);

ts=1858.87953 + tast0/365.25 + secfile(2:decim:x2,1);
oms=secfile(2:decim:x2,2);
onods=secfile(2:decim:x2,3);
es=secfile(2:decim:x2,4);
is=secfile(2:decim:x2,5);

tmax=8000;

figure(1);

subplot(221)
plot(t(t<tmax-11),ecc(t<tmax-11),'k');
hold on;
plot(ts(ts<tmax+11),es(ts<tmax+11),'k+');
xlabel('t (yr)')
ylabel('e');
axis tight;

subplot(223)
plot(t(t<tmax-20),inc(t<tmax-20),'k');
hold on;
plot(ts(ts<tmax+11),is(ts<tmax+11),'k+');
xlabel('t (yr)');
ylabel('I (deg)');
axis tight;

subplot(222)
plot(t(t<tmax-11),omega(t<tmax-11),'k');
hold on;
plot(ts(ts<tmax+11),oms(ts<tmax+11),'k+');
xlabel('t (yr)');
ylabel('\omega (deg)');
axis tight;

subplot(224)
plot(t(t<tmax-20),Omnod(t<tmax-20),'k');
hold on;
plot(ts(ts<tmax+11),onods(ts<tmax+11),'k+');
xlabel('t (yr)');
ylabel('\Omega (deg)');
axis tight;

print -deps secevol_bw.eps

figure(2);
plot(omega(t<tmax-11),ecc(t<tmax-11),'k');
hold on;
plot(oms(ts<tmax+11),es(ts<tmax+11),'k+');
xlabel('\omega (deg)');
ylabel('e');
axis tight;

print -deps omega_e_bw.eps

