% compare secular evolution with circular secular evolution
clear;close;

load mevfile.fla;
[x2 y2]=size(mevfile);
tast0=mevfile(1,1);
ts=1858.87953 + tast0/365.25 + mevfile(2:x2,1);
oms=mevfile(2:x2,2);
onods=mevfile(2:x2,3);
es=mevfile(2:x2,4);
is=mevfile(2:x2,5);

load mev_circ.fla;
[x1 y1]=size(mev_circ);
tcirc0=mev_circ(1,1);
step=2;
tc=1858.87953 + tcirc0/365.25 + mev_circ(2:step:x1,1);
omega=mev_circ(2:step:x1,2);
Omnod=mev_circ(2:step:x1,3);
ecc=mev_circ(2:step:x1,4);
inc=mev_circ(2:step:x1,5);

figure(1);
subplot(2,2,1)
plot(ts,es,'+r');%,'LineWidth',2);
hold on;
plot(tc,ecc,'xb');
xlabel('tempo')
ylabel('ecc');
%axis([0 30000 0.3 0.36]);

subplot(2,2,2)
plot(ts,oms,'+r');%,'LineWidth',2);
hold on;
plot(tc,omega,'xb');
xlabel('tempo');
ylabel('omega');
%axis([0 30000 200 500]);

subplot(2,2,3)
plot(ts,is,'+r');%,'LineWidth',2);
hold on;
plot(tc,inc,'xb');
xlabel('tempo');
ylabel('incl');
%axis([0 30000 12 15]);

subplot(2,2,4)
plot(ts,onods,'+r');%,'LineWidth',2);
hold on;
plot(tc,Omnod,'xb');
xlabel('tempo');
ylabel('Omega nodale');
%axis([0 30000 150 400]);

%suplabel('Asteroide (1620) Geographos','t');

figure(2)
subplot(221)
plot(ts,es,'r'); %'LineWidth',2);
xlabel('tempo (anni)');
ylabel('ecc');
hold on;
subplot(222)
plot(ts,is,'r'); %'LineWidth',2);
xlabel('tempo (anni)');
ylabel('inc');

%suplabel('Asteoride (433) Eros','t')

%print -depsc 433_secevol_omecc2.eps
%print -deps 433_secevol_bw_omecc2.eps

