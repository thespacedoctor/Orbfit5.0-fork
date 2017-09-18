% compare secular evolution with circular secular evolution
clear;
close all;
load ast_circ.fla;
[x1 y1]=size(ast_circ);
tdat0=ast_circ(1,1);
t=1858.87953 + tdat0/365.25 + ast_circ(2:x1,1);
omega=ast_circ(2:x1,2);
Omnod=ast_circ(2:x1,3);
ecc=ast_circ(2:x1,4);
inc=ast_circ(2:x1,5);


load mevfile.fla;
[x2 y2]=size(mevfile);
tast0=mevfile(1,1);
for i=1:7:x2-1
ts(i)=1858.87953 + tast0/365.25 + mevfile(i+1,1);
oms(i)=mevfile(i+1,2);
onods(i)=mevfile(i+1,3);
es(i)=mevfile(i+1,4);
is(i)=mevfile(i+1,5);
end

tmax=50000;

figure(1);

subplot(222)
plot(t(t<tmax-11),omega(t<tmax-11),'+b');
hold on;
plot(ts(ts<tmax+11),oms(ts<tmax+11),'xr');
xlabel('tempo');
ylabel('omega');
axis([0 20000 265 295]);

subplot(224)
plot(t(t<tmax-20),Omnod(t<tmax-20),'+b');
hold on;
plot(ts(ts<tmax+11),onods(ts<tmax+11),'xr');
xlabel('tempo');
ylabel('Omega nodale');
axis([0 20000 240 365]);

subplot(221)
plot(t(t<tmax-11),ecc(t<tmax-11),'+b');
hold on;
plot(ts(ts<tmax+11),es(ts<tmax+11),'xr');
xlabel('tempo')
ylabel('ecc');
axis([0 20000 0.3 0.7]);

subplot(223)
plot(t(t<tmax-20),inc(t<tmax-20),'+b');
hold on;
plot(ts(ts<tmax+11),is(ts<tmax+11),'xr');
xlabel('tempo');
ylabel('incl');
axis([0 20000 38 55]);

suplabel('Asteroide (1981) Midas','t');

%print -depsc 1981_secevol.eps
print -deps 1981_secevol_bw.eps

