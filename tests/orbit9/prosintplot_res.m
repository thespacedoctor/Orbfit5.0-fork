clear



% input data
sou=input('source? 1=source files 0=.mat file')
if sou==1
% load proper elements
   load v5core.syn
   pro=v5core;
% problem for N with multiopposition
N=pro(:,1);a=pro(:,2);e=pro(:,3);sinI=pro(:,4);  %H=pro(:,2);
n=pro(:,5);g=pro(:,6);s=pro(:,7);LCE=pro(:,8);  %My=pro(:,11);
   comment='resonance g+g5-2*g6aa; adapted synthetic proper elements'
   np=length(a)
   clear v5core pro
% inclination
   I=asin(sinI)*(180/pi);

% load uncertainties
   load v5core.sig
   sig=v5core;
% problem for N1 with multiopposition
   N1=sig(:,1);sa=sig(:,2);se=sig(:,3);ssI=sig(:,4);
   sn=sig(:,5);sg=sig(:,6);ss=sig(:,7);lamfit=sig(:,8);
   inconsistency=sum(N~=N1)
   np1=length(N1)
   if np~=np1
      comment=' inconsistenza tra syn e sig'
   end
   clear sig v5core

% frequencies for secular resonances
   g5=4.25749319;g6=28.24552984;g7=3.08675577;g8=0.67255084;
   s6=-26.34496354;s7=-2.99266093;s8=-0.69251386;
% from laskar 2004
   g3=17.368;g4=17.916;s3=-18.850;s4=-17.755;

% definition of chaotic, resonant
   bad=sa>3e-4|LCE>50;
   res=abs(g+g5-2*g6)<0.065;
   nres=sum(res)
   nonres=sum(~res)

% load initial phases
   load v5core.ang
   ang=v5core;
   N2=ang(:,1);peri=ang(:,3);libph=ang(:,4);sigph=ang(:,8);
   clear ang v5core


% save
   save synt
else
   load synt
end
%

fig=input('global figures? 1 yes 0 no')
if fig

%   a sin(I) plot
   figure(1)
   emin=input('minimum e for asinI plot?')
   emax=input('maximum e for asinI plot?')
   ff=e>emin&e<emax;sum(ff)
   hold off
   plot(a(ff&~res),sinI(ff&~res),'+k','LineWidth',3)
   hold on
   plot(a(ff&res),sinI(ff&res),'ok','LineWidth',3)
   xlabel('Proper semiamjor axis (au)')
   ylabel('Proper sine of inclination')
   title('Region of secular resonance g+g5-2*g6: o=resonant')
   %print -depsc fam_asinI.eps
   %print -djpeg MB_plot_aI.jpg

% a e plot
   figure(2)
   Imin=input('min sinI for ae plot?')
   Imax=input('max sinI for ae plot?')
   f=sinI>Imin&sinI<Imax;sum(f)
   hold off
   plot(a(f&~res),e(f&~res),'+k','LineWidth',3)
   hold on
   plot(a(f&res),e(f&res),'ok','LineWidth',3)
   xlabel('Proper semiamjor axis (au)')
   ylabel('Proper eccentricity')
   title('Region of secular resonance g+g5-2*g6: o=resonant')
   %print -depsc fam_ae.eps
   %print -djpeg MB_plot

% e sin(i) e plot
   figure(3)
   amin=input('min a for Ie plot?')
   amax=input('max a for Ie plot?')
   fa=a>amin&a<amax;
   hold off
   plot(sinI(fa&~res),e(fa&~res),'+k','LineWidth',3)
   hold on
   plot(sinI(fa&res),e(fa&res),'ok','LineWidth',3)
   xlabel('Proper sine of inclination')
   ylabel('Proper eccentricity')
   title('Region of secular resonance g+g5-2*g6: o=resonant')
   %print -depsc fam_Ie.eps
end

% resonance libration plot
  k=(e(res)-1).*cos(libph(res));
  h=(e(res)-1).*sin(libph(res));
  figure(4)
  hold off
  plot(h,k,'.r','LineWidth',3)
  xlabel('k=e*sin(phi)')
  ylabel('h=e*sin(phi)')
  title('Region of secular resonance g+g5-2*g6: polar plot')
