load clevmedia.out
[k1,k2]=size(clevmedia)
gr=clevmedia(3:k1,:);
ap(2)=clevmedia(1,1)
ap(3)=clevmedia(1,2)
ap(4)=clevmedia(1,3)
ap(5)=clevmedia(1,4)
ap(6)=clevmedia(1,5)
a=clevmedia(2,1)
emax=clevmedia(2,2)
Imax=clevmedia(2,3)*180/pi
ngro=clevmedia(2,4)
ngre=clevmedia(2,5)

j=gr(:,1);
i=gr(:,2);
om=gr(:,3);
e=gr(:,4);
r=gr(:,5);
aber=gr(:,6);
%ninterv=gr(:,7);
nev=gr(:,7);
clear clevmedia gr
[i1 i2]=size(e)
for ii=1:i1
   iii=i(ii)+1;
   jjj=j(ii)+1;
   rq(iii,jjj)=r(ii);
   rq(1+sqrt(i1)-iii,sqrt(i1)+jjj)=r(i1-ii+1);
%  eq(iii,jjj)=e(ii);
%  oq(iii,jjj)=om(ii)*180/pi;
   ev(iii)=e(ii);
   ov(jjj)=om(ii)*180/pi;
   ov1 = [ov,90+ov];
end

rql = [rq,rq];
evl = ev;
ovl = [ov1,180+ov1];

%figure(1)
%hold off
%meshc(ov1,ev,rq)
%mesh(ovl,evl,rql)
%surf(ovl,evl,rql)
%colormap(cool)
%shading interp
%xlabel('omega')
%ylabel('eccentricity')
%tit=['a=',num2str(a,4),' Imax=',num2str(Imax,4)]
%title(tit)

figure(1)
hold off
%contour(ov1,ev,rq,30)
 contour(ovl,evl,rql,30)
%colormap(cool)
%shading interp
xlabel('omega')
ylabel('eccentricity')
tit=['a=',num2str(a,4),' Imax=',num2str(Imax,4)]
title(tit)
for i=2:6
  ecosom=(1-ev.^2)*a/ap(i)-1;
  ee=ecosom./ev;
  f=(ee<1)&(ee>-1);
  ec=ev(f);
  oc=acos(abs(ee(f))).*(180/pi);
  sc=sum(f)
  eee(i)=(abs(ap(i) - a)/a);
  ec1 = [eee(i),ec];
  oc1 = [0,oc];
  ec2 = [eee(i),ec];
  oc2 = 180-oc1;
  ec3 = [eee(i),ec];
  oc3 = [180,oc+180];
  ec4 = [eee(i),ec];
  oc4 = 360-oc1;

  sc=sum(f)
  if sc>0
    hold on
    plot(oc1,ec1)
    plot(oc2,ec2)
    plot(oc3,ec3)
    plot(oc4,ec4)

    hold off
  end
%%%%%%%% LABELS FOR THE PLANETS %%%%%%%%%%%%
  if i == 2 
     text(oc3(1),ec3(1)-0.005,'V')
  elseif i == 3 
     text(oc3(1),ec3(1)-0.005,'E')
  elseif i == 4 
     text(oc3(1),ec3(1)-0.005,'M')
  elseif i == 5 
     text(oc3(1),ec3(1)-0.005,'J')
  elseif i == 6 
     text(oc3(1),ec3(1)-0.005,'S')
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end
