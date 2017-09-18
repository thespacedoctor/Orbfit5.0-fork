load singeq.mev
t=singeq(:,1);om=singeq(:,2);omnod=singeq(:,3);
e=singeq(:,4);I=singeq(:,5);nit=singeq(:,6);
dn=singeq(:,7);nnod=singeq(:,8);
clear singeq

load mevfile.fla
ts=mevfile(:,1);oms=mevfile(:,2);omnods=mevfile(:,3);
es=mevfile(:,4);Is=mevfile(:,5);nits=mevfile(:,6);
dns=mevfile(:,7);nnods=mevfile(:,8);
clear mevfile.fla

s1=length(ts);
for i=1:7:s1
    omse(i)=oms(i);ese(i)=es(i);	       
end

s2=length(t);
for i=1:2:s2
    omc(i)=om(i);ec(i)=e(i);    
end