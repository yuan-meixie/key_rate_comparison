function [objv] = qkd_figure(X,L_set,datapara)
% main function of decoy-state QKD
l=L_set; 
  
alfa=datapara(1);
eta_int = datapara(2);
N =datapara(3);
pdz= datapara(4);
pdx= datapara(5);
etabobz= datapara(6);
etabobx=datapara(7);
 
edx=0.02;  %loss-independent misalignment error in the X basis
edz=0.02; %loss-independent misalignment error in the Z basis
e0=1/2;

 
ebs=1e-10;
ebscor=1e-10;
ebssec=1e-10;

fec=1.1;
tab=10^(-alfa*l/10)*(10^(-eta_int/10));
etax=tab*etabobx;
etaz=tab*etabobz;

qx=X(:,1)/1000;
w=X(:,2)/1000;
v=X(:,3)/1000;
mu=X(:,4)/1000;
pmu=X(:,5)/1000;
pnu=X(:,6)/1000;
pw=X(:,7)/1000;
p0=1-pmu-pnu-pw;
qz=1-qx;


Q0z=1/2*(1+(1-pdx)^2)*(1-(1-pdz)^2);
Q0x=1/2*(1+(1-pdz)^2)*(1-(1-pdx)^2);
Qwx=1/2*(1-(1-pdx)^2*exp(-w*qx*etax))*(1+(1-pdz)^2*exp(-w*qz*etaz));
Qvx=1/2*(1-(1-pdx)^2*exp(-v*qx*etax))*(1+(1-pdz)^2*exp(-v*qz*etaz));
Qvz=1/2*(1-(1-pdz)^2*exp(-v*qz*etaz))*(1+(1-pdx)^2*exp(-v*qx*etax));
Qmux=1/2*(1-(1-pdx)^2*exp(-mu*qx*etax))*(1+(1-pdz)^2*exp(-mu*qz*etaz));
Qmuz=1/2*(1-(1-pdz)^2*exp(-mu*qz*etaz))*(1+(1-pdx)^2*exp(-mu*qx*etax));

EmuzQmuz=edz*Qmuz+(e0-edz)*1/2*(1-(1-pdz)^2)*exp(-mu*qz*etaz)*(1+(1-pdx)^2*exp(-mu*qx*etax));
EvzQvz=edz*Qvz+(e0-edz)*1/2*(1-(1-pdz)^2)*exp(-v*qz*etaz)*(1+(1-pdx)^2*exp(-v*qx*etax));
EwxQwx=edx*Qwx+(e0-edx)*1/2*(1-(1-pdx)^2)*exp(-w*qx*etax)*(1+(1-pdz)^2*exp(-w*qz*etaz));

n0z=N*p0*Q0z;
n0x=N*p0*Q0x;
nnuz=N*pnu*Qvz;
nnux=N*pnu*Qvx;
nmuz=N*pmu*Qmuz;
nmux=N*pmu*Qmux;
nz=nnuz+nmuz;

m0x=N*p0*e0*Q0x;
mwx=N*pw*EwxQwx;
mmuz=N*pmu*EmuzQmuz;
mvz=N*pnu*EvzQvz;
mz=mvz+mmuz;
 
[En0zl,En0zu] = xtoxbzmg(n0z,ebs);
[~,En0xu] = xtoxbzmg(n0x,ebs);
[Envzl,~] = xtoxbzmg(nnuz,ebs);
[Envxl,~] = xtoxbzmg(nnux,ebs);
[~,Enmuzu] = xtoxbzmg(nmuz,ebs);
[~,Enmuxu] = xtoxbzmg(nmux,ebs);
[Em0xl,~] = xtoxbzmg(m0x,ebs);

En0zzl=(pmu*exp(-mu)+pnu*exp(-v))*En0zl./p0;
En1zzl=(pmu*mu*exp(-mu)+pnu*v*exp(-v))*mu./(mu*v-v^2)*...
    (exp(v)*Envzl./pnu-v^2./mu^2*exp(mu)*Enmuzu./pmu...
    -(mu^2-v^2)./mu^2*En0zu./p0);
En1xxl=pw*w*exp(-w)*mu./(mu*v-v^2)*...
    (exp(v)*Envxl./pnu-v^2./mu^2*exp(mu)*Enmuxu./pmu...
    -(mu^2-v^2)./mu^2*En0xu./p0);

[n0zzl,~]=xbtoxzmg(En0zzl,ebs);
[n1zzl,~]=xbtoxzmg(En1zzl,ebs);
[n1xxl,~]=xbtoxzmg(En1xxl,ebs);


Em0wxl=exp(-w)*pw*Em0xl./p0;
[m0wxl,~]=xbtoxzmg(Em0wxl,ebs);
m1xxu=mwx-m0wxl;



e1zzu=m1xxu/n1xxl+gazmg(ebs,m1xxu/n1xxl,n1zzl,n1xxl);

eobs=mz./nz;
lamdec=nz*fec*H2(eobs);

ll=(n0zzl+n1zzl*(1-H2(e1zzu))-lamdec-6*log2(23/ebssec)-log2(2/ebscor))/N;

objv=-ll;
objv(e1zzu<0|eobs>1/2|m1xxu/n1xxl>1/2|e1zzu>1/2|n0zzl<0|n1xxl<0)=0;



end

