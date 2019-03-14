da=100;db=100;D=da+db;l0=da/2;
dx=50;
na=3.2;#slab A
nb=1;#slab B  
nc=2;
c=3e8;
w=0.6;
ka=w*1e7*na/c;
kb=w*1e7*nb/c;
kc=w*1e7*nc/c;
m=[[exp(im*ka*da) exp(-im*ka*da)];[ka*exp(im*ka*da) -ka*exp(-im*ka*da)]];
n=[[1 1];[kb -kb]];#[exp(im*kb*da),exp(-im*kb*da);kb*exp(im*kb*da),-kb*exp(-im*kb*da)];
o=[[exp(im*kb*(db-dx)) exp(-im*kb*(db-dx))];[kb*exp(im*kb*(db-dx)) -kb*exp(-im*kb*(db-dx))]];
p=[[1 1];[kc -kc]];
q=[[1+im*tan(kc*l0)/2 im*tan(kc*l0)/2];[-im*tan(kc*l0)/2 1-im*tan(kc*l0)/2]];
r=[[1 1];[kc -kc]];
s=[[1 1];[kb -kb]];
t=[[exp(im*kb*dx) exp(-im*kb*dx)];[kb*exp(im*kb*dx) -kb*exp(-im*kb*dx)]];
u=[[1 1];[ka -ka]];
T=inv(u)*t*inv(s)*r*q*inv(p)*o*inv(n)*m;
AAA=T[1,1]+T[2,2];
Temp0=AAA^2-4;
Temp1=(AAA-sqrt(Temp0))/2;
r0=abs((Temp1-T[1,1])/T[1,2]);
if r0<=1
    Temp1=(AAA-sqrt(Temp0))/2;
    r0=abs((Temp1-T[1,1])/T[1,2]);
    Arg=angle((Temp1-T[1,1])/T[1,2]);
else
    Temp1=(AAA+sqrt(Temp0))/2;
    r0=abs((Temp1-T[1,1])/T[1,2]);
    Arg=angle((Temp1-T[1,1])/T[1,2]);
end