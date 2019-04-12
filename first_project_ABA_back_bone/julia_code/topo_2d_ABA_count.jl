using Plots
using Unicode
using LaTeXStrings
gr()
GC.gc()
da=100;db=100;D=da+db;l0=da;
na=3.2;#slab A
nb=1;#slab B
nc=2;
c=3e8;
j=1;
w=0.00001:0.001:0.6;
dx=0.00001:0.1:D;
num1=length(w);
num2=length(dx);
r0=zeros(num2,num1);
Arg=zeros(num2,num1);
Em=zeros(num2,num1);
Temp2=zeros(num2,num1)+im*zeros(num2,num1);
#Threads.@threads for nn=1:num1
for w=0.00001:0.01:0.6;
    ka=w*1e7*na/c;
    kb=w*1e7*nb/c;
    kc=w*1e7*nc/c;
    i=1;
    for dx=0.00001:1:D;
        if dx<db
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
        else
            m=[[exp(im*ka*(da+db-dx)) exp(-im*ka*(da+db-dx))];[ka*exp(im*ka*(da+db-dx)) -ka*exp(-im*ka*(da+db-dx))]];
            n=[[1 1];[kc -kc]];
            o=[[1+im*tan(kc*l0)/2 im*tan(kc*l0)/2];[-im*tan(kc*l0)/2 1-im*tan(kc*l0)/2]];
            p=[[1 1];[kc -kc]];
            q=[[1 1];[ka -ka]];
            r=[[exp(im*ka*(dx-db)) exp(-im*ka*(dx-db))];[ka*exp(im*ka*(dx-db)) -ka*exp(-im*ka*(dx-db))]];
            s=[[1 1];[kb -kb]];
            t=[[exp(im*kb*db) exp(-im*kb*db)];[kb*exp(im*kb*db) -kb*exp(-im*kb*db)]];
            u=[[1 1];[ka -ka]];
            T=inv(u)*t*inv(s)*r*inv(q)*p*o*inv(n)*m;
        end
        AAA=T[1,1]+T[2,2];
        Temp0=AAA^2-4;
        Temp1=(AAA-sqrt(Temp0))/2;
        r0[i,j]=abs((Temp1-T[1,1])/T[1,2]);
        if r0[i,j]<1
            Temp1=(AAA-sqrt(Temp0))/2;
            r0[i,j]=abs((Temp1-T[1,1])/T[1,2]);
            Em[i,j]=abs(T[1,2]+(Temp1-T[1,1]));
            Arg[i,j]=angle((Temp1-T[1,1])/T[1,2]);
        else
            Temp1=(AAA+sqrt(Temp0))/2;
            r0[i,j]=abs((Temp1-T[1,1])/T[1,2]);
            Arg[i,j]=angle((Temp1-T[1,1])/T[1,2]);
            Temp2[i,j]=(Temp1-T[1,1])/T[1,2];
            Em[i,j]=abs(T[1,2]+(Temp1-T[1,1]));
        end
        i=i+1;
    end
    global j=j+1;
end
#heatmap(dx,w,r0')
#savefig("lc=da反射系数.pdf")
heatmap(dx,w,Arg')
xlabel!("Delta(mm)")
ylabel!("Normalized frequency")
title!("Arg of semi-infinite system")
heatmap(dx,w,r0')
xlabel!("Delta(mm)")
ylabel!("Normalized frequency")
title!("|r| of semi-infinite system")


#savefig("lc=da相位.pdf")
#heatmap(dx,w,Em')
