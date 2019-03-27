using Plots
using Unicode
using LaTeXStrings
gr()
da=100e-9;db=200e-9;D=da+db;l0=da;#单位m
Dt=D*1e9;
dx=0:0.1:Dt;
dx=dx*1e-9
na=2;#slab A
nb=1;#slab B
nc=2;
c=3e8;
w=0.1:0.01:15;
ka=w*1e15*na/c;
kb=w*1e15*nb/c;
kc=w*1e15*nc/c;
num1=length(w);
num2=length(dx);
r0=zeros(num2,num1);
Arg=zeros(num2,num1);
Em=zeros(num2,num1);
Temp2=zeros(num2,num1)+im*zeros(num2,num1);
Threads.@threads for nn=1:num1
#for nn=1:num1
    for mm=1:num2
        if mm<db/D*num2
            m=[[exp(im*ka[nn]*da) exp(-im*ka[nn]*da)];[ka[nn]*exp(im*ka[nn]*da) -ka[nn]*exp(-im*ka[nn]*da)]];
            n=[[1 1];[kb[nn] -kb[nn]]];#[exp(im*kb[nn]*da),exp(-im*kb[nn]*da);kb[nn]*exp(im*kb[nn]*da),-kb[nn]*exp(-im*kb[nn]*da)];
            o=[[exp(im*kb[nn]*(db-dx[mm])) exp(-im*kb[nn]*(db-dx[mm]))];[kb[nn]*exp(im*kb[nn]*(db-dx[mm])) -kb[nn]*exp(-im*kb[nn]*(db-dx[mm]))]];
            p=[[1 1];[kc[nn] -kc[nn]]];
            q=[[1+im*tan(kc[nn]*l0)/2 im*tan(kc[nn]*l0)/2];[-im*tan(kc[nn]*l0)/2 1-im*tan(kc[nn]*l0)/2]];
            r=[[1 1];[kc[nn] -kc[nn]]];
            s=[[1 1];[kb[nn] -kb[nn]]];
            t=[[exp(im*kb[nn]*dx[mm]) exp(-im*kb[nn]*dx[mm])];[kb[nn]*exp(im*kb[nn]*dx[mm]) -kb[nn]*exp(-im*kb[nn]*dx[mm])]];
            u=[[1 1];[ka[nn] -ka[nn]]];
            T=inv(u)*t*inv(s)*r*q*inv(p)*o*inv(n)*m;
        else
            m=[[exp(im*ka[nn]*(da+db-dx[mm])) exp(-im*ka[nn]*(da+db-dx[mm]))];[ka[nn]*exp(im*ka[nn]*(da+db-dx[mm])) -ka[nn]*exp(-im*ka[nn]*(da+db-dx[mm]))]];
            n=[[1 1];[kc[nn] -kc[nn]]];
            o=[[1+im*tan(kc[nn]*l0)/2 im*tan(kc[nn]*l0)/2];[-im*tan(kc[nn]*l0)/2 1-im*tan(kc[nn]*l0)/2]];
            p=[[1 1];[kc[nn] -kc[nn]]];
            q=[[1 1];[ka[nn] -ka[nn]]];
            r=[[exp(im*ka[nn]*(dx[mm]-db)) exp(-im*ka[nn]*(dx[mm]-db))];[ka[nn]*exp(im*ka[nn]*(dx[mm]-db)) -ka[nn]*exp(-im*ka[nn]*(dx[mm]-db))]];
            s=[[1 1];[kb[nn] -kb[nn]]];
            t=[[exp(im*kb[nn]*db) exp(-im*kb[nn]*db)];[kb[nn]*exp(im*kb[nn]*db) -kb[nn]*exp(-im*kb[nn]*db)]];
            u=[[1 1];[ka[nn] -ka[nn]]];
            T=inv(u)*t*inv(s)*r*inv(q)*p*o*inv(n)*m;
        end
        AAA=T[1,1]+T[2,2];
        Temp0=AAA^2-4;
        Temp1=(AAA-sqrt(Temp0))/2;
        r0[mm,nn]=abs((Temp1-T[1,1])/T[1,2]);
        if r0[mm,nn]<=1#导带取一种情况
            Temp1=(AAA-sqrt(Temp0))/2;
            r0[mm,nn]=abs((Temp1-T[1,1])/T[1,2]);
            Arg[mm,nn]=angle((Temp1-T[1,1])/T[1,2]);
        else#此时并不都是禁带，r0仍可能为1，所以取出r0等于1的值与上面导带的情况简并
            Temp1=(AAA+sqrt(Temp0))/2;
            r0[mm,nn]=abs((Temp1-T[1,1])/T[1,2]);
            if r0[mm,nn]<=1.0001&&r0[mm,nn]>=0.9999#1为int所以加一些小数位
                Temp1=(AAA-sqrt(Temp0))/2;
                r0[mm,nn]=abs((Temp1-T[1,1])/T[1,2]);
                Arg[mm,nn]=angle((Temp1-T[1,1])/T[1,2]);
            else
                Temp1=(AAA+sqrt(Temp0))/2;
                r0[mm,nn]=abs((Temp1-T[1,1])/T[1,2]);
                Arg[mm,nn]=angle((Temp1-T[1,1])/T[1,2]);
            end
        end
    end
end
#heatmap(dx,w,r0')
#savefig("lc=da反射系数.pdf")
heatmap(dx,w,Arg')
xlabel!("Delta(mm)")
ylabel!("Normalized frequency")
title!("Arg of semi-infinite system")
#heatmap(dx,w,r0')
#xlabel!("Delta(mm)")
#ylabel!("Normalized frequency")
#title!("|r| of semi-infinite system")


#savefig("lc=da相位.pdf")
#heatmap(dx,w,Em')
