%model D,Dr.Li 's paper + 2 column
% This program could simulate the paper of Li and Zww
clc,clear
%paramaters
da=100e-9;db=100e-9;D=da+db;lc=da / 3;
Dt=D*1e9;
dx0=0:0.1:Dt;
dx=dx0*1e-9;
na=3.2;
nb=1;
nc=2;
c=3e8;
w=0.1:0.01:15;
ka=w*1e15*na/c;
kb=w*1e15*nb/c;
kc=w*1e15*nc/c;
[r0, Arg] = ABA(da, db, lc, na, nb, nc, w, dx);
%PCLs(da, db, na, nb, w);


figure
colormap(jet);
imagesc(dx0,w,r0');
set(gca,'YDir','normal');  
title("|r| of semiinfinite")
xlabel("\Delta(nm)");ylabel("\omega(10^{15} Hz)");
hold on
PCLs(da, db, na, nb, w);
hold off
figure
colormap(jet);
imagesc(dx0,w,Arg');
set(gca,'YDir','normal');
title("Arg of semiinfinite")
xlabel("\Delta(nm)");ylabel("\omega(10^{15} Hz)");
hold on 
PCLs(da, db, na, nb, w);
hold off