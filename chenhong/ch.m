clc,clear
%paramaters
L=100;D=2*L;h0=L/2;
Delta=-0.8:0.001:0.8;
h=(1+Delta)*L;h1=L-h/2;h2=h;ne=2;c=3e8;
w=0:0.001:0.6;
k=w*1e7*ne/c;
%tic
%semiinfinite(ne,k,h0,w,Delta,h1,h2);
%toc
bloch(w,k,D,h0,Delta)
%onecell(k,h0,w,Delta,h1,h2)
