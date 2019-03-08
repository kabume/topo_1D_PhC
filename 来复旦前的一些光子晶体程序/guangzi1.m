na=2.10;nb=1.46;n1=1;n2=1;
c3=0;c1=asin(n1*sin(c3)/na);c2=asin(na*sin(c1)/nb);c4=asin(nb*sin(c2)/n2);
d1=1064;d2=400;d3=664;
f=4*pi*1e-7;
e=1e-9/(36*pi);
m=sqrt(e/f);
a=d2/(2*na);b=d3/(2*nb);
c=a+b;
for d=100:1600;
    Ba=2*pi*na*a*cos(c1)/d;
    Bb=2*pi*nb*b*cos(c2)/d;
    za=m*cos(c1)*na;zb=m*cos(c2)*nb;z1=f*cos(c3)*n1;z2=f*cos(c4)*n2;
    p1=cos(Bb);p2=-i*sin(Bb)/zb;p3=-i*zb*sin(Bb);p4=cos(Bb);
    P=[p1 p2;p3 p4];
    q1=cos(Ba);q2=-i*sin(Ba)/za;q3=-i*za*sin(Ba);q4=cos(Ba);
    Q=[q1 q2;q3 q4];
    O=Q*P;k=acos((O(1,1)+O(2,2))/2)/(a+b);
    s=d-99;
    h(1,s)=k;
end
d=100:1600;
x=h.*(c/pi);
y=d1./d;
plot(x,y,'k')
