function [ h ] = hconv_int(H1,H2,W,dx,T)

H=H1+H2;
Tf=(T+773.15)/2;
B=1/Tf;
g=9.81;
ro=357.45*Tf^(-1.004);
Cp=(10^(-13)*Tf^4 - 6*10^(-10)*Tf^3 + 10^(-6)*Tf^2 - 0.0004*Tf + 1.0613)*1000;
u=3*10^(-7)*Tf^0.7197;
k=0.0205*exp(0.001*Tf);
Gr=(H^3*ro^2*g*B*(773.15-T))/u^2;
Pr=(Cp*u)/k;
Ra=Gr*Pr;
if Ra<10^9
    Nu=0.59*Ra^0.25;
else
    Nu=0.1*Ra^0.333;
end
h=(k*Nu)/H;

dc=(W*dx)/(2*W+2*dx);
Gr=(dc^3*ro^2*g*B*(773.15-T))/u^2;
Pr=(Cp*u)/k;
Ra=Gr*Pr;
Nu=0.27*Ra^0.25;
h=2*h+(k*Nu)/dc;
h=h/3;
end

