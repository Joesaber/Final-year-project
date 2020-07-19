function [ Tp1,Tp2,Tc,Pb,Tw,Tout,X,w,Power ] = design4


% Dimensions: %
Rp=0.125;
tc=0.005;
tp=0.002;
H1=0.7;
H2=0.7;
Hext=2;
Wext=1.5;
W=0.8;
Wc=0.6;
s=0.15;
L=1.57;
tb=17;
dx=0.2;
B=0.424;
U_w_out=0.182;
R=3*10^(-3);

efficiency=0.3;


g=(L*10)/(dx*10);
n=ceil(g)+1;

Ta=773 ;
Tamb=298 ;

% Initial values:%
Tp1=zeros(n,1);
Tp1(1)=460;

Tp2=zeros(n,1);
Tp2(1)=307.15;

Tc=zeros(n,1);
Tc(1)=570.15;

Tw=zeros(n,1);
Tw(1)=577;

Tout=zeros(n,1);
Tout(1)=311.5;

X=zeros(n,1);
X(1)=0.56;

w=zeros(n,1);
w(1)=0.36;

x=zeros(n,1);

Pb=zeros(n,1);


Cp_a= 1093 ;

% Atmospheric pressure: %
Patm= 101325;

% Molecular weights: %
Mw=18.01528;
Ma=28.9647;
Pwv_a=10000;
Le=1;


for i =2:n 
     if i==n
       if rem(L,dx)~=0
           dx=rem(L,dx)*dx;
       end
     end
    
      % Convective heat transfer coefficients: %
    hc_p2_a=hconv_product(Rp,Tp2(i-1));
    hc_c_a=hconv_conveyer(dx,Wc,Rp,B,Tc(i-1));
    hc_w_a=hconv_int(H1,H2,W,dx,Tw(i-1));
    hc_out_amb=hconv_ext(Hext,Wext,dx,Tout(i-1));
    
    % Radiative heat trasnfer coefficients: %
    [ hr_w1_c , hr_w2_c , hr_w1_p2, hr_c_w1 , hr_c_w2 , hr_p2_w1] = Radiation(B,dx,W,Wc,H1,H2);
    
    x(i) = x(i-1) + dx;
    
    k_p= -0.5677 + 1.396*w(i-1) + 0.005131*(Tp2(i-1) - 273.15);
    
    % Water activity: %
    aw=(((100*X(i-1))/(exp(-0.0056*Tp2(i-1) +5.5)))^(-1/0.38)+1)^(-1);
    
    Psat=133.3*exp(18.3036-((3816.44)/(Tp2(i-1)-46.13)));
    
    % Water vapor pressure at the surface of the product: %
    Pwv_s= aw*Psat;
    
    
    % Convective mass transfer coefficient: %
    k_g= (( (Ma/(Mw*hc_p2_a)) * Patm * Cp_a * Le^(2/3) )^(-1))*3*10^(-2); 
    
    % Flux of water evaporation: %
    N= k_g*(Pwv_s - Pwv_a);
    
  
    
    
    Hvap=2442500*((647.3-Tp2(i-1))/349.1)^0.38; 
    
    a0=[Tp1(i-1);Tp2(i-1);Tc(i-1);Pb(i-1);Tw(i-1);Tout(i-1);X(i-1)];
    lsqOpts = optimoptions('fsolve','MaxFunEvals', 1e+30,'MaxIter',1e+30,'Display','iter','TolFun',1e-10,'tolX',1e-10);
    asol=fsolve(@(a) des4( a,X(i-1),Tp2(i-1),Tw(i-1),Tout(i-1),Tp1(i-1),Tc(i-1),L,tb,dx),a0,lsqOpts);
    
    
    Tp1(i)=asol(1);
    Tp2(i)=asol(2);
    Tc(i)=asol(3);
    Pb(i)=asol(4);
    Tw(i)=asol(5);
    Tout(i)=asol(6) ;
    X(i)=asol(7);
    Heat_transfers= hc_p2_a*W*dx*(Ta-Tp2(i)) + hc_c_a*(2-B)*W*dx*(Ta-Tc(i)) + hc_w_a*(2*H1 + 2*H2 +W )*dx*(Ta-Tw(i)) + hr_w1_c*(2*H1 + W)*dx*abs(Tc(i)^4 - Tw(i)^4) + hr_c_w1*(1-B)*W*dx*abs(Tw(i)^4 - Tc(i)^4) + hr_w2_c*(2*H2)*dx*abs(Tc(i)^4 - Tw(i)^4) + hr_c_w2*W*dx*abs(Tw(i)^4 - Tc(i)^4) + hr_p2_w1*B*W*dx*abs(Tw(i)^4 - Tp2(i)^4) + hr_w1_p2*(2*H1 + W)*dx*abs(Tp2(i)^4 - Tw(i)^4) + U_w_out*(((2*H1 + 2*H2)/W) +1)*W*dx*abs(Tw(i)-Tout(i)) + ((k_p)/(0.9*tp))*B*W*dx*abs(Tp1(i)-Tp2(i)) + (1/R)*B*W*dx*abs(Tc(i)-Tp1(i)) + abs(N*W*dx*B*Hvap);
    Pb(i)=Pb(i) + Heat_transfers;
    Pb(i)=Pb(i)/(efficiency);
    w(i)=(X(i)/(1+X(i)));
    
    
end
Power = sum(Pb);
figure
subplot(4,4,1);
stairs(x,Tp1);
xlabel('Distance');
ylabel('Lower product''s temperature(K)');

subplot(4,4,2);
stairs(x,Tp2);
xlabel('Distance');
ylabel('Top product''s temperature(K)');

subplot(4,4,3);
stairs(x,Tc);
xlabel('Distance');
ylabel('Conveyer''s temperature(K)');

subplot(4,4,4);
stairs(x,X);
xlabel('Distance');
ylabel('Water content');

subplot(4,4,5);
stairs(x,Tw);
xlabel('Distance');
ylabel('Wall''s temperature(K)');

subplot(4,4,6);
stairs(x,Tout);
xlabel('Distance');
ylabel('Outer wall''s temperature(K)');

subplot(4,4,[ 7 8 ]);
stairs(x,Pb);
xlabel('Distance');
ylabel('Power');


end

