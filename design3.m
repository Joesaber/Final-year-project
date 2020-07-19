function [ Tp1,Tp2,Tc,Pb,Tw,Tout,X,w,Bakingtime,Velocity,Length,Power,A,C] = design3


prompt = 'What is the desired production rate? ';
Q = input(prompt);
A=zeros(19,1);
C=zeros(19,1);

% Specific heat: %
Cp_c= 605.88 ;
Cp_a= 1093 ;

% Densities: %
ro_c= 7870 ;
ro_v=451.6;
ro_a=1.225;

% Beta: %
B=0.4242;

% Overall conductive heat transfer coefficient for walls: %
U_w_out=0.182;

% Contact resistance: %
R=3*10^(-3);

% Atmospheric pressure: %
Patm= 101325;

% Water vapor;s pressure in the air: %
Pwv_a=10000;

% Molecular weights: %
Mw=18.01528;
Ma=28.9647;

% Known temperatures: %
Ta=773 ;
Tamb=298 ;

% Lewis number:%
Le=1;

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
s=0.03;
f=0.15;

efficiency=0.3;


for j=2:20
    
    tb=j;
    L=(((Q*tb)/3600)*((2*Rp + f)*(2*Rp + s))/(Wc-s));
    dx=0.2;
    g=(L*10)/(dx*10); 
 n=ceil(g)+1;
Vc=L/tb;

q1=zeros(n,1);
q2=zeros(n,1);
q3=zeros(n,1);
q4=zeros(n,1);
q5=zeros(n,1);
q6=zeros(n,1);
q7=zeros(n,1);
q1(1)=460;
q2(1)=307.15;
q3(1)=570.15;
q5(1)=577;
q6(1)=311.5;
q7(1)=0.56;
z=zeros(n,1);
p=zeros(n,1);
p(1)=0.36;

for i =2:n 
   if i==n
       if rem(L,dx)~=0
           dx=rem(L,dx)*dx;
       end
   end
    z(i) = z(i-1) + dx;
    
     % Product's properties: %
    Cp_p= 837 + 3349*p(i-1) ;
    ro_p= 513 *(1 + p(i-1)) ;
    k_p= -0.5677 + 1.396*p(i-1) + 0.005131*(q2(i-1) - 273.15);
    
    
    % Convective heat transfer coefficients: %
    hc_p2_a=hconv_product(Rp,q2(i-1));
    hc_c_a=hconv_conveyer(dx,Wc,Rp,B,q3(i-1));
    hc_w_a=hconv_int(H1,H2,W,dx,q5(i-1));
    hc_out_amb=hconv_ext(Hext,Wext,dx,q6(i-1));
    
    % Radiative heat trasnfer coefficients: %
    [ hr_w1_c , hr_w2_c , hr_w1_p2, hr_c_w1 , hr_c_w2 , hr_p2_w1] = Radiation(B,dx,W,Wc,H1,H2);
    
    % Water activity: %
    aw=(((100*q7(i-1))/(exp(-0.0056*q2(i-1) +5.5)))^(-1/0.38)+1)^(-1);
    
    Psat=133.3*exp(18.3036-((3816.44)/(q2(i-1)-46.13)));
    
    % Water vapor pressure at the surface of the product: %
    Pwv_s= aw*Psat;
    
    
    % Convective mass transfer coefficient: %
    k_g= (( (Ma/(Mw*hc_p2_a)) * Patm * Cp_a * Le^(2/3) )^(-1))*3*10^(-2); 
    
    % Flux of water evaporation: %
    N= k_g*(Pwv_s - Pwv_a);
    
    Hvap=2442500*((647.3-q2(i-1))/349.1)^0.38; 
    
    a0=[q1(i-1);q2(i-1);q3(i-1);q4(i-1);q5(i-1);q6(i-1);q7(i-1)];
    
    lsqOpts = optimoptions('fsolve','MaxFunEvals', 1e+30,'MaxIter',1e+30,'Display','iter','TolFun',1e-10,'tolX',1e-10);
    asol=fsolve(@(a) des3( a,q7(i-1),q2(i-1),q5(i-1),q6(i-1),q1(i-1),q3(i-1),L,tb,dx),a0,lsqOpts);
    
    
    q1(i)=asol(1);
    q2(i)=asol(2);
    q3(i)=asol(3);
    q4(i)=asol(4);
    q5(i)=asol(5);
    q6(i)=asol(6);
    q7(i)=asol(7);
    p(i)=(q7(i)/(1+q7(i)));
    

    a5=hc_p2_a*W*dx*(Ta-q2(i)) + hc_c_a*(2-B)*W*dx*(Ta-q3(i)) + hc_w_a*(2*H1 + 2*H2 +W )*dx*(Ta-q5(i)) + hr_w1_c*(2*H1 + W)*dx*abs(q3(i)^4 - q5(i)^4) + hr_c_w1*(1-B)*W*dx*abs(q5(i)^4 - q3(i)^4) + hr_w2_c*(2*H2)*dx*abs(q3(i)^4 - q5(i)^4) + hr_c_w2*W*dx*abs(q5(i)^4 - q3(i)^4) + hr_p2_w1*B*W*dx*abs(q5(i)^4 - q2(i)^4) + hr_w1_p2*(2*H1 + W)*dx*abs(q2(i)^4 - q5(i)^4) + U_w_out*(((2*H1 + 2*H2)/W) +1)*W*dx*abs(q5(i)-q6(i)) + ((k_p)/(0.9*tp))*B*W*dx*abs(q1(i)-q2(i)) + (1/R)*B*W*dx*abs(q3(i)-q1(i)) + abs(N*W*dx*B*Hvap);
    
    q4(i)=(q4(i)+a5)/efficiency;

    
    if i==n
    A(j-1)=q7(i);
    C(j-1)=q2(i);
    end
    
    if i==n
        if (q2(i)<480) && (q2(i)>400) && (q7(i)<0.42) && (q7(i)>0.2)
            Bakingtime=tb;
            Velocity=Vc;
            Length=L;
            Tp1=q1;
            Tp2=q2;
            Tc=q3;
            Pb=q4;
            Tw=q5;
            Tout=q6;
            X=q7;
            x=z;
            w=p;
            Power = sum(Pb);
        end
    end
            
end
end
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