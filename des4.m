function [ fval ] = des4( A,X,Tp2,Tw,Tout,Tp1,Tc,L,tb,dx)
% Specific heat: %
Cp_c= 605.88 ;
Cp_a= 1093 ;

% Densities: %
ro_c= 7870 ;
ro_v=451.6;

% Contact resistance: %
R=3*10^(-3);

% Conveyer's velocity: %
Vc=L/tb;

% Atmospheric pressure: %
Patm= 101325;

% Molecular weights: %
Mw=18.01528;
Ma=28.9647;

% Beta: %
B=0.424;

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


% Overall conductive heat transfer coefficient for walls: %
U_w_out=0.182;

% Temperatures: %
Ta=773 ;
Tamb=298 ;




D=4*10^(-6);
Pwv_a=10000;
Le=1;

    
    
    % Water content (wet basis): %
    w=(X/(1+X));
    
    % Product's properties: %
    Cp_p= 837 + 3349*w ;
    ro_p= 513 *(1 + w) ;
    k_p= -0.5677 + 1.396*w + 0.005131*(Tp2 - 273.15);
    
    
    % Convective heat transfer coefficients: %
    hc_p2_a=hconv_product(Rp,Tp2);
    hc_c_a=hconv_conveyer(dx,Wc,Rp,B,Tc);
    hc_w_a=hconv_int(H1,H2,W,dx,Tw);
    hc_out_amb=hconv_ext(Hext,Wext,dx,Tout);
    
    % Radiative heat trasnfer coefficients: %
    [ hr_w1_c , hr_w2_c , hr_w1_p2, hr_c_w1 , hr_c_w2 , hr_p2_w1] = Radiation(B,dx,W,Wc,H1,H2);
    
    % Water activity: %
    aw=(((100*X)/(exp(-0.0056*Tp2 +5.5)))^(-1/0.38)+1)^(-1);
    
    Psat=133.3*exp(18.3036-((3816.44)/(Tp2-46.13)));
    
    % Water vapor pressure at the surface of the product: %
    Pwv_s= aw*Psat;
    
    
    % Convective mass transfer coefficient: %
    k_g= (( (Ma/(Mw*hc_p2_a)) * Patm * Cp_a * Le^(2/3) )^(-1))*3*10^(-2); 
    
    % Flux of water evaporation: %
    N= k_g*(Pwv_s - Pwv_a);
    
  
    
    
    Hvap=2442500*((647.3-Tp2)/349.1)^0.38; 
    
    e1=A(1);
    e2=A(2);
    e3=A(3);
    e4=A(4);
    e5=A(5);
    e6=A(6);
    e7=A(7);
    
    % e1=Tp1(i) ; e2=Tp2(i) ; e3=Tc(i) ; e4=Pb(i) ; e5=Tw(i) ; e6=Tout(i) ; e7=X(i);
    % Tp1=Tp1(i-1) ; Tp2=Tp2(i-1) ; Tc=Tc(i-1) ;
   
    
     fval(1,1) = ((ro_p*Cp_p*0.1*tp*Vc)/(dx))*Tp1 - (((ro_p*Cp_p*0.1*tp*Vc)/(dx)))*e1 - ((1/R) + ((k_p)/(0.9*tp)))*e1 + (1/R)*e3 + ((k_p)/(0.9*tp))*e2 ;
    
     fval(2,1) = ((ro_p*Cp_p*0.9*tp*Vc)/(dx))*Tp2 - (((ro_p*Cp_p*0.9*tp*Vc)/(dx)))*e2 -( ((k_p)/(0.9*tp)) +  hc_p2_a )*e2 - hr_p2_w1*(e2)^4 +((k_p)/(0.9*tp))*e1 + hc_p2_a*Ta + hr_p2_w1*(e5)^4 - N*Hvap;
    
     fval(3,1) = ((ro_c*Cp_c*tc*Vc)/(dx))*Tc - ((ro_c*Cp_c*tc*Vc)/(dx))*e3 -( + hc_c_a*(2-B) + (B/R) )*e3 - (hr_c_w1*(1-B) + hr_c_w2 )*(e3)^4 + hc_c_a*(2-B)*Ta + (B/R)*e1 + (hr_c_w1*(1-B)+hr_c_w2)*(e5)^4;
    
     fval(4,1) = hc_w_a*(((2*H1 + 2*H2)/W) +1)*e5 + hc_p2_a*B*e2 + hc_c_a*(2-B)*e3 + (1/(W*dx))*e4 + B*N*Hvap - ( hc_w_a*(((2*H1 + 2*H2)/W) +1) + hc_p2_a*B + hc_c_a*(2-B))*Ta;
    
     fval(5,1) = hc_w_a*(((2*H1 + 2*H2)/W) +1)*Ta + hr_w1_p2*( ((2*H1)/W) +1)*(e2)^4 + ( hr_w1_c*( ((2*H1)/W) +1) + hr_w2_c*((2*H2)/W))*(e3)^4 + U_w_out*(((2*H1 + 2*H2)/W) +1)*e6 - ( hr_w1_p2*( ((2*H1)/W) +1) + hr_w1_c*( ((2*H1)/W) +1) + hr_w2_c*((2*H2)/W) )*(e5)^4 - ( hc_w_a*(((2*H1 + 2*H2)/W) +1) + U_w_out*(((2*H1 + 2*H2)/W) +1) )*e5;
    
     fval(6,1) = hc_out_amb*Tamb + U_w_out*e5 - (hc_out_amb + U_w_out)*e6;
    
     fval(7,1) = -D*(ro_v/dx)*X + D*(ro_v/dx)*e7 + N;
     
     
     
end
