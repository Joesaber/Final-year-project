function [ hr_w1_c , hr_w2_c , hr_w1_p2, hr_c_w1 , hr_c_w2 , hr_p2_w1 ] = Radiation(B,dx,W,Wc,H1,H2)

% Emissivities: %

eps_w = 0.26;
eps_c = 0.26;
eps_p = 0.85;


% Stefan-Boltzmann constant: % 
si= 5.669*10^(-8);

S_w1 = 2*H1*dx + W*dx;
S_c = Wc*dx*(1-B);
S_p = Wc*dx*B;
S_w2 = 2*H2*dx;

% Shape factors: %

F_w2_c =0.12 ;

F_c_w2 =0.14 ;


hr_w1_c = si * (1/((1/eps_w)+((1-eps_c)/eps_c)*(S_w1/S_c)));


hr_w2_c = si * ( 1/((1/F_w2_c)+((1-eps_w)/eps_w)+((1-eps_c)/eps_c)*(S_w2/(Wc*dx))));


hr_w1_p2 = si * (1/((1/eps_w)+((1-eps_p)/eps_p)*(S_w1/S_p)));


hr_c_w1 = si * (1/((1/eps_c)+((1-eps_w)/eps_w)*(S_c/S_w1)));


hr_c_w2 = si * ( 1/((1/F_c_w2)+((1-eps_c)/eps_c)+((1-eps_w/eps_w)*((Wc*dx)/S_w2))));


hr_p2_w1 = si * (1/((1/eps_p)+((1-eps_w)/eps_w)*(S_p/S_w1)));



end
