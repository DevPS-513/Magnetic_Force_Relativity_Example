clear all
close all
clc
% This code solves for the force on a moving charg, or "cat" using magnetic equations
% and only electric field equation in conjunction with special relativity, it verifies, but does not validate
% the explaination given in the youtube video https://www.youtube.com/watch?v=1TKSfAkWWN0
% Author: Max B.W


% units
meters=1;
cm=(1/100)*meters;
c=299792457.8;          % speed of light

mu=4*pi*1e-7;
eps_o=(1/(mu*c*c));
format long e
%% Moving Cat frame electric force

I=10E3;             % current in amps
r=.05*cm;           % radius of wire
A_wire=pi*r^2;      % area of the wire
z=100;              % distance of cat from wire
q_cat=1;            % charge of the cat

rho_c=(8.94)/(cm^3);    % grams/cm^3
rho_c=10000;
Aw_c=63.546;            % g/mol of copper;
Av=6.0221413e+23;       % Atoms/mol;

mol_per_m_cubed=rho_c/Aw_c;
Atomic_density=mol_per_m_cubed*Av;
num_e=Atomic_density;
columb_density=Atomic_density*1.6E-19;


% drift velocity for a given current
v_d=I/(A_wire*columb_density);
v_cat=.001*c;

lam_o=columb_density*A_wire;

gamma_vd=sqrt(1-(v_d/c)^2);

gamma_cat_p=sqrt(1-((v_d+(v_cat-v_d))/c)^2);


gamma_cat_e=sqrt(1-((v_d)/c)^2);


gamma_diff=sqrt(1-((v_cat-v_d)/c)^2);

gamma_sum=sqrt(1-((v_cat+v_d)/c)^2);

  %1.000002417580733e+00
ratio_v_to_c=v_d/c;

%% ELECTRON FRAME
% linear charge density in the electron frame of reference
lam_neg_e_frame=lam_o*gamma_vd;
% linear charge density of positive ions seen by the electrons
lam_pos_e_frame=lam_o/gamma_vd;

% net charge in the electron frame
lam_tot_e_frame=lam_pos_e_frame-lam_neg_e_frame;
%% CAT FRAME
% linear charge density in the electron frame of reference
lam_neg_cat_frame=lam_o*(gamma_cat_e)/gamma_diff;
% linear charge density of positive ions seen by the electrons
lam_pos_cat_frame=lam_o/(gamma_cat_p);

% net charge in the electron frame
lam_tot_cat_frame=lam_pos_cat_frame-lam_neg_cat_frame;



E_theory=lam_o*gamma_vd*v_d^2/(2*pi*eps_o*z)

E_sim=lam_tot_e_frame/(2*pi*eps_o*z)



F_cat_relativity=lam_tot_cat_frame/(2*pi*eps_o*z)

% now force in rest fram due to magetic field


B=mu*I/(2*pi*z);
F_cat_rest=v_cat*B*q_cat

ratio=F_cat_relativity/F_cat_rest






