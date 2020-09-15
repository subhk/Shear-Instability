clear all;
clc;
% =========================================================
% This code solves non-dimensional form of non-Boussinesq 
% Taylor-Goldstein equation using finite difference matrix. 
% Formuations are taken from Shete&Guha - JFM.
% ========================================================= 
  global D1 D2 D4  
% ------------------------------------
% Input parameters :
% ------------------------------------
N  = 401;  % grid points
H  = 5;   % half length of the vertical domain
z  = linspace(-H,H,N)'; % origin at the middle of the domain
% -------------------------------------
% Differentiation matrix
% -------------------------------------
D1 = ddz_4(z);
D2 = ddz2_4(z);
D4 = ddz4_4(z);
% -------------------------------------
% Boundary conditions
% -------------------------------------
vel_bc = 0;  % 1-> no slip;    0-> stress free 
rho_bc = 0;  % 1-> insulating; 0-> fixed
% -------------------------------------
% Basic state profile
% -------------------------------------
% imposed shear profile
U = tanh(z);
DU = D1*U;  
D2U = D2*U;
% imposed stable stratified density profile
R = 3;          % Ratio of shear layer thickness to density thickess     
rho_ = -tanh(R*z);
Drho = D1*rho_;

invRe = 1/300;  %inverse of Reynolds number (1/Re)
invPr = 1/R^2;  %inverse of Prandtl number (1/Pr) 
J = 0.;
imode = 1;
alpha = 0.9;

[gamma_max]=TG_nonBoussinesq(U, D2U, J, z, Drho, alpha, invRe, invPr, imode, N, vel_bc, rho_bc);

gamma_max







