function [gamma_max] = TG_nonBoussinesq(U,D2U,J,z,rho_bz,alpha,invRe,invPr,imode,N, vel_bc, rho_bc)

    global D2 D2b D4
    
    del=mean(diff(z));

    ii = sqrt(-1);

% ----------------------------------------    
% implementing boundary condition
% ---------------------------------------

% Impenetrable BC on 2nd order derivative matrix for velocity
  D2(1,:) = 0;
  D2(1,1) = -2/del^2;
  D2(1,2) = 1/del^2;  

  D2(2,:) = 0; 
  D2(2,1) = (4/3)/del^2; 
  D2(2,2) = -(5/2)/del^2;
  D2(2,3) = (4/3)/del^2;
  D2(2,4) = -(1/12)/del^2;     

  D2(N,:) = 0;
  D2(N,N) = -2/del^2;
  D2(N,N-1) = 1/del^2;  

  D2(N-1,:)   = 0;
  D2(N-1,N)   = (4/3)/del^2;
  D2(N-1,N-1) = -(5/2)/del^2;
  D2(N-1,N-2) = (4/3)/del^2;
  D2(N-1,N-3) = -(1/12)/del^2;
  D2(N-1,N-4) = 0;  

% On 4th order derivative matrix for velocity -> no slip or, stress free condition
% no slip (vel_bc=1) && free slip (vel_bc=0)

  if vel_bc == 0   
    D4(1,:) = 0;                 D4(1,1) = 5/del^4;
    D4(1,2) = -4/del^4;          D4(1,3) = 1/del^4;

    D4(2,:) = 0;                 D4(2,1) = -(38/6)/del^4;
    D4(2,2) = (28/3)/del^4;      D4(2,3) = -(13/2)/del^4;
    D4(2,4) = 2/del^4;           D4(2,5) = -(1/6)/del^4;

    D4(3,:) = 0;                 D4(3,1) = 2/del^4;
    D4(3,2) = -(13/2)/del^4;     D4(3,3) = (28/3)/del^4;
    D4(3,4) = -(13/2)/del^4;     D4(3,5) = 2/del^4;
    D4(3,6) = -(1/6)/del^4;

    D4(N,:) = 0;                 D4(N,N)   = 5/del^4;
    D4(N,N-1) = -4/del^4;        D4(N,N-2) = 1/del^4;

    D4(N-1,:) = 0;               D4(N-1,N)   = -(38/6)/del^4;
    D4(N-1,N-1) = (28/3)/del^4;  D4(N-1,N-2) = -(13/2)/del^4;
    D4(N-1,N-3) = 2/del^4;       D4(N-1,N-4) = -(1/6)/del^4;

    D4(N-2,:) = 0;                 D4(N-2,1) = 2/del^4;
    D4(N-2,2) = -(13/2)/del^4;     D4(N-2,3) = (28/3)/del^4;
    D4(N-2,4) = -(13/2)/del^4;     D4(N-2,5) = 2/del^4;
    D4(N-2,6) = -(1/6)/del^4;
  end
 
  if vel_bc == 1   
    D4(1,:) = 0;                 D4(1,1) = 7/del^4;
    D4(1,2) = -4/del^4;          D4(1,3) = 1/del^4;

    D4(2,:) = 0;                 D4(2,1) = -(40/6)/del^4;
    D4(2,2) = (28/3)/del^4;      D4(2,3) = -(13/2)/del^4;
    D4(2,4) = 2/del^4;           D4(2,5) = -(1/6)/del^4;;

    D4(N,:) = 0;                 D4(N,N) = 7/del^4;
    D4(N,N-1) = -4/del^4;        D4(N,N-2) = 1/del^4;

    D4(N-1,:) = 0;               D4(N-1,N)   = -(40/6)/del^4;
    D4(N-1,N-1) = (28/3)/del^4;  D4(N-1,N-2) = -(13/2)/del^4;
    D4(N-1,N-3) = 2/del^4;       D4(N-1,N-4) = -(1/6)/del^4;
  end

% -----------------------------------------------  
% Boundary condition on density
% -----------------------------------------------

% Fixed-buoyancy boundary
  if rho_bc == 0 
    D2b(1,:) = 0;
    D2b(1,1) = -2/del^2;
    D2b(1,2) = 1/del^2;  

    D2b(2,:) = 0; 
    D2b(2,1) = (4/3)/del^2; 
    D2b(2,2) = -(5/2)/del^2;
    D2b(2,3) = (4/3)/del^2;
    D2b(2,4) = -(1/12)/del^2;     

    D2b(N,:) = 0;
    D2b(N,N) = -2/del^2;
    D2b(N,N-1) = 1/del^2;  

    D2b(N-1,:)   = 0;
    D2b(N-1,N)   = (4/3)/del^2;
    D2b(N-1,N-1) = -(5/2)/del^2;
    D2b(N-1,N-2) = (4/3)/del^2;
    D2b(N-1,N-3) = -(1/12)/del^2;
    D2b(N-1,N-4) = 0;
  end 

% Insulating boundaries for density
  if rho_bc == 1
    D2b(1,:) = 0;
    D2b(1,1) = -2/(3*del^2);
    D2b(1,2) = 2/(3*del^2);
 
    D2b(2,:) = 0;
    D2b(2,1) = 11/(9*del^2);
    D2b(2,2) = -89/(36*del^2);
    D2b(2,3) = 4/(3*del^2);
    D2b(2,4) = -1/(12*del^2);

    D2b(N,:) = 0;
    D2b(N,N) = -2/(3*del^2);
    D2b(N,N-1) = 2/(3*del^2);

    D2b(N-1,:) = 0;
    D2b(N-1,1) = 11/(9*del^2);
    D2b(N-1,2) = -89/(36*del^2);
    D2b(N-1,3) = 4/(3*del^2);
    D2b(N-1,4) = -1/(12*del^2);
  end


    I  = eye(N);
    
    L  = D2  - alpha^2*I;
    Lb = D2b - alpha^2*I;
    LL = D4 - 2*alpha^2*D2 + alpha^4*I;
  
    A11 = diag(U)*L - diag(D2U) + (ii*invRe)/(alpha)*LL;
    A12 = -J*I;
    
    A21 = -diag(rho_bz);
    A22 = diag(U) + (ii*invRe*invPr)/(alpha)*Lb;
    
    A = [A11 A12; A21 A22];
    
    B11 = L;
    B12 = 0*I;
    
    B21 = 0*I;
    B22 = I;
    
    B = [B11 B12; B21 B22];
    
    % Solving the eigenvalue problem
    [eig_n,lambda] = eig(A, B);  % Eq. (25) in manual

    % the complex phase-speeds are the eigenvalues
    lambda=diag(lambda);
    cr=real(lambda); % real phase-speed
    ci=imag(lambda); % imaginary phase-speed (positive means growing disturbance)

    % Sort eigvalues, from most unstable to least unstable
    [~,ind]=sort(-ci);
    ci=ci(ind);
    cr=cr(ind); 
    eig_n=eig_n(:,ind); % eigenvectors sorted accordingly  

   % growth-rate of most unstable mode(s)
    gamma_max=alpha*ci(1:imode); % outputs only the least stable mode
end
