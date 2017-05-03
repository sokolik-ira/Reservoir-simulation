% 1D Elliptic FD & FV Solver
% Incompressible flow

clear all
close all
%% Input data and grid configuration:
% num of interfaces = num of cells + 1

L = 1; % Length [m]
N = 20; % Number of grid cells
DX = L / N; % Dx [m]
x = linspace(DX/2, L - DX/2, N); %Grid centers
xi = linspace(0, L, N + 1); %inteface locations

pw = [1, 0]; % well pressure
PI = [1000, 1000]; % productivity index
cellno = [1,N]; % location of wells

Lambda  = ones(N,1); % homogenuous reservoir
% Lambda  = 100.*(10.^rand(N,1)); % heterogenuous reservoir;

dt = 0.01;
Nt = 10;%round(t_fin/dt);

A = zeros(N,N);

%% define IC
% compute porosity and density

rho_o = 1; % initial density
c1 = 1;
phi_o = 0.3; % initial density
cr = 1;
p_o = 0; % initial pressure
% p = 1:0.5:10.5;
p = zeros(N,1);

[rho, drho_dp] = ComputeDensity( p, rho_o, p_o, c1);
[phi, dphi_dp] = ComputePorosity( p, phi_o, p_o, cr);
T = Trans(rho.*Lambda,DX,N);
q = ComputeWellFluxes( pw, p, PI, Lambda, cellno, N );

%% time loop
for Nt = 1:100
    %Neuton loop
    Converged = 0;
    rho_n = rho;
    phi_n = phi;
    % Find P_n+1
    Residual = ComputeResidual(rho_n, phi_n, rho, phi, T, dt, p, q, N)
    
    while Converged == 0
        % Build Jacobian J = A + C
        % Assemble A:
        % cell 1:
        A = zeros (N,N);
        A(1,1) = T(1) + T(2);
        A(1,2) = -T(2);
        % cell N:
        A(N,N) = T(N+1) + T(N);
        A(N,N-1) = -T(N);
        % cell i:
        for i = 2:N-1
            A(i,i) = T(i+1)+T(i);
            A(i,i-1) = -T(i);
            A(i,i+1) = -T(i+1);
        end
        
        % Assemble C:
        vec = (rho .* dphi_dp + phi .* drho_dp)/dt;
        C = diag(vec);
        
        % Add wells
        
        W = zeros(N,N);
        for i = 1:length(PI)
            W(cellno(i),cellno(i)) = W(cellno(i),cellno(i)) + Lambda(cellno(i))*PI(i)*rho(cellno(i));
        end
        
        J = A + C+ W;
        
        dp = J\Residual;
        p = p + dp;
        
        
        % Update fluid and rock propertied
        [rho, drho_dp] = ComputeDensity( p, rho_o, p_o, c1);
        [phi, dphi_dp] = ComputePorosity( p, phi_o, p_o, cr);
        T = Trans(rho.*Lambda,DX,N);
        q = ComputeWellFluxes( pw, p, PI, Lambda, cellno, N );
        
        % Compure Residual
        Residual = ComputeResidual(rho_n, phi_n, rho, phi, T, dt, p, q, N) ;
        
        % Convergence criteria
        if norm(Residual)<1e-6
            Converged = 1;
        end
    end
    
    plot(p); hold on;
end


