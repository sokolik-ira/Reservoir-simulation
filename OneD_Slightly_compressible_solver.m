% 1D Elliptic FD & FV Solver
% Incompressible flow

clear all
close all
%% Input data and grid configuration:
% num of interfaces = num of cells + 1

L = 100; % Length [m]
N = 20; % Number of grid cells
DX = L / N; % Dx [m]
x = linspace(DX/2, L - DX/2, N); %Grid centers
xi = linspace(0, L, N + 1); %inteface locations

% boundaty conditions

wells = 1; % 1 - switch on wells; 0 - apply Dirichlet conditions
PL = 1; % Left Pressure BC [Pa]
PR = 0; % Right Pressure BC[Pa]
pw = [1, 0]; % well pressure
PI = [10, 10]; % productivity index
cellno = [1,N]; % location of wells

% Lambda = zeros(N,1); % Initialize Lambda, filled later

% generate mobilities:

Lambda  = 100.*(10.^rand(N,1)); % heterogenuous reservoir;
Lambda_hom  = ones(N,1); % homogenuous reservoir; 
% Lambda  = 100*ones(N,1); % homogenuous reservoir; 

LambdaH = Harmonic(Lambda,DX,N);
T = Trans(Lambda,DX,N);

% Plot Lambdas
% figure(1)
% subplot(1,3,1)
% plot(x,Lambda, xi, LambdaH)
% ylabel('Lambda')
% xlabel('X')
% title('Mobility')

% Plot permeability and harmonicly averaged fields

q = zeros(N,1); % source function
BC = zeros(N,1); % Pressure boundary conditions

%% Slightly compressible

Ceff = 0.1;
phi = 0.2;

t_fin = 0.5;
dt = 0.000001;
Nt = round(t_fin/dt);

C = eye(N,N)*phi*Ceff/dt;

explicit = 0;
P_t = zeros(N,10);
% P_t(:,1) = 0; 

A = zeros(N,N);
q = zeros(N,1); % source function

% cell 1:
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

% hange A and q according to wells or 
switch(wells)
    case(0)
q = q - BC;
    case(1)        
[A,q] = AddWells( A, q, Lambda, pw, PI, cellno );
end 
      
switch(explicit) 
     case(0)
       for j = 2:1000000 
            P_t(:,j) = (C+A)\(q+C*P_t(:,j-1));
       end
     case(1)   
       for j = 2:1000000
            P_t(:,j) = C\(q-A*P_t(:,j-1)+C*P_t(:,j-1));
       end
end
%%
for i = 1:10000:100000
    plot(x,P_t(:,i))
    hold on
 end








