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
PL = 1; % Left Pressure BC [Pa]
PR = 0; % Right Pressure BC[Pa]
pw = [1, 0]; % well pressure
PI = [10, 10]; % productivity index
cellno = [1,N]; % location of wells
wells = 1; % 1 - switch on wells; 0 - apply Dirichlet conditions

% Lambda = zeros(N,1); % Initialize Lambda, filled later

% generate mobilities:

Lambda  = 100.*(10.^rand(N,1)); % heterogenuous reservoir;
Lambda_hom  = ones(N,1); % homogenuous reservoir; 
% Lambda  = 100*ones(N,1); % homogenuous reservoir; 

LambdaH = Harmonic(Lambda,DX,N);
T = Trans(Lambda,DX,N);

% Plot Lambdas
figure(1)
subplot(1,3,1)
plot(x,Lambda, xi, LambdaH)
ylabel('Lambda')
xlabel('X')
title('Mobility')

% Plot permeability and harmonicly averaged fields

% figure1 = figure(1);
% axes1 = axes('Parent',figure1,'FontSize',14);
% plot(x,Lambda,'Parent',axes1,'Marker','.','LineStyle','-','color',[0 0 0],'DisplayName','L');
% hold on;
% plot(xi,LambdaH,'Parent',axes1,'Marker','*','LineStyle','-', 'color',[1 0 0],'DisplayName','{\lambda}^H');
% legend show;
% xlabel('x'); 
% ylabel('\lambda  and  {\lambda}^H'); 

q = zeros(N,1); % source function
BC = zeros(N,1); % Pressure boundary conditions

%% solver

A = zeros(N,N);
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

BC(1) = -T(1)*PL;
BC(N) = -T(N+1)*PR;

% hange A and q according to wells or 
switch(wells)
    case(0)
q = q - BC;
    case(1)        
[A,q] = AddWells( A, q, Lambda, pw, PI, cellno );
end        

%solver
p = A\q;


subplot(1,3,2)
plot(x,p)
ylabel('Pressure')
xlabel('X')
title('Pressure distribution')

%% calculate flow velocity

U = zeros (N+1, 1);
U(1) = - 2*LambdaH(1)*( p(1) - PL)/DX;
U(N+1) = - 2*LambdaH(N+1)*( PR - p(N))/DX;
for i = 2:N
U(i) = - LambdaH(i) * (p(i) - p(i-1))/DX;
end

subplot(1,3,3)
plot(xi,U)
ylabel('Velocity')
ylim([2 3])
xlabel('X')
title('Flow velocity')










