function [ R ] = ComputeResidual( rho_n, phi_n, rho, phi, Lambda, dt, p, q, N, DX ) % V for finite volume
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% rho, phi at ineration nu
R = zeros(N,1);

T = Trans(Lambda,DX,N);

R(1) = rho(1)*q(1) - (rho(1)*phi(1) - rho_n(1)*phi_n(1))/dt - ((T(1) + T(2))*p(1) - T(2)*p(2));
for i = 2:N-1
R(i) = rho(i)*q(i) - (rho(i)*phi(i) - rho_n(i)*phi_n(i))/dt - (T(i)*(p(i)-p(i-1)) + T(i+1)*(p(i)-p(i+1)));
end
R(N) = rho(N)*q(N) - (rho(N)*phi(N) - rho_n(N)*phi_n(N))/dt - ((T(N)+T(N+1))*p(N) - T(N)*p(N-1));

end

