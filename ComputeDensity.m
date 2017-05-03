function [ rho, drho_dp ] = ComputeDensity( p, rho_o, p_o, c1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% rho       = zeros(N, 1);
rho  = rho_o .* exp(c1.*(p-p_o));
drho_dp = rho_o .* c1 .* exp(c1*(p-p_o));

end

