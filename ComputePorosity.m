function [ phi, dphi_dp ] = ComputePorosity( p, phi_o, p_o, cr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% rho       = zeros(N, 1);
phi  = phi_o .* exp(cr .* (p-p_o));
dphi_dp = phi_o .* cr .* exp(cr .*(p-p_o));

end