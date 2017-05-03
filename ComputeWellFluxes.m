function [ q ] = ComputeWellFluxes( pw, p, PI, lambda, cellno, N )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
q = zeros(N,1);
for i = 1: length(pw) % number of wells
q(cellno(i)) =  PI(i)*lambda(cellno(i))*(pw(i)-p(cellno(i)));
end
end

