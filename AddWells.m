function [ A, q ] = AddWells( A, q, lambda, pw, PI, cellno )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
for i = 1: length(pw) % number of wells
A(cellno(i), cellno(i)) = A(cellno(i), cellno(i)) + PI(i)*lambda(cellno(i));
q(cellno(i)) = q(cellno(i)) + PI(i)*lambda(cellno(i))*pw(i);
end

end

