function [T] = Trans(Lambda,DX,N)
T       = zeros(N+1, 1);

% calculate harmonic:
LambdaH = Harmonic(Lambda,DX,N);
% calculate transmissibility:
T(1) = 2*LambdaH(1)./(DX^2);
T(2:N) = LambdaH(2:N)./(DX^2);
T(N+1) = 2*LambdaH(N+1)./(DX^2);
end