function [LambdaH] = Harmonic(Lambda,DX,N)
LambdaH = zeros(N+1, 1);

% calculate harmonic:
LambdaH(1)   = Lambda(1);
LambdaH(N+1) = Lambda(N);
LambdaH(2:N) = 2*Lambda(2:N).*Lambda(1:N-1)./(Lambda(2:N)+Lambda(1:N-1));

end