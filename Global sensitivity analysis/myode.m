%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rafael da Silva Raqueti
% Bruna Silveira Pavlack
% João Pedro Canisso Valese Norenberg
% Luccas Pereira Miguel
% Uncertainty Quantification
% Global sensitivity analysis
%  - myode.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script defines the system of first-order nonlinear differential
% equations.

function dY = myode(t,Y,C,X)

% t is a time variable;
% Y = {A,Z} is the computational model response in which:
%   - A is the internal variable output;
%   - Z is the histeretic output.
% C(t) is the input component of the right Cauchy-Green tensor:
%   - Use C = griddedInterpolant(t,dataC) to perform time interpolation.
% X = {c,k,\alpha,\gamma,\delta,\nu} is the set of parameters.

dY = zeros(2,1);
dY(1,1) = (C(t)-X(2)*Y(1,1)-Y(2,1))/X(1);
dY(2,1) = X(3)*dY(1,1) - ...
          X(4)*abs(dY(1,1))*abs(Y(2,1))^(X(6)-1)*Y(2,1) - ...
          X(5)*dY(1,1)*abs(Y(2,1))^X(6);
end