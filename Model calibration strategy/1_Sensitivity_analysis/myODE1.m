%%=======================================================================%%
%  Author: Rafael da Silva Raqueti
%  Advisor: Samuel da Silva
%  On the Calibration of Reduced-Order Models to Describe the 
%  Viscoelasticity in Steady-State Rolling Tires 
%  Methodology >> 1_Sensitivity_analysis
%    myODE1.m
%%=======================================================================%%

function dY = myODE1(t,Y,C,X)
%
%  myODE1
%    This function defines the Reduced-Order Model that is a system of
%    First-Order Nonlinear Ordinary Differential Equations to simulate
%    viscoelastic internal variables.
%
%  USAGE dY = myODE1(t,Y,C,X)
%__________________________________________________________________________
%  OUTPUTS
%    dY : State vector derivative [\dot{A}, \dot{Z}] where:
%      - \dot{A} is the viscoelastic internal variable derivative;
%      - \dot{Z} is the histeretic output derivative.
%__________________________________________________________________________
%  INPUTS
%    t : Time vector;
%    Y : State vector [A, Z] where:
%      - A is the viscoelastic internal variable;
%      - Z is the histeretic output.
%    C : Component of the Right Cauchy-Green tensor:
%      - C = griddedInterpolant(t,data) to perform interpolation.
%    X : Vector of unknown parameters:
%      - X(1) = c;
%      - X(2) = k;
%      - X(3) = \alpha;
%      - X(4) = \gamma;
%      - X(5) = \delta;
%      - X(6) = \nu.
%__________________________________________________________________________
dY = zeros(2,1);
dY(1,1) = (C(t)-X(2)*Y(1,1)-Y(2,1))/X(1);
dY(2,1) = X(3)*dY(1,1)-...
          X(4)*abs(dY(1,1))*abs(Y(2,1))^(X(6)-1)*Y(2,1)-...
          X(5)*dY(1,1)*abs(Y(2,1))^X(6);
end