function xprime=colamod(t,X,U) 
%
% colamod - This is a nonlinear model of a distillation column with
%           NT-1 theoretical stages including a reboiler (stage 1) plus a
%           total condenser ("stage" NT). The liquid flow dynamics are
%           modelled by a simple linear relationship.
%           Model assumptions: Two components (binary separation); constant
%           relative volatility; no vapor holdup; one feed and two products;
%           constant molar flows (same vapor flow on all stages); 
%           total condenser
%
%           The model is based on column A in Skogestad and Postlethwaite
%           (1996). The model has 82 states.
%
% Inputs:    t    - time in [min].
%            X    - State, the first 41 states are compositions of light
%                   component A with reboiler/bottom stage as X(1) and 
%                   condenser as X(41). State X(42)is holdup in reboiler/
%                   bottom stage and X(82) is hold-up in condenser. 
%            U(1) - reflux L,
%            U(2) - boilup V,
%            U(3) - top or distillate product flow D,
%            U(4) - bottom product flow B,
%            U(5) - feed rate F,
%            U(6) - feed composition, zF.
%            U(7) - feed liquid fraction, qF.
%
% Outputs:   xprime - vector with time derivative of all the states 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------------------------------
% The following data need to be changed for a new column.
% These data are for "colmn A".
% Number of stages (including reboiler and total condenser: 
    NT=41; 
% Location of feed stage (stages are counted from the bottom):
    NF=21;
% Relative volatility
    alpha=1.5;
% Nominal liquid holdups
    M0(1)=0.5;      	% Nominal reboiler holdup (kmol)
    i=2:NT-1; M0(i)=0.5*ones(1,NT-2);% Nominal stage (tray) holdups (kmol)
    M0(NT)=0.5;      	% Nominal condenser holdup (kmol)
% Data for linearized liquid flow dynamics (does not apply to reboiler and condenser):
    taul=0.063;     	% time constant for liquid dynamics (min)
    F0=1;	 	% Nominal feed rate (kmol/min) 
    qF0 = 1; 		% Nominal fraction of liquid in feed 
    L0=2.70629;     	% Nominal reflux flow (from steady-state data)
    L0b=L0 + qF0*F0;	% Nominal liquid flow below feed (kmol/min)
    lambda=0;		% Effect of vapor flow on liquid flow ("K2-effect")
    V0=3.20629;V0t=V0+(1-qF0)*F0;% Nominal vapor flows - only needed if lambda is nonzero 
% End data which need to be changed
%------------------------------------------------------------

% Splitting the states
x=X(1:NT)';                          % Liquid composition from btm to top
M=X(NT+1:2*NT)';                     % Liquid hold up from btm to top

% Inputs and disturbances
LT = U(1);                            % Reflux
VB = U(2);                            % Boilup
D  = U(3);                            % Distillate
B  = U(4);                            % Bottoms

F  = U(5);                            % Feedrate
zF = U(6);                            % Feed composition
qF = U(7);                            % Feed liquid fraction

% THE MODEL

% Vapor-liquid equilibria
i=1:NT-1;    y(i)=alpha*x(i)./(1+(alpha-1)*x(i));

% Vapor Flows assuming constant molar flows
i=1:NT-1;    V(i)=VB*ones(1,NT-1);
i=NF:NT-1;   V(i)=V(i) + (1-qF)*F;

% Liquid flows assuming linearized tray hydraulics with time constant taul
% Also includes coefficient lambda for effect of vapor flow ("K2-effect").
i=2:NF;      L(i) = L0b + (M(i)-M0(i))./taul + lambda.*(V(i-1)-V0);
i=NF+1:NT-1; L(i) = L0  + (M(i)-M0(i))./taul + lambda.*(V(i-1)-V0t);
L(NT)=LT;

% Time derivatives from  material balances for 
% 1) total holdup and 2) component holdup

% Column
i=2:NT-1;
dMdt(i) = L(i+1)         - L(i)       + V(i-1)         - V(i);
dMxdt(i)= L(i+1).*x(i+1) - L(i).*x(i) + V(i-1).*y(i-1) - V(i).*y(i);

% Correction for feed at the feed stage
% The feed is assumed to be mixed into the feed stage
dMdt(NF) = dMdt(NF)  + F;
dMxdt(NF)= dMxdt(NF) + F*zF;

% Reboiler (assumed to be an equilibrium stage)
dMdt(1) = L(2)      - V(1)      - B;
dMxdt(1)= L(2)*x(2) - V(1)*y(1) - B*x(1);

% Total condenser (no equilibrium stage)
dMdt(NT) = V(NT-1)         - LT       - D;
dMxdt(NT)= V(NT-1)*y(NT-1) - LT*x(NT) - D*x(NT);

% Compute the derivative for the mole fractions from d(Mx) = x dM + M dx
i=1:NT;
dxdt(i) = (dMxdt(i) - x(i).*dMdt(i) )./M(i);

% Output
xprime=[dxdt';dMdt'];

