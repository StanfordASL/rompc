function [f,g]=cola4_lin(X,u)
%
% cola4_lin - This function designed for use with 'cola_linearize.m' to
%        create a linear model of a distillation column.
%        The number '4' denotes that there are 4 inputs (L,V,D,B)
%        plus the 2 disturbances. 
%        For more details about the column model see: colamod.m
%
%     x - states (liquid composition and liquid hold up)
%     u - inputs (reflux LT, boilup VB, distillate D, bottoms B) and
%          disturbances (feedrate F and feed composition zF)
%
%  NOTE: Everything here is row vectors rather than columnn vectors


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NT=41;

% Splitting the states
x=X(1:NT);                          % Liquid composition
M=X(NT+1:2*NT);                     % Liquid hold up

% Inputs
U(1)=u(1); % L
U(2)=u(2); % V
U(3)=u(3); % D
U(4)=u(4); % B
U(5)=u(5); % F
U(6)=u(6); % zF

% Another possible disturbance is the fraction of liquid in the feed,
% but it is not used here:
qF=1.0; U(7)=qF;

% This variable is not used
t=0;

xprime=colamod(t,X',U');

% Output
f=xprime';
g=[x(NT),x(1),M(NT),M(1)];



