function [A,B,C,D]=cola_linearize(func,x0,u0)
%
% cola_linearize - This function linearizes xprime=f(x,u); y=g(x,u)
% It is written as a utility program for the distillation models,
% but it is a general purpose routine.
%    Written by: Atle C. Christiansen 
%    
%     SYNTAX: [A,B,C,D]=cola_linearize('func',Xinit',Uinit');
%
%     INPUTS:
%     func  - nonlinear model; must be on the form: [f,g]=func(x,u)
%     Xinit - row vector containing initial values for the states.
%     Uinit - row vector containing inputs.
%
%     OUTPUTS:
%     The outputs are the matrixes from the linear model:
%
%     xprime = A*x + B*u
%     y      = C*x + D*u
%
%     x - states
%     u - inputs and disturbances             
%     y - outputs
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n=max(size(x0));
r=max(size(u0));
delta=1e-6;

for i=1:n
  xh=x0; xh(i)=xh(i)+delta;
  xl=x0; xl(i)=xl(i)-delta;
  [fh,gh]=feval(func,xh,u0);
  [fl,gl]=feval(func,xl,u0);
  A(:,i)=(fh'-fl')/(2*delta);
  C(:,i)=(gh'-gl')/(2*delta);
end


for i=1:r
  uh=u0; uh(i)=uh(i)+delta;
  ul=u0; ul(i)=ul(i)-delta;
  [fh,gh]=feval(func,x0,uh);
  [fl,gl]=feval(func,x0,ul);
  B(:,i)=(fh'-fl')/(2*delta);
  D(:,i)=(gh'-gl')/(2*delta);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A1=A;B1=B;C1=C;D1=D;
delta=delta/2;


n=max(size(x0));
r=max(size(u0));


for i=1:n
  xh=x0; xh(i)=xh(i)+delta;
  xl=x0; xl(i)=xl(i)-delta;
  [fh,gh]=feval(func,xh,u0);
  [fl,gl]=feval(func,xl,u0);
  A(:,i)=(fh'-fl')/(2*delta);
  C(:,i)=(gh'-gl')/(2*delta);
end


for i=1:r
  uh=u0; uh(i)=uh(i)+delta;
  ul=u0; ul(i)=ul(i)-delta;
  [fh,gh]=feval(func,x0,uh);
  [fl,gl]=feval(func,x0,ul);
  B(:,i)=(fh'-fl')/(2*delta);
  D(:,i)=(gh'-gl')/(2*delta);
end

A=(4*A-A1)/3;
B=(4*B-B1)/3;
C=(4*C-C1)/3;
D=(4*D-D1)/3;      


