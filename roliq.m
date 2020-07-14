function rhol = roliq(pressure,T)
 
% A simple liquid dens model wich takes into pressure varations vs. pressure
% is implemented. P0 is the atmosperic pressure. D0 is density at surface
% conditions

  rou0 = 10.829;% 
  T0 = 68;
 
  beta=-0.0001535;
  alpha=2.7384e-06;
  gamma=-7.469e-07;
 
  rhol = rou0*exp(alpha*pressure + beta*(T-T0)+gamma*(T-T0)^2 );