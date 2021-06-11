
function [rho_g] = lenz_project(antenna_length, theta0, epsilon_r)
  nr = sqrt(epsilon_r);
  z = antenna_length*cos(theta0);
  rho = antenna_length*sin(theta0);
  a = (nr-1)/(nr+1);
  b = 2*rho/(nr+1);
  c = -(rho^2 + z*z*(nr*nr)/(nr*nr-1));
  rho_g = max(roots([a b c]))
end
  