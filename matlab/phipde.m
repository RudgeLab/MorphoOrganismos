function f = phipde(inphi, r, w, h, alpha);
%---
% Compute grad(r).grad(phi) + phi*del2(r)
% RHS of differential eqn for phi = speed
% given a radial direction vector (grad(r))
%
% Eqn is:
%   mu(x,y) = grad(r).grad(phi) + phi*del2(r)
%
%---

% Pad phi with zeros for boundary condition
%phi = zeros(w,h);
%phi(:,2:end) = reshape(inphi, w, h-1);
phi = reshape(inphi, w, h);

r = reshape(r, w, h);

[drdx, drdy] = gradient(r);
[dpdx, dpdy] = gradient(phi); %upwind_grad(drdx, drdy, phi);

% Boundary conditions phi(:,1) = 0
dpdy(1,:) = phi(2,:)*0.5;
%phi(1,:) = 0;

f = (drdx.*dpdx + drdy.*dpdy) + phi.*del2(r)*4; % Note: the *4 is stupid matlab because del2 returns 1/4 laplacian!
f = f + alpha*phi; % regularisation
f = f(:);
