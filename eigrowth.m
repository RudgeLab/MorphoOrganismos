function eigrowth(domain, pattern, dt, n);

[w,h] = size(domain);

% Initialise distance map = level set func phi
phi0 = bwdist(domain>0) - bwdist(1-(domain>0));
phi = zeros(w,h,n);
phi(:,:,1) = phi0;

% Initialise pattern u
u = zeros(w,h,n);
pattern = pattern/(max(pattern(:)));
u(:,:,1) = pattern;

% peak in spectrum
d0 = 0.5;
% slope of spectrum
kd = 100; 

for i=2:n+1;
    % Update phi by growing with speed u.*u
    phi(:,:,i) = levelSet(phi(:,:,i-1), u(:,:,i-1) - min(min(u(:,:,i-1))), dt, 1);
    %phi(:,:,i) = levelSet(phi(:,:,i-1), ones(w,h), dt, 1);
    
    % Compute eigenvectors and eigenvalues (100 closest to d0
    [V,D,G] = lapeigs(phi(:,:,i)<=0, 50, d0);
    d = diag(D);    % ordered list of eigenvalues
    
    % Project previous pattern onto eigenvectors
    uprev = u(:,:,i-1);
    ww = V(G(G>0),:)'*uprev(G>0);
    
    % Compute dispersion relation values, quadratic
    lam = 1 - kd*(d-d0).*(d-d0);
    
    % New pattern is sum of eigenvectors weighted by lam*w
    uu = G;
    uu(G>0) = V(G(G>0),:)*(exp(lam*dt).*ww);
    uu = uu/max(uu(:));
    u(:,:,i) = reshape(real(uu),w,h);
    
    % Plotting
    subplot(1,2,1); imagesc(phi(:,:,i)<0); axis image; colorbar;
    subplot(1,2,2); imagesc(u(:,:,i)); axis image; colorbar;
    pause(0.1);
end;