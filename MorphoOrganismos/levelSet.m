function L = levelSet(L0, speed0, dt, n);
% Evolve level set in L using speed by n time steps dt

s = size(L0);
w = s(1); h = s(2);
L = zeros(w,h,n);
L(:,:,1) = L0;
speed = speed0;
for i=2:n+1;
    [gy,gx] = gradient(L(:,:,i-1));
    %[gx,gy] = upwind_grad(igvx, igvy, L(:,:,i-1));
    g = sqrt(gx.*gx+gy.*gy);
    mask = L(:,:,i-1)<=0;
    speed = mask.*speed0;
    [bw,idx] = bwdist(speed);
    speed(find(mask==0)) = speed(idx(find(mask==0))); % Put background speed as that of nearest boundary point
    L(:,:,i) = L(:,:,i-1) - g.*speed*dt;
    
    % plotting
    %subplot(1,2,1); imagesc(L(:,:,i)<0); colorbar;
    %subplot(1,2,2); imagesc(g.*(L(:,:,i)>0).*(L(:,:,i)<1)); colorbar;
    %pause(0.1);
end;

L = L(:,:,end);