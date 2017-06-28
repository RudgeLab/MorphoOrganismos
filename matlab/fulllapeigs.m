function [V,D] = fulllapeigs(w,h,n,v0);
% Compute eigenvectors of full matrix laplacian

L = buildLaplacian(w,h);
opts = struct();
opt.v0 = v0; %zeros(length(L));
[V,D] = eig(L); %eigs(L,n,0,opts);
return;

figure;
for i=1:n;
    imagesc(reshape(V(:,i),w,h));
    colorbar; axis image;
    pause;
end;
