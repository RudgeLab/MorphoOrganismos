function [V,D] = fulllapeigs(w,h,n);
% Compute eigenvectors of full matrix laplacian

L = buildLaplacian(w,h);
[V,D] = eigs(L,n,'lm');

figure;
for i=1:n;
    imagesc(reshape(V(:,i),w.h));
    colorbar; axis image;
    pause;
end;
