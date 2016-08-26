 function [V,D,G] = lapeigs(im, n);

% Number grid inside image>0
G = double(im)./double(max(im(:)));
p = find(G);
G(p) = (1:length(p))';

% Discrete laplacian on grid, with u=0 on boundary
A = delsq(G);

% Solve eigs of laplacian 
[V,D] = eigs(A,n,'sm');
return;

figure; 
for i=1:n;
    U = G;
    U(G>0) = full(V(G(G>0),i));
    imagesc(U);
    pause;
end;