function [V,D,G,A] = lapeigs(im, n, eval);

% Number grid inside image>0
G = double(im)./double(max(im(:)));
p = find(G);
G(p) = (1:length(p))';

% Discrete laplacian on grid, with u=0 on boundary
A = delsq_noflux(G);

% Solve eigs of laplacian 
if isempty(n);
    n = length(A)-2;
end;
opts = struct();
opts.isreal = true;
opts.issym = true;
[V,D] = eigs(A,n,eval,opts);
return;

figure; 
for i=1:n;
    U = G;
    U(G>0) = full(V(G(G>0),i));
    imagesc(U);
    pause;
end;
