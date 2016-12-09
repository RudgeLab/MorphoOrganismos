function circleEigs(w,h,cx,cy);

im = zeros(w,h);
[x,y] = meshgrid([1:w],[1:h]);
R = sqrt((x-cx).*(x-cx) + (y-cy).*(y-cy));

figure; 
hold on;
for r=1:4;
    % Circle of increasing radius r
    bim = R<r*5;
    %imagesc(bim);
    
    % Compute eigenvectors of laplacian on this domain
    [V,D,G] = lapeigs(bim, 20, 0);
    %figure;
    for i=1:6;
        U = G;
        k = length(D)-i+1;
        U(G>0) = full(V(G(G>0),k));
        subplot(5,6,i+r*6);
        imagesc(U, [-0.01,0.01]); axis image; colorbar; title(sprintf('k^2 = %03.2e',D(k,k)));
    end;
    %plot(diag(D));
    pause(0.1);
end;