function circleEigs(w,h,cx,cy);

im = zeros(w,h);
[x,y] = meshgrid([1:w],[1:h]);
R = sqrt((x-cx).*(x-cx) + (y-cy).*(y-cy));

figure; 
hold on;
for r=10:14;
    % Circle of increasing radius r
    bim = R<r;
    %imagesc(bim);
    
    % Compute eigenvectors of laplacian on this domain
    [V,D,G,A] = lapeigs(bim, 20, 0);
    %figure;
    for i=1:5;
        U = G;
        U(G>0) = full(V(G(G>0),end-i+1));
        subplot(5,5,i+(r-10)*5);
        imagesc(U); axis image; colorbar; title(D(i,i));
    end;
    %plot(diag(D));
    pause(0.1);
end;