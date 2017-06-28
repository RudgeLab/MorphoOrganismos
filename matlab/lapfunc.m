function L=lapfunc(im, bim);
% Compute laplacian of im on domain bim>0, with zero flux BCs 
im = reshape(im, size(bim));
im = im.*(bim>0); % make sure zero outside domain
L = del2(im); % This has zero BC
%return;
L=-L(:);
return;

% Gradient of domain is +ve at left/top boundary, -ve at right/bottom
[gy,gx] = gradient(1*(bim>0)); % Yes Matlab returns them in the wrong order
gx = gx./max(gx(:));
gy = gy./max(gy(:));

% At the left boundary we want L(x,y) = L(x,y) + im(x+1,y)*0.25
% At the top boundary we want L(x,y) = L(x,y) + im(x,y+1)*0.25

% At the right boundary we want L(x,y) = L(x,y) + im(x-1,y)*0.25
% At the bottom boundary we want L(x,y) = L(x,y) + im(x,y-1)*0.25

% Hence:
% In x-axis
idx = find((gx~=0).*(bim>0));
L(idx) = L(idx) + im(idx+gx(idx))*0.25;

% In y-axis
idx = find((gy~=0).*(bim>0));
[ix,iy] = ind2sub(size(gy),idx);
idx2 = sub2ind(size(gy),ix,iy+gy(idx));
L(idx) = L(idx) + im(idx2)*0.25;
L=L(:);
