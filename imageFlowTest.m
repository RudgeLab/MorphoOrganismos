ims=testim;
%ims(:,:,1)=sum(ims,3);
%ims(:,:,2)=0;

H = fspecial('gaussian',3,3);
%ims = imfilter(double(ims), H, 'symmetric');

figure;
s = size(ims);
w = s(1); h = s(2);
xidx = 2:w-1;
yidx = 2:h-1;

gxpos = zeros(w,h,3);
gxneg = zeros(w,h,3);
gypos = zeros(w,h,3);
gyneg = zeros(w,h,3);

u = zeros(w,h);

% Compute an eigenvector of initial shape for reference
[V,D,G] = lapeigs(ims(:,:,1)>0.2, 100, 1.25);
d = diag(D);
%idx = find((d<1.5).*(d>1)); % indices of eigenvalues 1.5<lambda<2
l = 1-(d-1.25).*(d-1.25);
idx = find(l>0); % quadratic dispersion
%Vm = mean(V(:,idx),2); % linear sum of eigenvectors weighted equally (flat dispersion relation)
Vm = V(:,idx)*exp(10*l(idx));
Uprev = G;
Uprev(G>0) = full(Vm(G(G>0))); % Just mapping sparse thing onto grid
Vprev = V(:,idx);
dprev = d(idx);;

for i=1:1000;
    du = 10*del2(u) + ims(:,:,3) - u;
    u = u + 0.01*du;
    
    % Growth
    bims = ims(:,:,1)>0.2; %sum(ims>0.2,3)>0;    
    dbims = bwdist(bims) - bwdist(1-bims);
    %dbims = imfilter(double(dbims), H, 'symmetric'); % smooth kinks in distance map due to pixelisation of starting shape
    %dbims = imfilter(double(dbims), H, 'symmetric');
    v = 1./(1+dbims.*dbims.*dbims.*dbims/100);
    %v = 1 - (dbims/min(dbims(:))).^2;
    
    % growth regulation by colour
    gims = imfilter(ims, H, 'symmetric');
    %red = ims(:,:,1).*u./(1+u);
    %blue = 0.75*ims(:,:,3); %0.75*ims(:,:,3).*u./(1+u);
    red = ims(:,:,1);
    blue = ims(:,:,2);
    %gfac = blue; %(red*0.1 + blue*1)./(red+blue+0.1);
    
    % Compute eigenvectors of laplacian on domain and use to compute speed
    [V,D,G] = lapeigs(ims(:,:,1)>0.2, 100, 1.25);
    d = diag(D);
    l = 1-(d-1.25).*(d-1.25);
    idx = find(l>0); % quadratic dispersion
    %idx = find((d<1.5).*(d>1)); % indices of eigenvalues 1.5<lambda<2
    nv = length(idx);
    
    % Find multiplicity of eigenvalues (assume only 2x)
    %mul = find(dprev==d(idx+1));
    
    
    % construct matrix of eigenvectors in range
    Uall = zeros(w*h,nv);
    for vv=1:nv;
        Uall(G>0,vv) = full( V(G(G>0),idx(vv)) );
    end;
    % Project previous pattern onto new eigenspace
    Uproj = Uall*(exp(10*l(idx)).*(Uall'*Uprev(:)));
    Uproj = reshape(Uproj,w,h);
    %Uproj = Uproj./max(Uproj(:));
    %Uprev = Uproj; % reset previous pattern to new
    Vprev = V(:,idx);
    dprev = d(idx);
    [gux,guy] = gradient(Uproj/max(Uproj(:)));
    gfac = (gux.*gux+guy.*guy)*1e1; % totally adhoc speed function
    %v = v.*gfac;
    
    % Radial direction vector from distance map
    [gvy,gvx] = gradient(dbims); %(sin(x/30).*cos(y/50))');
    igvx = zeros(w,h,3);igvy = zeros(w,h,3);
    igvx(:,:,1) = gvx;
    igvx(:,:,2) = gvx;
    igvx(:,:,3) = gvx;
    igvy(:,:,1) = gvy;
    igvy(:,:,2) = gvy;
    igvy(:,:,3) = gvy;
    
    [gx,gy] = upwind_grad(igvx, igvy, ims);
    
    % Copy velocity to r,g,b channels
    vv = zeros(w,h,3);vv(:,:,1)=v; vv(:,:,2)=v; vv(:,:,3)=v;

    % Scale velocity to give uniform strain=divergence
    %velx = vv.*igvx; vely = vv.*igvy;
    %mu = divergence(vely(:,:,1), velx(:,:,1));
    %vv = vv./reshape(repmat(mu,1,3),w,h,3);

    % Magnitude of gradient (density gradient==pressure)
    %diffims = -vv.*(gx.*gx + gy.*gy);

    diffims = vv.*(igvx.*gx + igvy.*gy); % this part is advection    
    % advection of mass == img(:,:,3)
    % advection + diffusion of signals == img(:,:,1:2)
    
    ims = ims - 0.1*diffims;
    %ims(:,:,3) = ims(:,:,3) + 0.1*(ims(:,:,2) - mu*ims(:,:,3));
    
    %ims(:,:,1) = ims(:,:,1).*bims;
    %ims(:,:,2) = ims(:,:,2).*bims;
    %ims(:,:,3) = ims(:,:,3).*bims;

    subplot(1,3,1);
    imagesc(Uproj);colorbar; axis image; title('Uproj');
    subplot(1,3,2);
    plot(d(idx), l(idx), '.');
    %imagesc(v);colorbar; axis image; title('Uprev');
    Uprev = Uproj;
    %imshow(255*ims/max(ims(:)));
    %plot(ims(33,:,1)); %colorbar; axis image;
    subplot(1,3,3);
    %imagesc(diffims(:,:,1)); colorbar; axis image;
    imagesc((ims(:,:,1)>0.2)); colorbar; axis image;
    %imagesc(dbims); colorbar; axis image;
    fname = sprintf('frame%04d.png',i);
    %imwrite(ims, fname, 'png');
    
    disp(diag(D)); % print eigenvalues
    pause(0.1); 
end;