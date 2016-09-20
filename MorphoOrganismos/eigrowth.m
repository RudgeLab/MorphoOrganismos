% w=101; h=101; cx=50; cy=50;
% [x,y] = meshgrid([1:w],[1:h]);
% R = sqrt((x-cx).*(x-cx) + (y-cy).*(y-cy));
% eigrowth(R-20<0, (R-20<0).*1, 0.1, 100)
% eigrowth(R-5<0, (R-5<0).*1, 0.1, 100)

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

%
time = [];          %--> 'time' for plotting (indeed time = loop number)
circ = []; ecc = []; Ar = [];
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
    
    %Get some morpho measures
    phii=(phi(:,:,i)<0);
    
    rp = regionprops(phii, 'all');
    rl = length(rp);

    for j=1:rl;         % --> just in case shape get splitted
        c = 4*pi*rp(j).Area/(rp(j).Perimeter^2);
        circ = [circ, c];       %% --> circle=1, ohter, should be <1.
        ecc = [ecc, rp(j).Eccentricity];  %% circle=0 , ellipse ]0,1[
        MaxAx=rp(j).MajorAxisLength;
        MinAx=rp(j).MinorAxisLength;
        Aratio=MinAx/MaxAx;
        Ar=[Ar,Aratio];
    end;
    time=[time,i-1];
    
    % Plotting

    subplot(2,2,1); imagesc(phi(:,:,i)<0); axis image; colorbar;
    subplot(2,2,2); imagesc(u(:,:,i)); axis image; colorbar;
    subplot(2,2,3); plot(time,circ,time,ecc,time,Ar); xlabel('time');legend('circ','ecc','Asp. ratio','Location','northeastoutside');
    pause(0.1);
    %pause()
        
end;

%save results
mfiles=dir('*.mat');
k=[];

for i=1:length(mfiles)
    k = strfind(mfiles(i).name,'results');
    
    if isempty(k)==0
        fnum(i)=str2num(strtok(strtok(mfiles(i).name,'results')),'.');
    end  
end

if isempty(k)==1
mkdir('Figures');  
savefig([pwd,'\Figures\figure1.fig'])
save('results1.mat','phi','u','circ','ecc','Ar');

else
save(['results',num2str(max(fnum)+1)],'phi','u','circ','ecc','Ar');
savefig([pwd,'\Figures\figure',num2str(max(fnum)+1),'.fig'])
end


    