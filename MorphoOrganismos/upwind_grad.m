function [gx,gy] = upwind_grad(igvx, igvy, ims);
    % Compute image (ims) gradient by upwind differencing with respect to
    % velocity (gvx,gvy)

    s = size(ims);
    w = s(1); h = s(2); 
    if length(s)>2;
        d = s(3);
    else;
        d = 1;
    end;
    xidx = 2:w-1;
    yidx = 2:h-1;

    gxpos = zeros(w,h,d);
    gxneg = zeros(w,h,d);
    gypos = zeros(w,h,d);
    gyneg = zeros(w,h,d);

    gxneg(xidx,:,:) = ims(xidx,:,:)-ims(xidx-1,:,:);
    gxpos(xidx,:,:) = ims(xidx+1,:,:)-ims(xidx,:,:);
    gyneg(:,yidx,:) = ims(:,yidx,:)-ims(:,yidx-1,:);
    gypos(:,yidx,:) = ims(:,yidx+1,:)-ims(:,yidx,:);
    
    % Choose upwind direction of image gradient based on velocity in each axis
    gx = 0.5*(1-sign(igvx)).*gxpos + 0.5*(sign(igvx)+1).*gxneg;
    gy = 0.5*(1-sign(igvy)).*gypos + 0.5*(sign(igvy)+1).*gyneg;