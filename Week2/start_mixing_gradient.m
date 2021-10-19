clearvars;
dst = double(imread('Images/rainbow_dst.png'));
src = double(imread('Images/rainbow_src.png')); % flipped girl, because of the eyes
[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;


%masks to exchange: Eyes
mask_src=logical(imread('Images/rainbow_src_mask.png'));
mask_dst=logical(imread('Images/rainbow_dst_mask.png'));

for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
%     drivingGrad_i = G3_DiBwd(src(:,:,nC), param.hi) - G3_DiFwd(src(:,:,nC), param.hi);
%     drivingGrad_j = G3_DjBwd(src(:,:,nC), param.hj) - G3_DjFwd(src(:,:,nC), param.hj);
    drivingGrad_i = G3_DiFwd(src(:,:,nC), param.hi);
    drivingGrad_i_aux = zeros(size(src(:,:,1)));
    drivingGrad_i_aux(mask_dst(:)) = drivingGrad_i(mask_src(:));

    drivingGrad_j = G3_DjFwd(src(:,:,nC), param.hj);
    drivingGrad_j_aux = zeros(size(src(:,:,1)));
    drivingGrad_j_aux(mask_dst(:)) = drivingGrad_j(mask_src(:));

    magnitudeSrc = sqrt(drivingGrad_i.^2+drivingGrad_j.^2);
    magnitudeSrcAux = zeros(size(src(:,:,1)));
    magnitudeSrcAux(mask_dst(:)) = magnitudeSrc(mask_src(:));
    driving_on_src = drivingGrad_i + drivingGrad_j;

    driving_on_dst = zeros(size(src(:,:,1))); 

%     drivingGrad_i_dst = G3_DiBwd(dst(:,:,nC), param.hi) - G3_DiFwd(dst(:,:,nC), param.hi);
%     drivingGrad_j_dst = G3_DjBwd(dst(:,:,nC), param.hj) - G3_DjFwd(dst(:,:,nC), param.hj);
    drivingGrad_i_dst = G3_DiFwd(dst(:,:,nC), param.hi);
    drivingGrad_j_dst = G3_DjFwd(dst(:,:,nC), param.hj);
    magnitudeDst = sqrt(drivingGrad_i_dst.^2+drivingGrad_j_dst.^2);
    magnitudeDstAux = zeros(size(src(:,:,1)));
    magnitudeDstAux(mask_dst(:)) = magnitudeSrc(mask_dst(:));
    driving_on_dst_aux = drivingGrad_i_dst + drivingGrad_j_dst;
        
%     driving_mix = cat(3, driving_on_src(mask_src(:)), driving_on_dst_aux(mask_dst(:)));
%     [~, idx] = max(abs(driving_mix), [], 3);
    
    drivingGrad_i_dst(magnitudeSrcAux > magnitudeDstAux) = drivingGrad_i_aux(magnitudeSrcAux > magnitudeDstAux);
    drivingGrad_j_dst(magnitudeSrcAux > magnitudeDstAux) = drivingGrad_j_aux(magnitudeSrcAux > magnitudeDstAux);
%     driving_on_dst(mask_dst(:)) = driving_on_src(mask_src(:)).*(idx==1) + driving_on_dst_aux(mask_dst(:)).*(idx==2);
    
    driving_on_dst = G3_DiBwd(drivingGrad_i_dst(:,:), param.hi) + G3_DjBwd(drivingGrad_j_dst(:,:), param.hj);
%     driving_on_dst = drivingGrad_i_dst + drivingGrad_j_dst;
    param.driving = driving_on_dst;

    dst1(:,:,nC) = G3_Poisson_Equation_Axb(dst(:,:,nC), mask_dst,  param);
end

imshow(dst1/256)