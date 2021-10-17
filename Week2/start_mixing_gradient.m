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
    drivingGrad_i = G3_DiBwd(src(:,:,nC), param.hi) - G3_DiFwd(src(:,:,nC), param.hi);
    drivingGrad_j = G3_DjBwd(src(:,:,nC), param.hj) - G3_DjFwd(src(:,:,nC), param.hj);

    driving_on_src = drivingGrad_i + drivingGrad_j;

    driving_on_dst = zeros(size(src(:,:,1))); 

    drivingGrad_i_dst = G3_DiBwd(dst(:,:,nC), param.hi) - G3_DiFwd(dst(:,:,nC), param.hi);
    drivingGrad_j_dst = G3_DjBwd(dst(:,:,nC), param.hj) - G3_DjFwd(dst(:,:,nC), param.hj);
    driving_on_dst_aux = drivingGrad_i_dst + drivingGrad_j_dst;
        
    driving_mix = cat(3, driving_on_src(mask_src(:)), driving_on_dst_aux(mask_dst(:)));
    [~, idx] = max(abs(driving_mix), [], 3);
        
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_src(:)).*(idx==1) + driving_on_dst_aux(mask_dst(:)).*(idx==2);
    
    param.driving = driving_on_dst;

    dst1(:,:,nC) = G3_Poisson_Equation_Axb(dst(:,:,nC), mask_dst,  param);
end

imshow(dst1/256)