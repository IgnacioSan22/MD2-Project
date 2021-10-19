clearvars;
% dst = double(imread('Images/rainbow_dst.png'));
% src = double(imread('Images/rainbow_src.png'));
dst = double(imread('Images/lena.png'));
src = double(imread('Images/girl.png'));
[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;


%masks to exchange: Eyes
% mask_src=logical(imread('Images/rainbow_src_mask.png'));
% mask_dst=logical(imread('Images/rainbow_dst_mask.png'));
mask_src=logical(imread('Images/mask_src_mouth.png'));
mask_dst=logical(imread('Images/mask_dst_mouth.png'));

for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
    Src_Grad_i = G3_DiFwd(src(:,:,nC), param.hi);
    Src_Grad_j = G3_DjFwd(src(:,:,nC), param.hj);

    driving_Src_Grad_i = zeros(size(src(:,:,1)));
    driving_Src_Grad_i(mask_dst(:)) = Src_Grad_i(mask_src(:));
    driving_Src_Grad_j = zeros(size(src(:,:,1)));
    driving_Src_Grad_j(mask_dst(:)) = Src_Grad_j(mask_src(:));
    magnitudeSrc = sqrt(driving_Src_Grad_i.^2 + driving_Src_Grad_j.^2);

    Dst_Grad_i = G3_DiFwd(dst(:,:,nC), param.hi);
    Dst_Grad_j = G3_DjFwd(dst(:,:,nC), param.hj);

    driving_Dst_Grad_i = zeros(size(src(:,:,1)));
    driving_Dst_Grad_i(mask_dst(:)) = Dst_Grad_i(mask_dst(:));
    driving_Dst_Grad_j = zeros(size(src(:,:,1)));
    driving_Dst_Grad_j(mask_dst(:)) = Dst_Grad_j(mask_dst(:));
    magnitudeDst = sqrt(driving_Dst_Grad_i.^2 + driving_Dst_Grad_j.^2);
    
    driving_Src_Grad_i(magnitudeDst > magnitudeSrc) = driving_Dst_Grad_i(magnitudeDst > magnitudeSrc);
    driving_Src_Grad_j(magnitudeDst > magnitudeSrc) = driving_Dst_Grad_j(magnitudeDst > magnitudeSrc);    
    
    driving_on_dst = zeros(size(src(:,:,1))); 

    driving_on_dst = - G3_DiBwd(driving_Src_Grad_i(:,:), param.hi) - G3_DjBwd(driving_Src_Grad_j(:,:), param.hj);
    param.driving = driving_on_dst;

    dst1(:,:,nC) = G3_Poisson_Equation_Axb(dst(:,:,nC), mask_dst,  param);
end

imshow(dst1/256)