clearvars;
dst = double(imread('Images/lena.png'));
src = double(imread('Images/girl.png')); % flipped girl, because of the eyes
[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;


%masks to exchange: Eyes
mask_src=logical(imread('Images/mask_src_eyes.png'));
mask_dst=logical(imread('Images/mask_dst_eyes.png'));

for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
    drivingGrad_i = G3_DiBwd(src(:,:,nC), param.hi) - G3_DiFwd(src(:,:,nC), param.hi);
    drivingGrad_j = G3_DjBwd(src(:,:,nC), param.hj) - G3_DjFwd(src(:,:,nC), param.hj);

    driving_on_src = drivingGrad_i + drivingGrad_j;
    
    driving_on_dst = zeros(size(src(:,:,1)));   
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_src(:));
    
    param.driving = driving_on_dst;

    dst1(:,:,nC) = G3_Poisson_Equation_Axb(dst(:,:,nC), mask_dst,  param, false);
end

%Mouth
%masks to exchange: Mouth
mask_src=logical(imread('Images/mask_src_mouth.png'));
mask_dst=logical(imread('Images/mask_dst_mouth.png'));
for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
    %Compute the gradients with a forward and backward pass of the kernel
    drivingGrad_i = G3_DiBwd(src(:,:,nC), param.hi) - G3_DiFwd(src(:,:,nC), param.hi);
    drivingGrad_j = G3_DjBwd(src(:,:,nC), param.hj) - G3_DjFwd(src(:,:,nC), param.hj);

    driving_on_src = drivingGrad_i + drivingGrad_j;
    
    driving_on_dst = zeros(size(src(:,:,1)));
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_src(:));
    
    param.driving = driving_on_dst;

    %Perform poisson editing
    dst1(:,:,nC) = G3_Poisson_Equation_Axb(dst1(:,:,nC), mask_dst,  param, false);
end

imshow(dst1/256)