function [ phi ] = G3_ChanVeseIpol_GDExp( I, phi_0, mu, nu, eta, lambda1, lambda2, tol, epHeaviside, dt, iterMax, reIni )
%Implementation of the Chan-Vese segmentation following the explicit
%gradient descent in the paper of Pascal Getreur "Chan-Vese Segmentation".
%It is the equation 19 from that paper

%I     : Gray color image to segment
%phi_0 : Initial phi
%mu    : mu lenght parameter (regularizer term)
%nu    : nu area parameter (regularizer term)
%eta   : epsilon for the total variation regularization
%lambda1, lambda2: data fidelity parameters
%tol   : tolerance for the sopping criterium
% epHeaviside: epsilon for the regularized heaviside. 
% dt     : time step
%iterMax : MAximum number of iterations
%reIni   : Iterations for reinitialization. 0 means no reinitializacion

[ni,nj]=size(I);
hi=1;
hj=1;

filename = 'noised_cone.gif';
h = figure;
phi=phi_0;
phi_extend = zeros(ni+2, nj+2);
dif=inf;
nIter=0;
while dif>tol && nIter<iterMax
    
    phi_old=phi;
    nIter=nIter+1;        
    
    
    %Fixed phi, Minimization w.r.t c1 and c2 (constant estimation)
    H = 0.5 * (1 + (2/pi) .* atan(phi ./ epHeaviside));
%     H = (phi > 0);
%     Hout = phi <= 0;
    c1 = sum(I .* H,"all") ./ sum(H,"all"); %TODO 1: Line to complete
    c2 = sum(I .* (1-H),"all") ./ sum(1-H,"all"); %TODO 2: Line to complete
%     c1 = sum(I .* H,"all") ./ sum(H,"all"); %TODO 1: Line to complete
%     c2 = sum(I .* Hout,"all") ./ sum(Hout,"all"); %TODO 2: Line to complete

    %Boundary conditions
%     phi(:,1) = phi(:,2);
%     phi(:,end) = phi(:,end-1);
%     phi(1,:) = phi(2,:);
%     phi(end,:) = phi(end-1,:);
    phi_extend(2:end-1,2:end-1) = phi;
    phi_extend(1,:)   = phi_extend(end-1,:); %TODO 3: Line to complete
    phi_extend(end,:) = phi_extend(2,:); %TODO 4: Line to complete

    phi_extend(:,1)   = phi_extend(:,end-1); %TODO 5: Line to complete
    phi_extend(:,end) = phi_extend(:,2); %TODO 6: Line to complete

    
    %Regularized Dirac's Delta computation
    delta_phi = sol_diracReg(phi, epHeaviside);   %notice delta_phi=H'(phi)	
    
    %derivatives estimation
    %i direction, forward finite differences
    phi_iFwd  = DiFwd(phi_extend(:,:), hi); %TODO 7: Line to complete
    phi_iBwd  = DiBwd(phi_extend(:,:), hi); %TODO 8: Line to complete
    
    %j direction, forward finitie differences
    phi_jFwd  = DjFwd(phi_extend(:,:), hj); %TODO 9: Line to complete
    phi_jBwd  = DjBwd(phi_extend(:,:), hj); %TODO 10: Line to complete
    
    %centered finite diferences
    phi_icent   = (phi_iFwd + phi_iBwd) / 2; %TODO 11: Line to complete
    phi_jcent   = (phi_jFwd + phi_jBwd) / 2; %TODO 12: Line to complete
    
    %A and B estimation (A y B from the Pascal Getreuer's IPOL paper "Chan
    %Vese segmentation
    A = mu ./ sqrt(eta^2 + (phi_iFwd).^2 + (phi_jcent).^2); %TODO 13: Line to complete
    B = mu ./ sqrt(eta^2 + (phi_icent).^2 + (phi_jFwd).^2); %TODO 14: Line to complete


    %%Equation 22, for inner points
    A1 = A(2:end-1,2:end-1) .* phi_extend(3:end,2:end-1); 
    A2 = A(1:end-2,2:end-1) .* phi_extend(1:end-2,2:end-1);
    B1 = B(2:end-1,2:end-1) .* phi_extend(2:end-1,3:end);
    B2 = B(2:end-1,1:end-2) .* phi_extend(2:end-1,1:end-2);
%     A1 = A(2:end-1,2:end-1) .* phi(3:end,2:end-1);
%     A2 = A(1:end-2,2:end-1) .* phi(1:end-2,2:end-1);
%     B1 = B(2:end-1,2:end-1) .* phi(2:end-1,3:end);
%     B2 = B(2:end-1,1:end-2) .* phi(2:end-1,1:end-2);

%     c1Factor = lambda1 * (I(2:end-1,2:end-1) - c1).^2;
%     c2Factor = lambda2 * (I(2:end-1,2:end-1) - c2).^2;
    c1Factor = lambda1 * (I - c1).^2;
    c2Factor = lambda2 * (I - c2).^2;
    divisionFactor = A(2:end-1,2:end-1) + A(1:end-2,2:end-1) + ...
                     B(2:end-1,2:end-1) + B(2:end-1,1:end-2);
    
    phi = (phi_old + dt .* delta_phi .* (A1 + A2 + B1 + B2 - nu - c1Factor + c2Factor) ) ...
        ./ (1 + dt * delta_phi .* divisionFactor); 
    %TODO 15: Line to complete
    
    %Reinitialization of phi
    if reIni>0 && mod(nIter, reIni)==0
        indGT = phi >= 0;
        indLT = phi < 0;
        
        phi=double(bwdist(indLT) - bwdist(indGT));
        
        %Normalization [-1 1]
        nor = min(abs(min(phi(:))), max(phi(:)));
        phi=phi/nor;
    end
  
    %Diference. This stopping criterium has the problem that phi can
    %change, but not the zero level set, that it really is what we are
    %looking for.
%     dif = sqrt(mean(sum( (phi(:) - phi_old(:)).^2 )));
    dif = sqrt(sum( (phi(:) - phi_old(:)).^2 ) / (ni*nj));
         
   if mod(nIter, 25) == 0 || nIter == 1
        
        %Plot the level sets surface
        subplot(1,2,1) 
        %The level set function
        surfc(phi,'LineStyle','-')  %TODO 16: Line to complete 
        hold on
        %The zero level set over the surface
%         contour(phi, 'LineColor', 'blue'); %TODO 17: Line to complete
        hold off
        title('Phi Function');
        
        %Plot the curve evolution over the image
        subplot(1,2,2)
        imagesc(I);        
        colormap gray;
        hold on;

        contour(phi > 0, 'LineColor', 'red') %TODO 18: Line to complete
        title(nIter)
    
        axis off;
        hold off
        drawnow;
        pause(.0001); 
        fprintf('iter %03d, error: %.6f\n', nIter, dif);
        % Capture the plot as an image 
          frame = getframe(h); 
          im = frame2im(frame); 
          [imind,cm] = rgb2ind(im,256); 
          % Write to the GIF File 
          if nIter == 1 
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
          else 
              imwrite(imind,cm,filename,'gif','WriteMode','append'); 
          end 
    end
end

phi = phi >= 0;