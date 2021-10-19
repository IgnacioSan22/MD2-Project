function [u] = G3_Poisson_Equation_Axb(f, dom2Inp, param)
%this code is not intended to be efficient. 

[ni, nj] = size(f);

%We add the ghost boundaries (for the boundary conditions)
f_ext = zeros(ni+2, nj+2);
f_ext(2:end-1, 2:end-1) = f;
dom2Inp_ext = zeros(ni+2, nj+2);
dom2Inp_ext(2:end-1, 2:end-1) = dom2Inp;

%Store memory for the A matrix and the b vector    
nPixels = (ni+2) * (nj+2); %Number of pixels

%We will create A sparse, this is the number of nonzero positions
% 4 for each pixel in the mask
% 2 times row and column size, beacuse ghost boundaries only use one neighbour
nonZeroPos = nPixels + 4*nnz(dom2Inp) + 2*(ni+2) + 2*(nj+2);

%idx_Ai: Vector for the nonZero i index of matrix A
%idx_Aj: Vector for the nonZero j index of matrix A
%a_ij: Vector for the value at position ij of matrix A
%Initialize vectors to fill later in the loops
idx_Ai = zeros(nonZeroPos,1);
idx_Aj = zeros(nonZeroPos,1);
a_ij = zeros(nonZeroPos,1);

b = zeros(nPixels,1);

%Vector counter
idx = 1;

%North side boundary conditions
i = 1;
for j=1:nj+2
    %from image matrix (i,j) coordinates to vectorial (p) coordinate
    p = (j-1)*(ni+2)+i;
    
    
    %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
    %vector b
    idx_Ai(idx)=p; 
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx = idx+1;
    
    idx_Ai(idx) = p;
    idx_Aj(idx) = p+1;
    a_ij(idx) = -1;   
    idx = idx+1;
            
    b(p) = 0;
end

%South side boundary conditions
i = ni+2;
for j=1:nj+2
    %from image matrix (i,j) coordinates to vectorial (p) coordinate
    p = (j-1)*(ni+2)+i;
    
    %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
    %vector b
    %TO COMPLETE 2
    idx_Ai(idx) = p; 
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx = idx+1;
    
    idx_Ai(idx) = p;
    idx_Aj(idx) = p-1;
    a_ij(idx) = -1;   
    idx = idx+1;
            
    b(p) = 0;
    
end

%West side boundary conditions
j = 1;
for i=1:ni+2
    %from image matrix (i,j) coordinates to vectorial (p) coordinate
    p = (j-1)*(ni+2)+i;

    %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
    %vector b
    %TO COMPLETE 3
    idx_Ai(idx) = p; 
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx = idx+1;
    
    idx_Ai(idx) = p;
    %Move the vertical size of the image to change to the next column,
    % the image displayed as vector is natural rowwise ordered
    idx_Aj(idx) = p + (ni+2);
    a_ij(idx) = -1;   
    idx = idx+1;
            
    b(p) = 0;
end

%East side boundary conditions
j = nj+2;
for i=1:ni+2
    %from image matrix (i,j) coordinates to vectorial (p) coordinate
    p = (j-1)*(ni+2)+i;
    
    %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
    %vector b
    %TO COMPLETE 4
    idx_Ai(idx) = p; 
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx = idx+1;
    
    idx_Ai(idx) = p;
    %Move the vertical size of the image to change to the previous column,
    % the image displayed as vector is natural rowwise ordered
    idx_Aj(idx) = p - (ni+2);
    a_ij(idx) = -1;   
    idx = idx+1;
            
    b(p) = 0;
    
end

%Inner points

for j=2:nj+1
    for i=2:ni+1
     
        %from image matrix (i,j) coordinates to vectorial (p) coordinate
        p = (j-1)*(ni+2)+i;
                                            
        if (dom2Inp_ext(i,j) == 1) %If we have to inpaint this pixel
            % 4V(x,y) - V(x+1,y) - V(x-1,y) - V(x,y+1) - V(x,y-1) = 0
            %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
            %vector b
            %TO COMPLETE 5
            idx_Ai(idx) = p; 
            idx_Aj(idx) = p; 
            a_ij(idx) = 4;
            idx = idx+1;

            idx_Ai(idx) = p;
            idx_Aj(idx) = p+1;
            a_ij(idx) = -1;
            idx = idx+1;

            idx_Ai(idx) = p;
            idx_Aj(idx) = p-1;
            a_ij(idx) = -1;
            idx = idx+1;

            idx_Ai(idx) = p;
            idx_Aj(idx) = p+(ni+2);
            a_ij(idx) = -1;
            idx = idx+1;

            idx_Ai(idx) = p;
            idx_Aj(idx) = p-(ni+2);
            a_ij(idx) = -1;
            idx = idx+1;
                
            %If driving exists, do Poisson editting
            if (isfield(param, 'driving'))
                b(p) = param.driving(i-1,j-1);
            else
                b(p) = 0;
            end
    
        else %we do not have to inpaint this pixel 
        % Just set the original value 𝑉(x,y)=𝑈(x,y) at each (x,y) in A
            %Fill Idx_Ai, idx_Aj and a_ij with the corresponding values and
            %vector b
             %TO COMPLETE 6
            idx_Ai(idx) = p; 
            idx_Aj(idx) = p; 
            a_ij(idx) = 1;
            idx = idx + 1;
            b(p) = f_ext(i, j); 
            
        end       
    end
end

%A is a sparse matrix, so for memory requirements we create a sparse
%matrix
%TO COMPLETE 7
A = sparse(idx_Ai, idx_Aj, a_ij, nPixels, nPixels); %??? and ???? is the size of matrix A

%Solve the system of equations
% x = mldivide(A,b);
x = G3_gradient_descent(A, f_ext(:), b, 0.001);

%From vector to matrix
u_ext = reshape(x, ni+2, nj+2);

%Eliminate the ghost boundaries
u = full(u_ext(2:end-1, 2:end-1));
    
