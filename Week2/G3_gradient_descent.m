function [optimX] = G3_gradient_descent(A,X,b, errorThreshold)

  residual = b - A*X;
  error = norm(residual);
  it = 0;
  while (error > errorThreshold)
    %New learning rate for the iteration
    learnRate = (residual.'*residual) / (residual.' * A * residual);
    %New solution
    X = X + learnRate * residual;
    %New residual and error
    residual = b - A*X;
    error = norm(residual);
    fprintf('iter %03d, error: %.6f\n', it, error);
    it = it + 1;
  end
  
  optimX = X;
  
end