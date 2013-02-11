function N = vectnorm (M)
% Take L2 norm of all row vectors of matrix
  N = sqrt(sum(abs(M).^2,2));
end

