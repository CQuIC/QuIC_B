function c = bgrape_trace_matmul(A,B)
% taken from DYNAMO package
% 
% Computes trace(A*B) efficiently
% Utilizes the identity: trace(A*B) == sum(sum(transpose(A).*B)
% left side is O(n^3) to compute, right side is O(n^2)

c = sum(sum((A.') .* B));
end
