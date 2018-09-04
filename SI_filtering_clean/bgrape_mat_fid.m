function [ fid ] = bgrape_mat_fid( mat1, mat2, subspace_dim )

mat1_conj = ctranspose(mat1);

fid = (1/subspace_dim^2) * abs(bgrape_trace_matmul(mat1_conj,mat2))^2;

end

