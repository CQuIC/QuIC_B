function [ outFid,uniEigFid,uniFid ] = checkUniFid_fn( uni_final, uni_target )
%CHECKSTATEFID_FN Summary of this function goes here
%   Detailed explanation goes here

%load(filename);

% uni_final=opt_params.uni_final;
% uni_target=opt_params.target_uni;

%Find ordered eigenvectors
[final_vectors,final_values]=eig(uni_final);
[target_vectors,target_values]=eig(uni_target);

%Order the eigenvalues by phase
phase_final=zeros(1,16);
phase_target=zeros(1,16);
for i=1:1:16
    phase_final(1,i)=phase(final_values(i,i));
    phase_target(1,i)=phase(target_values(i,i));
end

%Order the matrices containing eigenvectors and eigenvalues
[sorted_final,index_f]=sort(phase_final);
[sorted_target,index_t]=sort(phase_target);

final_vectors=final_vectors(:,index_f);
target_vectors=target_vectors(:,index_t);

sorted_final_values=final_values(index_f,index_f);
sorted_target_values=target_values(index_t,index_t);

%Create map between each unitary and its representation in its
%eigen-basis
A=zeros(16);
B=zeros(16);
I=eye(16);
for j=1:1:16
    for i=1:1:16
        A(j,i)=ctranspose(target_vectors(:,j))*I(:,i);
        B(j,i)=ctranspose(final_vectors(:,j))*I(:,i);
    end
end

%Create the map between bases
eigen_map=ctranspose(B)*A;

% mapped_uni_final should look like the target unitary
mapped_uni_final=ctranspose(eigen_map)*ctranspose(B)*sorted_final_values*B*eigen_map;

% Grab the unitary fidelity
uniEigFid=((1/16)^2)*abs(trace(ctranspose(sorted_final_values)*sorted_target_values))^2;

% Calculate the cost function
costFun = 0;
for l=1:1:16
    costFun = costFun + abs(sorted_final_values(l,l)-sorted_target_values(l,l))^2;
end

% Test whether the fidelity of a state v1 after evolution is above 0.995
% v1=rand(16,1);
% v1=v1/norm(v1);
% 
% (mapped_uni_final*v1)-(uni_target*v1);
% 
% outFid=(norm(ctranspose(mapped_uni_final*v1)*(uni_target*v1)))^2;

% Calculate the mapped unitary fidelity (for comparison)
outFid = ((1/16)^2)*abs(trace(ctranspose(mapped_uni_final)*uni_target))^2;

uniFid = ((1/16)^2)*abs(trace(ctranspose(uni_final)*uni_target))^2;

end

