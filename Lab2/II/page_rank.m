function [numer_indeksu, Edges, I, B, A, b, r] = page_rank()
numer_indeksu = 193249;
%L1 = 4         mod(L1, 7)+1, = 4%7+1=4+1=5
%L2 = 2         mod(L2, 7)+1  = 2%7+1=2+1=3

Edges = [1,1,2,2,2,3,3,3,3,4,4,5,5,6,6,7,8;
         4,6,3,4,5,5,6,7,8,5,6,4,6,4,7,6,5];
I = speye(8);

row = Edges(2,:);
column = Edges(1,:);
B = sparse(row,column, 1);

N = size(B, 2);

sum_column = sum(B);
A = spdiags(1./sum_column', 0, N, N);

d = 0.85;
b = ((1 - d) / N) * ones(N, 1);

M = I-d*B*A;
r_vector = M\b;
r = reshape(r_vector, [], 1);
end