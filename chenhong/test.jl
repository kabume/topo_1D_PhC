
using SharedArrays
using Distributed

addprocs(36)

A = SharedArray{Float64}(2000, 2000);
rr = CartesianIndices(size(A));
@distributed for k in rr
    A[k] = 1;
end
