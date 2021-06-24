module DemoSPCA

include("./SPCA.jl")

using .SparsePCA
using LinearAlgebra
using Printf: @printf
using DelimitedFiles: readdlm

maxiter = 1e4
tol = 1e-5
n = 4  # columns
d = 500  # dimension
m = 1000  # sample size
mu = 0.1 * ones(n, 1)
type = 0
lambda = 1
f_palm = 1e5

a = readdlm("A.csv", ',', Float64, '\n')
_, _, v = svd(a)
x0 = v[:, 1:n]

results = spca_amanpg(a, mu, lambda, n, x0, x0, f_palm, verbose=true)
# iter, f_amanpg, sparsity, timediff, x, y_man = spca_amanpg(a, mu, lambda, n, x0, x0, f_palm, verbose=true)
# printf("%i iterations with final value %.4f, sparsity %.4f, timediff %.4f.", iter, f_amanpg, sparsity, time)

end