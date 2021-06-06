using ITensors

i = Index(4,"i")
j = Index(4,"j")
k = Index(4,"k")
l = Index(4,"l")
T = randomITensor(i,j,k,l)
U,S,V = svd(T,i,k)   # compute SVD with (i,k) as row indices (indices of U)
@show hasinds(U,i,k) # = true
@show hasinds(V,j,l) # = true
@show T ≈ U*S*V      # = true
@show hasinds(S,j)
