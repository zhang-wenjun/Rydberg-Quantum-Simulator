##
using ITensors

##
# In order to define tensors in the ITensors package, one needs to first define index.
# That is, legs of tensors.
local_dim = 2 # The dimension of this leg
sample_index = Index(local_dim, "description") # Define a "leg" with description

# The built-in function `siteinds` produces a set of index of N lattice sites, 
# i.e., "physical legs of tensors" for every sites,
# with well-defined dimension given by the description "S=1/2", and descriptions
N = 2
sites = siteinds("S=1/2",N)
@show sites

##
# The built-in function `randomMPS` creates an MPS on `sites` with bond_dimension
bond_dimension = 2
psi1 = randomMPS(sites, bond_dimension)
# The built-in function `productMPS` creates a product state in the form of MPS
# given by `states`
# Note that since `Julia` is a static language like C++ (not like Python),
# We must specify its data type at definition
states = ["Up" for n=1:N]
psi2 = productMPS(ComplexF64, sites, states)
# We can see an MPS instance contain N tensors, each with one physical leg
# and to bond or link legs, except for the sites at boundary
@show psi1, psi2
# We can access these tensors like an array
@show psi2[1]

##
# We can add two MPSs like ordinary vectors. The built-in `MPS` class overrides the
# operator `add!` (`+`) which lets the MPSs remain in a well-defined orthogonal gauge.
# Note that in `Julia`, functions whose names ended with `!` means they can make
# changes to the input arguments.
psi = (psi1 + psi2) * (0.5^0.5)
@show psi

##
# We can calculate the entanglement entropy between the subsystems (1,2,...,b) 
# and (b+1,...,N) like this
psi = psi2
b = 1
# Move the orthogonality center to `b`
orthogonalize!(psi, b)
# Perform at site `b`
# 
# `linkind(MPS, site_number)` 
# extracts the link or bond Index connecting the MPS tensor on site j to site j+1.
# `siteind(MPS, site_number)` 
# extracts the site Index of the MPS tensor on site j.
# `svd(Tensor, (Index1, Index2, ...))`
# performs svd on Tensor, with (Index1, Index2, ...) grouped together, 
# and all other Index grouped together
U, S, V = svd(psi[b], linkind(psi, b))
SvN = 0.0
for n = 1:dim(S, 1)
    p = S[n,n]^2
    SvN -= p * log2(p)
end
@show SvN

##
# We can also create MPS on our own
print("====================================================\n")
function MPSfromTensor(P::ITensor)
    sites = inds(P)
    N = length(sites)
    tl = Vector{ITensor}(undef, N)
    Pr = copy(P)
    bond = Index(0)
    for s in 1:N-1
        @show s
        if s == 1
            M,S,Pr = svd(Pr, sites[s])
        else
            @show Pr
            @show bond
            M,S,Pr = svd(Pr, (bond, sites[s]))
        end
        temp_link = commonind(Pr, S)
        bond = Index(dim(temp_link), "Link,l=$s")
        @show bond
        tl[s] = M * S * δ(bond, temp_link)
        Pr *= δ(bond, temp_link)
    end
    tl[N] = Pr
    return MPS(tl)
end

N = 3
sites = siteinds("S=1/2",N)
P = ITensor(sites[1], sites[2], sites[3])
P[sites[1]=>1, sites[2]=>1, sites[3]=>1] = P[sites[1]=>2, sites[2]=>2, sites[3]=>2] = 0.5^0.5

psi = MPSfromTensor(P)

b = 1
orthogonalize!(psi, b)
U,S,V = svd(psi[b], (linkind(psi, b), siteind(psi,b+1)))
SvN = 0.0
for n = 1:dim(S, 1)
p = S[n,n]^2
SvN -= p * log2(p)
end

@show SvN