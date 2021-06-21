#############################################
# A simple DMRG code for pedagogical reason #
# 2D version                                #
# I did NOT use the `ITensors` built-in     #
# DMRG related functions, but utilized its  #
# functions for tensors, MPS, and MPO.      #
# For the usage of these, see `ex_MPS.jl`   #
# and `ex_MPO.jl`                           #
#                                           #
# Created by: Wenjun Zhang                  #
#             zhangwenjun1998@gmail.com     #
#             University of Chicago         #
#         on: June 20, 2021                 #
#############################################


# Import packages
using ITensors
using Plots
using ColorSchemes
include("G:/我的云端硬盘/Grad/Rydberg-Quantum-Simulator/latticemodels.jl")
##
# Define lattice
Nx = 5
Ny = 5
N = Nx*Ny
max_nnorder = 2
sites = siteinds("S=1/2", N)
latt = square_lr(Nx, Ny; max_nnorder=max_nnorder)
# Define Hamiltonian
# This example Hamiltonian describes a Rydberg system
# whose ground state has a striated order
Omega = 2π * 4.3
Delta = 2*Omega
Rb = 1.6
ampo = AutoMPO()
# Single-qubit terms
for j = 1:N
    ampo .+=  Omega,     "Sx", j
    ampo .+= -Delta, "ProjUp", j
end
# Interaction terms
for cplt in latt.cplterms
    cplrange = parse(Float64, cplt.type[5:end])
    s1 = cplt.s[1]
    s2 = cplt.s[2]
    Vij = Omega*(Rb/cplrange)^6
    ampo .+= Vij, "ProjUp", s1, "ProjUp", s2
end
H = MPO(ampo,sites)
# Define initial state
psi0 = randomMPS(sites, 2)
# Set sweep parameters
nsweep = 10    # number of sweeps
maxdim = 200   # maximum bond dimension kept during sweeps
cutoff = 1E-10 # cutoff for Schmidt coefficients (singular values)
# Store the contraction for projected local 
# effective Hamiltonians
# 1  2 ... i i+1            N-1 N  
# o--o--o-      -o--o--o--o--o--o <psi|
# |  |  |  |  |  |  |  |  |  |  |
# o--o--o--o--o--o--o--o--o--o--o H
# |  |  |  |  |  |  |  |  |  |  |
# o--o--o-      -o--o--o--o--o--o |psi>
# 
# TNs[i] with i = 1,2,...,lpos
# stores the network
# 1 ... i
# o--o--o-- 
# |  |  |  
# o--o--o--
# |  |  |  
# o--o--o--
# TNs[i] with i = rpos,...,N
# stores the network
#   i i+1 ...  N-1 N
# --o--o--o--o--o--o
#   |  |  |  |  |  |
# --o--o--o--o--o--o
#   |  |  |  |  |  |
# --o--o--o--o--o--o
# 
# lpos and rpos denotes the parts of networks
# which have already been contracted
lpos = 0
rpos = N + 1
TNs = Vector{ITensor}(undef, N)

# Begin sweeps
# copy the initial state
psi = copy(psi0)
# and move the canonical center to the first site
orthogonalize!(psi, 1)
energy = 0.0
for sw = 1:nsweep
    println("Begin sweep $sw:")
    sw_time = @elapsed begin
    maxtruncerr = 0.0
    maxlocaldim = 0
    # Sweep towards the right
    for b in 1:N-1
        # Get local effective Hamiltonian
        # Contract the network to the left
        while lpos < b - 1
            if lpos < 1
                TNs[1] = psi[1]*H[1]*dag(prime(psi[1]))
                lpos = 1
            else
                TNs[lpos+1] = TNs[lpos]*
                    psi[lpos+1]*H[lpos+1]*dag(prime(psi[lpos+1]))
                lpos += 1
            end
        end
        # Contract the network to the right
        while rpos > b + 2
            if rpos > N
                TNs[N] = psi[N]*H[N]*dag(prime(psi[N]))
                rpos = N
            else
                TNs[rpos-1] = TNs[rpos]*
                    psi[rpos-1]*H[rpos-1]*dag(prime(psi[rpos-1]))
                rpos -= 1
            end
        end
        lpos = b - 1
        rpos = b + 2
        Heff = H[b] * H[b+1]
        b > 1 && (Heff *= TNs[b-1])
        b+1 < N && (Heff *= TNs[b+2])
        # Solve local eigen problem
        f = eigen(Heff; ishermitian = true)
        # f.D is the diagonalized matrix
        # f.V is the right-eigenvectors
        # f.Vt is the left-eigenvectors
        # H*f.V≈f.D*f.Vt or dag(f.V)*f.D*f.Vt≈H
        # i--H--i'
        # ==
        # i--dag(f.V)--ie ~ ie--f.D--ie' ~ ie'--f.Vt--i'
        # where ~ means contraction
        ie = commonind(f.D, f.V)
        edim = dim(ie)
        maxlocaldim = max(maxlocaldim, edim)
        energy = f.D[ie=>edim,ie'=>edim]
        ext = ITensor(ie)
        ext[ie=>edim] = 1
        phi = f.V * ext
        # Decompose the two-site tensor phi
        # To put the canonical center on the 
        # left tensor or the right tensor
        # When sweeping towards right, we put 
        # it on the right tensor, which means 
        # the left tensor should be in orthogonal basis.
        ortho = "left" 
        indices = inds(psi[b])
        # Rename the link index for future convenience 
        # and compatibility with `ITensors` built-in functions
        L, R, spec = factorize(phi, indices;
                        tags = tags(linkind(psi, b)), 
                        ortho = ortho,
                        cutoff = cutoff,
                        maxdim = maxdim)
        psi[b] = L
        psi[b+1] = R

        # Truncation error
        maxtruncerr = max(maxtruncerr, spec.truncerr)
    end
    # Sweep towards the right
    for b in N-1:-1:1
        # Get local effective Hamiltonian
        # Contract the network to the left
        while lpos < b - 1
            if lpos < 1
                TNs[1] = psi[1]*H[1]*dag(prime(psi[1]))
                lpos = 1
            else
                TNs[lpos+1] = TNs[lpos]*
                    psi[lpos+1]*H[lpos+1]*dag(prime(psi[lpos+1]))
                lpos += 1
            end
        end
        # Contract the network to the right
        while rpos > b + 2
            if rpos > N
                TNs[N] = psi[N]*H[N]*dag(prime(psi[N]))
                rpos = N
            else
                TNs[rpos-1] = TNs[rpos]*
                    psi[rpos-1]*H[rpos-1]*dag(prime(psi[rpos-1]))
                rpos -= 1
            end
        end
        lpos = b - 1
        rpos = b + 2
        Heff = H[b] * H[b+1]
        b > 1 && (Heff *= TNs[b-1])
        b+1 < N && (Heff *= TNs[b+2])
        # Solve local eigen problem
        f = eigen(Heff; ishermitian = true)
        # f.D is the diagonalized matrix
        # f.V is the right-eigenvectors
        # f.Vt is the left-eigenvectors
        # H*f.V≈f.D*f.Vt or dag(f.V)*f.D*f.Vt≈H
        # i--H--i'
        # ==
        # i--dag(f.V)--ie ~ ie--f.D--ie' ~ ie'--f.Vt--i'
        # where ~ means contraction
        ie = commonind(f.D, f.V)
        edim = dim(ie)
        maxlocaldim = max(maxlocaldim, edim)
        energy = f.D[ie=>edim,ie'=>edim]
        ext = ITensor(ie)
        ext[ie=>edim] = 1
        phi = f.V * ext
        # Decompose the two-site tensor phi
        ortho = "right"
        indices = inds(psi[b])
        L, R, spec = factorize(phi, indices;
                        tags = tags(linkind(psi, b)),
                        ortho = ortho,
                        cutoff = cutoff,
                        maxdim = maxdim)
        psi[b] = L
        psi[b+1] = R

        # Truncation error
        maxtruncerr = max(maxtruncerr, spec.truncerr)
    end
    end
    maxlinkdimpsi = maxlinkdim(psi)
    println("    energy=$energy maxlinkdim=$maxlinkdimpsi maxerr=$maxtruncerr time=$sw_time")
    println("    maximum local dimension: $maxlocaldim")
    # println("After sweep $sw energy=$energy maxlinkdim=$maxlinkdimpsi maxerr=$maxtruncerr")
end

##
function meas_dens(psi, s)
    psi = orthogonalize(psi, s)
    site = siteind(psi, s)
    n = scalar(dag(prime(psi[s], "Site"))*op("ProjUp", site)*psi[s])
    return real(n)
end

dens_list = [meas_dens(psi, n) for n = 1:N]
plot(latt, dens_list)