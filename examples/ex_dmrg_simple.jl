#############################################
# A simple DMRG code for pedagogical reason #
#                                           #
# I did NOT use the `ITensors` built-in     #
# DMRG related functions, but utilized its  #
# functions for tensors, MPS, and MPO.      #
# For the usage of these, see `ex_MPS.jl`   #
# and `ex_MPO.jl`                           #
#                                           #
# Created by: Wenjun Zhang                  #
#             zhangwenjun1998@gmail.com     #
#             University of Chicago         #
#         on: June 18, 2021                 #
# Modified:                                 #
#   June 20, 2021                           #
#     - Debugged the local eigensolver      #
#############################################


# Import packages
using ITensors
using Plots
##
# Define lattice
N = 10
# For the usage of `siteinds` and MPS related functions,
# see `examples/ex_MPS.jl`
sites = siteinds("S=1/2", N)
# Define Hamiltonian
# This example Hamiltonian describes a Rydberg system
# whose ground state has a Z3 order
Omega = 2π * 2
Delta = 4*Omega
Rb = 2.432
# For the usage of `AutoMPO`, see `examples/ex_MPO.jl`
ampo = AutoMPO()
# Only keep the `maxnnorder` nearest neighbor coupling
# E.g., setting `maxnnorder=1` keeps the nearest neighbor
# coupling and the next-nearest neighbor coupling
maxnnorder = 4
for j = 1:N
    # Single-qubit terms
    ampo .+=  Omega,     "Sx", j
    ampo .+= -Delta, "ProjUp", j
    # Interaction terms
    for i = j+1:min(j+maxnnorder, N)
        Vij = Omega*(Rb/(i - j))^6
        ampo .+= Vij, "ProjUp", j, "ProjUp", i
    end
end
H = MPO(ampo, sites)
# Define initial state
psi0 = randomMPS(sites, 2)
# Set sweep parameters
nsweep = 10    # number of sweeps
maxdim = 200   # maximun bond dimension kept during sweeps
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
    sw_time = @elapsed begin
    maxtruncerr = 0.0
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
    println("After sweep $sw energy=$energy maxlinkdim=$maxlinkdimpsi maxerr=$maxtruncerr time=$sw_time")
    # println("After sweep $sw energy=$energy maxlinkdim=$maxlinkdimpsi maxerr=$maxtruncerr")
end

##
function meas_Sz(psi, n)
    psi = orthogonalize(psi,n)
    sn = siteind(psi, n)
    Sz = scalar(dag(prime(psi[n],"Site"))*op("Sz",sn)*psi[n])
    return real(Sz)
end

Sz_list = [meas_Sz(psi, n) for n = 1:N]
plot(1:N, Sz_list)