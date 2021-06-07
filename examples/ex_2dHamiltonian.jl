using ITensors
using ColorSchemes
##
include("latticemodels.jl")

##
Ny = 6
Nx = 4

N = Nx*Ny

sites = siteinds("S=1/2", N)

# Obtain a LatticeModel2D struct
# which defines a lattice
latt = square_lattice_wjz(Nx, Ny; yperiodic = false)

# Define the Ising spin Hamiltonian on this lattice
ampo = AutoMPO()
for b in latt.locterms
    ampo .+= "Sz", b.s[1], "Sz", b.s[2]
end
H = MPO(ampo, sites)

##

# Initialize wavefunction to a random MPS
# of bond-dimension 10
psi0 = randomMPS(sites, 10)

sweeps = Sweeps(10)
maxdim!(sweeps,20,60,100,100,200,400,800)
cutoff!(sweeps,1E-8)
@show sweeps

energy,psi = dmrg(H,psi0,sweeps)

##
function meas_Sz(psi, n)
    psi = orthogonalize(psi,n)
    sn = siteind(psi, n)
    Sz = scalar(dag(prime(psi[n],"Site"))*op("Sz",sn)*psi[n])
    return real(Sz)
end

Sz_list = [meas_Sz(psi, n) for n = 1:N]
plot(latt, Sz_list)