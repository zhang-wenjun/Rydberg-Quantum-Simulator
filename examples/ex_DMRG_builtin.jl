##
using ITensors
using Plots
##

N = 10
sites = siteinds("S=1/2", N)
ampo = AutoMPO()
for j = 1:N-1
    ampo .+= "Sz",j,"Sz",j+1
end
H = MPO(ampo,sites)

sweeps = Sweeps(10) # number of sweeps is 5
maxdim!(sweeps,10,20,100,100,200) # gradually increase states kept
cutoff!(sweeps,1E-10) # desired truncation error

psi0 = randomMPS(sites, 2)

energy, psi = dmrg(H,psi0,sweeps)

function meas_Sz(psi, n)
    psi = orthogonalize(psi,n)
    sn = siteind(psi, n)
    Sz = scalar(dag(prime(psi[n],"Site"))*op("Sz",sn)*psi[n])
    return real(Sz)
end

Sz_list = [meas_Sz(psi, n) for n = 1:N]
plot(1:N, Sz_list)