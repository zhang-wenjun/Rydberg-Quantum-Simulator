using Plots: display
##
using ITensors
using Plots
using ColorSchemes

println("=======================================================")
println("Import packages: Complete!")

##
println("=======================================================")
println("Compiling: Complete!")
# All in MHz
N = 13
sites = siteinds("S=1/2", N)
Omega = 2π * 2
Delta_list = range(-2π * 4, stop=2π * 8, step=2π * 0.2)
# Rb_list = range(0.5, 4, step=0.5)
# Delta_list = [10*Omega]
Rb_list = [1.513, 2.432, 3.027]
# Rb_list = [3.6]


# Measurements
function meas_dens(psi, n)
    psi = orthogonalize(psi,n)
    sn = siteind(psi, n)
    Sz = scalar(dag(prime(psi[n],"Site"))*op("ProjUp",sn)*psi[n])
    return real(Sz)
end

# Calculation
dens_list = Array{Float64}(undef, length(Delta_list), length(Rb_list), N)
for id = 1:length(Delta_list)
    for ib = 1:length(Rb_list)
        # Simulation settings
        trunc = N
        if ib == 1
            sweeps = Sweeps(20) # number of sweeps is 20
            maxdim!(sweeps,200) # gradually increase states kept
            cutoff!(sweeps,1E-10) # desired truncation error
        elseif ib == 2
            sweeps = Sweeps(80) # number of sweeps is 20
            maxdim!(sweeps,200) # gradually increase states kept
            cutoff!(sweeps,1E-10) # desired truncation error
        else
            sweeps = Sweeps(200) # number of sweeps is 20
            maxdim!(sweeps,200) # gradually increase states kept
            cutoff!(sweeps,1E-10) # desired truncation error
        end
        
        Delta = Delta_list[id]
        Rb = Rb_list[ib]
        println("Begin calculation for Δ=", Delta, ", Rb=", Rb)
        ampo = AutoMPO()
        for j = 1:N
            # Single-qubit terms
            ampo .+=  Omega,     "Sx", j
            ampo .+= -Delta, "ProjUp", j
            # Interaction terms
            for i = j+1:min(j+trunc, N)
                Vij = Omega*(Rb/(i - j))^6
                ampo .+= Vij, "ProjUp", j, "ProjUp", i
            end
        end
        H = MPO(ampo,sites)

        psi0 = randomMPS(sites, 2)

        energy, psi = dmrg( H,psi0,sweeps)

        dens_list[id, ib, :] = [meas_dens(psi, n) for n = 1:N]
        # plt = plot(1:N, dens_list[id, ib])

        # display(plt)
    end
end

##
x_axis = collect(Delta_list/(2π))
y_axis = collect(1:N)

c_val = hcat(dens_list[:, 1, :])
f1 = heatmap(x_axis, y_axis, c_val');

c_val = hcat(dens_list[:, 2, :])
f2 = heatmap(x_axis, y_axis, c_val');

c_val = hcat(dens_list[:, 3, :])
f3 = heatmap(x_axis, y_axis, c_val');

plot(f1, f2, f3, layout = (3, 1))