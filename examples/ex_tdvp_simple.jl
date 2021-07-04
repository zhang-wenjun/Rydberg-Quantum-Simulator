using ITensors
##
function simpletdvp(H::MPO, psi0::MPS, ttot::Float64, dt::Float64; kwargs...)
    # Get sweep parameters
    # maximun bond dimension kept during sweeps
    maxdim = get(kwargs, :maxdim, 200)   
    # cutoff for Schmidt coefficients (singular values)
    cutoff = get(kwargs, :cutoff, 1E-10) 

    N = length(psi0)
    Nt = Int(ceil(ttot/dt)) + 1

    # Store the contraction for projected local 
    # effective Hamiltonians
    lpos = 0
    rpos = N + 1
    TNs = Vector{ITensor}(undef, N)

    # Begin sweeps
    psi_t = Vector{MPS}(undef, Nt+1)
    # copy the initial state
    psi_t[1]  = copy(psi0)
    # and move the canonical center to the first site
    orthogonalize!(psi_t[1], 1)
    energy = 0.0
    for nt = 1:Nt
        sw_time = @elapsed begin
        maxtruncerr = 0.0
        tcur = nt * dt
        psi = copy(psi_t[nt])
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

            # Two site forward evolution
            Heff2 = H[b] * H[b+1]
            b > 1 && (Heff2 *= TNs[b-1])
            b+1 < N && (Heff2 *= TNs[b+2])
            phi = psi[b] * psi[b+1]
            gate = exp(-0.5im * dt * Heff2)
            # The difference from gate*phi is:
            # the result of apply(gate, phi) has 
            # unprimed indices
            phi = apply(gate, phi)
            ortho = "left" 
            indices = inds(psi[b])
            A, B, spec = factorize(phi, indices;
                            tags = tags(linkind(psi, b)), 
                            ortho = ortho,
                            cutoff = cutoff,
                            maxdim = maxdim)
            psi[b] = A
            if b+1 != N
                # One-site backwards evolution
                # Update canonical center
                TNs[b] = A*H[b]*dag(prime(A))
                b > 1 && (TNs[b] *= TNs[b-1])
                Heff1 = TNs[b] * H[b+1]
                b+1 < N && (Heff1 *= TNs[b+2])
                lpos = b
                gate = exp(0.5im * dt * Heff1)
                psi[b+1] = apply(gate, B)
            end

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

            # Two site forward evolution
            Heff2 = H[b] * H[b+1]
            b > 1 && (Heff2 *= TNs[b-1])
            b+1 < N && (Heff2 *= TNs[b+2])
            phi = psi[b] * psi[b+1]
            gate = exp(-0.5im * dt * Heff2)
            # The difference from gate*phi is:
            # the result of apply(gate, phi) has 
            # unprimed indices
            phi = apply(gate, phi)
            ortho = "right" 
            indices = inds(psi[b])
            A, B, spec = factorize(phi, indices;
                            tags = tags(linkind(psi, b)), 
                            ortho = ortho,
                            cutoff = cutoff,
                            maxdim = maxdim)
            psi[b+1] = B

            if b != 1
                # One-site backwards evolution
                # Update canonical center
                TNs[b+1] = B*H[b+1]*dag(prime(B))
                b+1 < N && (TNs[b+1] *= TNs[b+2])
                Heff1 = TNs[b+1] * H[b]
                b > 1 && (Heff1 *= TNs[b-1])
                rpos = b + 1
                gate = exp(0.5im * dt * Heff1)
                psi[b] = apply(gate, A)
            end

            # Truncation error
            maxtruncerr = max(maxtruncerr, spec.truncerr)
        end
        psi_t[nt+1] = copy(psi)
        end
        maxlinkdimpsi = maxlinkdim(psi)
        println("After sweep $nt maxlinkdim=$maxlinkdimpsi maxerr=$maxtruncerr time=$sw_time")
    end
    return psi_t
end

##
function meas_Sz(psi, n)
    psi = orthogonalize(psi,n)
    sn = siteind(psi, n)
    Sz = scalar(dag(prime(psi[n],"Site"))*op("Sz",sn)*psi[n])
    return real(Sz)
end

function meas_Sx(psi, n)
    psi = orthogonalize(psi,n)
    sn = siteind(psi, n)
    Sz = scalar(dag(prime(psi[n],"Site"))*op("Sx",sn)*psi[n])
    return real(Sz)
end
##
include("exactdiag.jl")
## 
# Test
# Define a Hamiltonian
N = 8
sites = siteinds("S=1/2", N)
alpha = 1.5
ampo = AutoMPO()
for i = 1:N-1
    # ampo .+= "Sx", i, "Sx", i+1
    for j = i+1:N
        d = abs(i-j)
        ampo .+= π/d^alpha, "Sx", i, "Sx", j
        ampo .+= π/d^alpha, "Sy", i, "Sy", j
    end
    # ampo .+= 2π*1, "ProjUp", i, "ProjUp", i+1
    # ampo .+= 2π, "Sx", i
end
# ampo .+= "Sx", N
H = MPO(ampo, sites);

##
ttot = 1.0
dt = 0.02
cutoff = 1e-9
psi0 = randomMPS(sites, 2)
##
psi_t_tdvp = simpletdvp(H, psi0, ttot, dt; cutoff = cutoff);
##
println("====================================")
ITensors.disable_warn_order()
include("exactdiag.jl")
psi_t_ed = ed_evol(H, psi0, ttot, dt);

##
t_list = Vector(0:dt:ttot)
obs_tdvp = [meas_Sz(psi_t_tdvp[nt], n) for n = 1:N, nt = 1:length(t_list)]
obs_ed = [meas_Sz(psi_t_ed[nt], n) for n = 1:N, nt = 1:length(t_list)]
h1 = heatmap(1:N, t_list, obs_tdvp);
h2 = heatmap(1:N, t_list, obs_ed);
plot(h1, h2, layout=(1,2))
##
plot(t_list, obs_tdvp[6, :])
plot!(t_list, obs_ed[6, :])