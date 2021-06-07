using ITensors
# The example from ITensors website 
# http://itensor.org/docs.cgi?vers=julia&page=getting_started/mps_time_evolution

##
N = 10
cutoff = 1E-8

tau = 0.1
ttotal = 20.0

# Compute the number of steps to do
Nsteps = Int(ttotal/tau)

# Make an array of 'site' indices
s = siteinds("S=1/2", N)

# Make gates (1,2),(2,3),(3,4),...
gates = ITensor[]
for j = 1:N-1
    s1 = s[j]
    s2 = s[j+1]
    hj = op("Sz",s1) * op("Sz",s2)
    Gj = exp(-tau/2 * hj)
    push!(gates, Gj)
end
# hN = - op("Sz",s[N])
# GN = exp(-tau/2 * hN)
# push!(gates, GN)
# Include gates in reverse order too
# (N,N-1), (N-1,N-2),...
append!(gates, reverse(gates))

# Function that measures <Sz> on site n
function meas_Sz(psi, n)
    psi = orthogonalize(psi,n)
    sn = siteind(psi, n)
    Sz = scalar(dag(prime(psi[n],"Site"))*op("Sz",sn)*psi[n])
    return real(Sz)
end

# Initialize psi to be a product state (alternating up and down)
# psi = productMPS(s, n -> isodd(n) ? "Up" : "Dn")

# Initializing a random MPS
psi = randomMPS(s, 2)

t = 0.0

# Do the time evolution by applying the gates
# for Nsteps steps
for step = 1:Nsteps
    psi = apply(gates, psi; cutoff=cutoff)
    t += tau
end

psi *= 1/norm(psi)

Sz_list = [meas_Sz(psi, n) for n = 1:N]
plot(1:N, Sz_list)
