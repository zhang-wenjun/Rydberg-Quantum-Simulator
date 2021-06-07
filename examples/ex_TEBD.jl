using ITensors
# The example from ITensors website 
# http://itensor.org/docs.cgi?vers=julia&page=getting_started/mps_time_evolution

##
let
    N = 100
    cutoff = 1E-8
    tau = 0.1
    ttotal = 5.0

    # Compute the number of steps to do
    Nsteps = Int(ttotal/tau)

    # Make an array of 'site' indices
    s = siteinds("S=1/2", N)

    # Make gates (1,2),(2,3),(3,4),...
    gates = ITensor[]
    for j = 1:N-1
    s1 = s[j]
    s2 = s[j+1]
    hj = op("Sz",s1) * op("Sz",s2) +
         1/2 * op("S+",s1) * op("S-",s2) +
         1/2 * op("S-",s1) * op("S+",s2)
    Gj = exp(-1.0im * tau/2 * hj)
    push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1), (N-1,N-2),...
    append!(gates, reverse(gates))

    # Function that measures <Sz> on site n
    function measure_Sz(psi, n)
        psi = orthogonalize(psi,n)
        sn = siteind(psi, n)
        Sz = scalar(dag(prime(psi[n],"Site"))*op("Sz",sn)*psi[n])
        return real(Sz)
    end

    # Initialize psi to be a product state (alternating up and down)
    psi = productMPS(s, n -> isodd(n) ? "Up" : "Dn")

    # Center position
    c = div(N, 2)

    # Compute and print initial <Sz> value
    t = 0.0
    Sz = measure_Sz(psi, c)
    println("$t $Sz")

    # Do the time evolution by applying the gates
    # for Nsteps steps
    for step = 1:Nsteps
        psi = apply(gates, psi; cutoff=cutoff)
        t += tau
        Sz = measure_Sz(psi, c)
        println("$t $Sz")
    end

    return
end