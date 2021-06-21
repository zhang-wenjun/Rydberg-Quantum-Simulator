# `PH` 
# the projected Hamiltonian
# used to keep track of the tensor network contraction
# for local effective Hamiltonian

# `checkdone!`
# if it has reached minimum number of sweeps specified, 
# and the energy difference between two sweeps is lower 
# than the tolerance specified, then tells the program
# no more sweeps are needed by returning true

# `sweepnext(N)` 
# gives an iterable, equivalent to local sites 
# of (1,2), (2,3), ..., (N-1,N), (N-1, N), ..., (2,3), (1,2)
# That is, going to right and then back left is called a sweep
# ha == 1 means it's going to right, 
# which corresponds to ortho == "left", 
# meaning the updated local tensor should be in left-canonical form
# ha == 2 means it's going to left, 
# which corresponds to ortho == "right", 
# meaning the updated local tensor should be in left-canonical form
# b == 1 and ha == 2 together mean the sweep is over

psi = copy(psi0)
N = length(psi)
orthogonalize!(psi,1)
position!(PH, psi, 1)
energy = 0.0
for sw=1:nsweep(sweeps)
    maxtruncerr = 0.0
    for (b, ha) in sweepnext(N)
        # solve local eigen-problem
        position!(PH, psi, b)
        phi = psi[b] * psi[b+1]
        vals, vecs = eigsolve(PH, phi, 1, eigsolve_which_eigenvalue;
                            ishermitian = ishermitian,
                            tol = eigsolve_tol,
                            krylovdim = eigsolve_krylovdim,
                            maxiter = eigsolve_maxiter)
        energy, phi = vals[1], vecs[1]
        #?
        ortho = ha == 1 ? "left" : "right"
        # update local tensors and bonds
        spec = replacebond!(psi, b, phi; maxdim = maxdim(sweeps, sw),
                                       mindim = mindim(sweeps, sw),
                                       cutoff = cutoff(sweeps, sw),
                                       eigen_perturbation = drho,
                                       ortho = ortho,
                                       normalize = true,
                                       which_decomp = which_decomp,
                                       svd_alg = svd_alg)

        # truncation error
        maxtruncerr = max(maxtruncerr,spec.truncerr)

        sweep_is_done = (b==1 && ha==2)
        measure!(obs; energy=energy,
                        psi=psi,
                        bond = b,
                        sweep = sw,
                        half_sweep = ha,
                        spec=spec,
                        outputlevel=outputlevel,
                        sweep_is_done=sweep_is_done)
    end
    isdone = checkdone!(obs;energy=energy,
                            psi=psi,
                            sweep=sw,
                            outputlevel=outputlevel) 

    isdone && break
end
