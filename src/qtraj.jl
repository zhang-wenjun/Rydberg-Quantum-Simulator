using Base: Float64
# Calculate the time evolution of open systems
# with MPS initial states by using quantum
# trajectories based on TDVP.
#
# See PRB 94, 165116 (2016) for a mathematical
# description of TDVP.
# See Quantum Sci. Technol. 4, 013001 (2019)
# for an introduction to numerical simulation
# with quantum trajectories.
#
# Created by: Wenjun Zhang
#             The University of Chicago
#             zhangwenjun1998@gmail.com
#             Jul 6, 2021
# Modified from `src/tdvp.jl`


import ITensors: check_hascommoninds, permute, makeR!, makeL!, @debug_check, @timeit_debug, @printf, orthocenter

#############################
# Wrappers of TDVP function #
#############################

"""
    tdvp(H::MPO,psi0::MPS,ttot::Real,tau::Real;kwargs...)
                    
Use the time-dependent variational principle (TDVP) algorithm
to evolve a matrix product state (MPS) under the Hamiltonian H,
represented as a matrix product operator (MPO). The MPS `psi0` 
is the initial state at ``t = 0``.
The Observer system is supported.

Returns:
* `psi_t::Vector{MPS}` - evolution of the state at ``t=n*tau``,
  `length(psi_t) = Int(ceil(ttot/tau))+1`, `psi_t[1] = psi0`
"""
function tdvp(H::MPO, psi0::MPS, 
              ttot::Real, tau::Real; 
              kwargs...)
  check_hascommoninds(siteinds, H, psi0)
  check_hascommoninds(siteinds, H, psi0')
  # Permute the indices to have a better memory layout
  # and minimize permutations
  H = permute(H, (linkind, siteinds, linkind))
  PH = ProjMPO(H)
  return tdvp(PH,psi0,ttot,tau;kwargs...)
end


"""
    tdvp(Hs::Vector{MPO},psi0::MPS,ttot::Real,tau::Real;kwargs...)
                    
Use the time-dependent variational principle (TDVP) algorithm
to evolve a matrix product state (MPS) under the Hamiltonian H,
represented as a matrix product operator (MPO). The MPS `psi0` 
is the initial state at ``t = 0``.
The `Observer` system is supported.

This version of `tdvp` accepts a representation of H as a
Vector of MPOs, Hs = [H1,H2,H3,...] such that H is defined
as H = H1+H2+H3+...
Note that this sum of MPOs is not actually computed; rather
the set of MPOs [H1,H2,H3,..] is efficiently looped over at 
each step of the TDVP algorithm when evolving the MPS.

Returns:
* `psi_t::Vector{MPS}` - evolution of the state at ``t=n*tau``,
  `length(psi_t) = Int(ceil(ttot/tau))+1`, `psi_t[1] = psi0`
"""
function tdvp(Hs::Vector{MPO}, psi0::MPS, 
              ttot::Real, tau::Real; 
              kwargs...)
  for H in Hs
    check_hascommoninds(siteinds, H, psi0)
    check_hascommoninds(siteinds, H, psi0')
  end
  Hs .= permute.(Hs, Ref((linkind, siteinds, linkind)))
  PHS = ProjMPOSum(Hs)
  return tdvp(PHS,psi0,ttot,tau;kwargs...)
end


###############################
# Calculate local Hamiltonian #
###############################

"""
    localHam(PH, psi0::MPS, s1::Int, s2::Int=s1)
                    
Calculate the local Hamiltonian of site 
`s1`, `s1+1`, ..., `s2`, given by Hamiltonian 
`PH::Union{ProjMPO, ProjMPOSum}` and MPS `psi0`
The default value of `s2` is `s1`.
"""
function localHam(PH::ProjMPO, s1::Int, s2::Int=s1)
  H = PH.H[s1]
  for s in s1+1:s2
    H *= PH.H[s]
  end
  s1 > 1 && (H *= PH.LR[s1-1])
  s2 < length(PH) && (H *= PH.LR[s2+1]) # previously incorrect as s2+1 < length(PH)
  return H
end

function localHam(PH::ProjMPOSum, s1::Int, s2::Int=s1)
  return sum([localHam(ph, s1, s2) for ph in PH.pm])
end


#############################
# The actual OSTDVP codes   #
#############################

function app_jump!(psi::MPS, 
                   c_ops::Array{ITensor, Nops}; 
                   kwargs...) where Nops
  cutoff = get(kwargs, :cutoff, 1e-8)
  sweep_is_done = get(kwargs, :sweep_is_done, false)

  if !sweep_is_done
    return
  end

  normfac = norm(psi)
  uw_probs = Float64[]
  c_ops_list = ITensor[]
  nops = 0
  if normfac < rand(Float64)
    for k = 1:N
      orthogonalize!(psi, k)
      sk = siteind(psi, k)
      for cop in c_ops[k]
        psi_pr = cop*psi[k]
        push!(uw_probs, scalar(dag(prime(psi_pr,"Site"))*psi_pr))
        push!(c_ops_list, cop)
      end
    end
    probs = cumsum(uw_probs ./ sum(uw_probs))
    rv = rand(Float64)
    for v = 1:length(probs)
      if probs[v] >= rv
        apply(c_ops_list[v], psi; cutoff=cutoff)
        break
      end
    end
  end
  psi *= 1/normfac
end

# TODO: specify ishermitian=false, normalize=false
function ostdvp(PH, psi0::MPS, ttot::Real, tau::Real,
                c_ops::Array{ITensor, Nops}; 
                kwargs...) where Nops
  if length(psi0) == 1
    error("`tdvp` currently does not support system sizes of 1. You can diagonalize the MPO tensor directly with tools like `LinearAlgebra.eigen`, `KrylovKit.eigsolve`, etc.")
  end

  @debug_check begin
    # Debug level checks
    # Enable with ITensors.enable_debug_checks()
    checkflux(psi0)
    checkflux(PH)
  end

  which_decomp::Union{String, Nothing} = get(kwargs, :which_decomp, nothing)
  svd_alg::String = get(kwargs, :svd_alg, "divide_and_conquer")
  obs = get(kwargs, :observer, NoObserver())
  outputlevel::Int = get(kwargs, :outputlevel, 1)

  maxDim = get(kwargs, :maxdim, 200)
  minDim = get(kwargs, :mindim, 1)
  cutOff = get(kwargs, :cutoff, 1e-8)

  N = length(psi0)
  Nt = Int(ceil(ttot/tau))
  psi_t = Vector{MPS}(undef, Nt+1)
  psi_t[1] = copy(psi0)

  if !isortho(psi_t[1]) || orthocenter(psi_t[1]) != 1
    orthogonalize!(psi_t[1],1)
  end
  @assert isortho(psi_t[1]) && orthocenter(psi_t[1]) == 1

  position!(PH, psi_t[1], 1)

  for sw=1:Nt
    sw_time = @elapsed begin
    maxtruncerr = 0.0
    psi = copy(psi_t[sw])
    for (b, ha) in sweepnext(N)
      @debug_check begin
        println("    Evolving sweep $sw, block ($b, $ha)")
        checkflux(psi)
        checkflux(PH)
      end
      # ha == 1 means sweeping towards the right
      # sow we need to let the left tensor be orthogonal
      ortho = ha == 1 ? "left" : "right"

      @timeit_debug timer "ostdvp: position!" begin
      PH.nsite = 2
      position!(PH, psi, b)
      end

      @debug_check begin
        checkflux(psi)
        checkflux(PH)
      end

      @timeit_debug timer "ostdvp: psi[b]*psi[b+1]" begin
      phi = psi[b] * psi[b+1]
      end

      @timeit_debug timer "ostdvp: two-site forward evolution" begin
      Heff = localHam(PH, b, b+1)
      gate = exp(-0.5im * tau * Heff)
      phi = apply(gate, phi)
      end

      @debug_check begin
        checkflux(phi)
      end

      @timeit_debug timer "ostdvp: replacebond!" begin
      spec = replacebond!(psi, b, phi; maxdim = maxDim,
                                       mindim = minDim,
                                       cutoff = cutOff,
                                       ortho = ortho,
                                       normalize = false,
                                       which_decomp = which_decomp,
                                       svd_alg = svd_alg)
      end
      maxtruncerr = max(maxtruncerr,spec.truncerr)

      @timeit_debug timer "ostdvp: one-site backwards evolution" begin
      if (ha == 1 && b+1 != N) || (ha == 2 && b != 1)
        b1 = ha == 1 ? b+1 : b
        phi1 = psi[b1]
        PH.nsite = 1
        position!(PH, psi, b1)
        Heff = localHam(PH, b1)
        gate = exp(0.5im * tau * Heff)
        phi1 = apply(gate, phi1)
        psi[b1] = phi1
      end
      end

      @debug_check begin
        checkflux(psi)
        checkflux(PH)
      end

      if outputlevel >= 2
        @printf("Sweep %d, half %d, bond (%d,%d)\n",sw,ha,b,b+1)
        @printf("  Truncated using cutoff=%.1E maxdim=%d mindim=%d\n",
                cutOff,maxDim,minDim)
        @printf("  Trunc. err=%.2E, bond dimension %d\n",spec.truncerr,dim(linkind(psi,b)))
      end

      sweep_is_done = (b==1 && ha==2)

      # Begin: stochastically apply jump operators
      @timeit_debug timer "ostdvp: stochastically apply jump operators" begin
      app_jump!(psi, c_ops; sweep_is_done = sweep_is_done, cutoff=cutoff)
      end
      # End: stochastically apply jump operators

      measure!(obs; psi = psi,
                    bond = b,
                    sweep = sw,
                    half_sweep = ha,
                    spec = spec,
                    outputlevel = outputlevel,
                    sweep_is_done = sweep_is_done)
    end
    psi_t[sw+1] = copy(psi)
    end
    if outputlevel >= 1
      @printf("After sweep %d (t = %.2f) maxlinkdim=%d maxerr=%.2E time=%.3f\n",
              sw, sw*tau, maxlinkdim(psi), maxtruncerr, sw_time)
    end
  end
  return psi_t
end

