using Core: apply_type
mutable struct ExactDiag
    H::ITensor
    V::ITensor
    D::ITensor
    Vt::ITensor
    dim::Int
end

# TODO: Check if H has site indices
# TODO: Check if H is hermitian
ExactDiag(H::MPO; kwargs...) = ExactDiag(prod(H); kwargs)

function ExactDiag(H::ITensor; kwargs...)
    eigsolver = get(kwargs, :eigsolver, eigen)
    if eigsolver == eigen
        f = eigen(H; ishermitian=true)
        ie = commonind(f.D, f.V)
        return ExactDiag(H, f.V, f.D, f.Vt, dim(ie))
    else
        error("Currently $eigsolver is not supported!")
    end
end

function evol_opr(ed::ExactDiag, dt::Real)
    u = ITensor(ComplexF64, inds(ed.D))
    [u[n,n] = exp(-1im*dt*ed.D[n,n]) for n in 1:ed.dim]
    return dag(ed.V)*u*ed.Vt
end

evol_opr(H::MPO, tau::Float64; kwargs...) = evol_opr(ExactDiag(H; kwargs), tau)

function ed_evol(H::MPO, psi0::MPS, 
                 ttot::Float64, tau::Float64; 
                 kwargs...)
    # println("ttot=$ttot, dt=$tau")
    # println("Getting evolution operator")
    # exp(-1im*tau)
    # ed = ExactDiag(H)
    # u = ITensor(ComplexF64, siteinds(H, 1))
    # [u[n,n] = exp(-1im*tau) for n in 1:2]
    # # [u[n,n] = exp(-1im*tau*ed.D[n,n]) for n in 1:ed.dim]
    # # evo = dag(ed.V)*u*ed.Vt
    # # evo = evol_opr(ed, tau)
    # println("ttot=$ttot, tau=$tau")
    # println("Get evolution operator")
    evo = exp(-1im*tau*prod(H))
    Nt = Int(ceil(ttot/tau)) + 1
    psi_t = Vector{MPS}(undef, Nt+1)
    psi_t[1] = copy(psi0)
    for n = 1:Nt
        println("Evolve at step $n")
        psi_t[n+1] = apply(evo, psi_t[n])
    end
    return psi_t
end

# Solve local eigen problem
# f = eigen(H; ishermitian = true)
# f.D is the diagonalized matrix
# f.V is the right-eigenvectors
# f.Vt is the left-eigenvectors
# H*f.V≈f.D*f.Vt or dag(f.V)*f.D*f.Vt≈H
# i--H--i'
# ==
# i--dag(f.V)--ie ~ ie--f.D--ie' ~ ie'--f.Vt--i'
# where ~ means contraction