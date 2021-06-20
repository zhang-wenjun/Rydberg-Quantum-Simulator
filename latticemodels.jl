# A better way to realize 2D lattice for DMRG 
# calculation.
#
# Created by: Wenjun Zhang
#             University of Chicago
#             zhangwenjun1998@gmail.com
#             June 7, 2021
#
# Modified on: 20:15 (GMT+8), June 7, 2021
# * Add `links` to `LatticeModel2D`
#
# Modified on: 21:21 (GMT+8), June 16, 2021
# * Change name of `LocalTerm` to `CplTerm`

import Plots: plot, plot!, annotate!

###########################
# CplTerm               # 
###########################

"""
A CplTerm is a struct which represents
a single local term in a geometrical lattice
defining a physical model such as a quantum 
Hamiltonian.

CplTerm has the following data fields:
* s::Vector{Int} -- numbers of the sites involved
* type::String   -- optional description of term type
"""
struct CplTerm
    s::Vector{Int}
    type::String
end

"""
    CplTerm(s1::Int,s2::Int)

    CplTerm(s1::Int,s2::Int,
                type::String="")

    CplTerm(s::Vector{Int})

    CplTerm(s::Vector{Int},
                type::String="")

Construct a CplTerm struct by
specifying just the numbers of sites
involved, or additional details such
as an optional type string.
"""
function CplTerm(s1::Int, s2::Int,
                 type::String="")
    return CplTerm([s1, s2], type)
end

function CplTerm(s::Vector{Int})
    return CplTerm(s, "")
end



###########################
# Link                    # 
###########################

"""
A Link is a struct which represents
a link to be plotted when showing a lattice. 

Link has the following data fields:
* xs::Vector{Float64} -- x coordinates of the two ends of the link
* ys::Vector{Float64} -- y coordinates of the two ends of the link
"""
struct Link
    xs::Vector{Float64}
    ys::Vector{Float64}
end



###########################
# LatticeModel2D          # 
###########################

"""
A LatticeModel2D is a struct which represents
a lattice structure.

LatticeModel2D has the following data fields:
* xs::Vector{Float64} -- ``x``-coordinates of sites
* ys::Vector{Float64} -- ``y``-coordinates of sites
* cplterms::Vector{CplTerm} -- array of local terms
* (Deprecated) Nx::Int -- total number of sites 
                            in ``x`` direction
* (Deprecated) Ny::Int -- total number of sites 
                            in ``y`` direction
"""
struct LatticeModel2D
    xs::Vector{Float64}
    ys::Vector{Float64}
    cplterms::Vector{CplTerm}
    links::Vector{Link}
end

# TODO: add an initializer function

#
"""
    plot(latt::LatticeModel2D)

Use plot to show a LatticeModel
and its dimension-lowering path.
"""
function plot(latt::LatticeModel2D)
    p = plot()
    # Plot links
    for n = 1:length(latt.links)
        plot!(p, latt.links[n].xs, latt.links[n].ys,
              linestyle = :solid,
              linecolor = :black,
              label = "")
    end
    # Add label to sites
    a = latt.ys[2] - latt.ys[1]
    N = length(latt.xs)
    maxpad = length(string(N))
    annotate!(p, latt.xs .+ maxpad*0.1*a, latt.ys .+ 0.12*a, 
              string.(Vector(1:N), pad=maxpad), 
              color = :blue)
    
    plot!(p, latt.xs, latt.ys,
          markershape = :circle,
          markercolor = :black,
          markersize = 6,
          linestyle = :dash,
          linecolor = :red,
          label = "",
          aspect_ratio = :equal)
    return p
end
#

"""
    plot(latt::LatticeModel2D, 
         loc_obs::Vector{Float64})

Use plot to show a LatticeModel
and its dimension-lowering path, 
together with local observables shown in color.
"""
function plot(latt::LatticeModel2D, 
              loc_obs::Vector{Float64};
              kwargs...)
    cmap = get(kwargs, :cmap, :bwr)
    plot_label = get(kwargs, :plot_label, false)
    plot_path = get(kwargs, :plot_path, false)
    p = plot()
    # Plot links
    for n = 1:length(latt.links)
        plot!(p, latt.links[n].xs, latt.links[n].ys,
              linestyle = :solid,
              linecolor = :black,
              label = "")
    end
    if plot_label
        # Add label to sites
        a = latt.ys[2] - latt.ys[1]
        N = length(latt.xs)
        maxpad = length(string(N))
        annotate!(p, latt.xs .+ maxpad*0.1*a, latt.ys .+ 0.12*a, 
                string.(Vector(1:N), pad=maxpad), 
                color = :blue)
    end
    
    lt = plot_path ? :line : :scatter
    plot!(p, latt.xs, latt.ys,
          marker_z = loc_obs,
          markershape = :circle,
          markercolor = :match,
          markersize = 6,
          seriescolor = cmap,
          linestyle = :dash,
          linetype = lt,
          linecolor = :grey,
          label = "",
          aspect_ratio = :equal)
    return p
end

################################
# Functions to create lattices #
################################

"""
    square(Nx::Int, Ny::Int;
           kwargs...)::LatticeModel2D

Return a LatticeModel2D corresponding to the 
two-dimensional square lattice of dimensions 
(Nx,Ny). By default the lattice has open 
boundaries, but can be made periodic in the 
y direction by specifying the keyword argument 
`yperiodic=true`.
"""
function square(Nx::Int, Ny::Int;
                kwargs...)::LatticeModel2D
    yperiodic = get(kwargs, :yperiodic, false)
    yperiodic = yperiodic && (Ny > 2)
    a = get(kwargs, :lattconst, 1)

    N = Nx*Ny
    Nbond = 2N-Ny + (yperiodic ? 0 : -Nx)
    cplterms = Vector{CplTerm}(undef, Nbond)
    links = Vector{Link}(undef, 2N-Nx-Ny)
    xs = Vector{Float64}(undef, N)
    ys = Vector{Float64}(undef, N)
    b = 0
    l = 0
    for n = 1:N
        xs[n] = ((x = div(n-1, Ny) + 1) - 1) * a
        ys[n] = ((y = mod(n-1, Ny) + 1) - 1) * a
        if x < Nx
            # Create a $x$-bond
            cplterms[b+=1] = CplTerm(n, n+Ny, "x")
            links[l+=1] = Link([xs[n],xs[n]+a], [ys[n],ys[n]])
        end
        if Ny > 1
            if y < Ny
                # Create a $y$-bond
                cplterms[b+=1] = CplTerm(n, n+1, "y");
                links[l+=1] = Link([xs[n],xs[n]], [ys[n],ys[n]+a])
            end
            if yperiodic && y==1
                # Create a $y$-bond due to the periodic 
                # boundary condition
                cplterms[b+=1] = CplTerm(n, n+Ny-1, "y")
            end
        end
    end
    return LatticeModel2D(xs, ys, cplterms, links)
end

#
# TODO: add Wilson loops
"""
    toric_code(Nx::Int, Ny::Int;
               kwargs...)::LatticeModel2D

Return a LatticeModel2D corresponding to the 
toric-code (dual-square) lattice of dimensions 
(Nx,Ny), thus containing ``2 N_x N_y`` sites. 
By default the lattice has open boundaries, 
but can be made periodic in the y direction 
by specifying the keyword argument `yperiodic=true`.
"""
function toric_code(Nx::Int, Ny::Int;
                    kwargs...)::LatticeModel2D
    yperiodic = get(kwargs, :yperiodic, false)
    yperiodic = yperiodic && (Ny > 2)
    a = get(kwargs, :lattconst, 1)

    N = 2*Nx*Ny
    # Nplaq = Nx*Ny
    # Nstar = Nx*Ny
    cplterms = Vector{CplTerm}(undef, N)
    xs = Vector{Float64}(undef, N)
    ys = Vector{Float64}(undef, N)
    links = Vector{Link}(undef, N)
    for n = 1:N
        x = div(n-1, Ny) + 1
        y = mod(n-1, Ny) + 1

        xs[n] = 0.5*a * (x-1)
        ys[n] = a * (y-1) + (x % 2 == 0 ? 0 : 0.5*a)

        if x % 2 == 1
            links[n] = Link([xs[n], xs[n]], 
                            [ys[n]-0.5*a, ys[n]+0.5*a])
            # Create a plaquette term
            s = Int[]
            push!(s, n)
            push!(s, n+Ny)
            if x < Nx
                push!(s, n+2Ny)
            end
            if y < Ny
                push!(s, n+Ny+1)
            elseif yperiodic
                push!(s, n+1)
            end

            cplterms[n] = CplTerm(s, "plaq")
        else
            links[n] = Link([xs[n]-0.5*a, xs[n]+0.5*a], 
                            [ys[n], ys[n]])
            # Create a star term
            s = Int[]
            push!(s, n)
            push!(s, n-Ny)
            if x != 2
                push!(s, n-2Ny)
            end
            if y != 1
                push!(s, n-Ny-1)
            elseif yperiodic
                push!(s, n-1)
            end

            cplterms[n] = CplTerm(s, "star")
        end
    end
    return LatticeModel2D(xs, ys, cplterms, links)
end




#####################################
# Lattices with long-range coupling #
#####################################

"""
    lattice_lr(xs::Vector{Float64},
               ys::Vector{Float64};
               kwargs...)::LatticeModel2D

Return a LatticeModel2D corresponding to the 
given (ordered) list of coordinates of sites.
Two-body interaction terms are created and 
sorted according to their interaction ranges.


This function is useful when creating 
arbitrary lattice without a specific geometry.
Since the geometry is arbitrary, making it 
`yperiodic` is not supported.
"""
function lattice_lr(xs::Vector{Float64},
                    ys::Vector{Float64};
                    kwargs...)::LatticeModel2D
    N = length(xs) == length(ys) ? length(xs) :
          error("The lengths of `xs` and `ys` are different!")

    links = get(kwargs, :links, Link[])

    cplterms = CplTerm[]
    for n1 = 1:N
        for n2 = n1+1:N
            cp_range = sqrt((xs[n1] - xs[n2])^2 + (ys[n1] - ys[n2])^2)
            push!(cplterms, CplTerm([n1, n2], "cpr=$cp_range"))
        end
    end
    get_cpr(lct) = parse(Float64, lct.type[5:end])
    sort!(cplterms, by=get_cpr)
    return LatticeModel2D(xs, ys, cplterms, links)
end


#
"""
    square_lr(Nx::Int, Ny::Int;
              kwargs...)::LatticeModel2D

Return a LatticeModel2D corresponding to the 
two-dimensional square lattice of dimensions 
(Nx,Ny). By default the lattice has open 
boundaries, but can be made periodic in the 
y direction by specifying the keyword argument 
`yperiodic=true`.

Long-range couplings are included, with maximum
coupling range specified by the keyword argument
`max_nnorder`.
"""
function square_lr(Nx::Int, Ny::Int;
                  kwargs...)::LatticeModel2D
    # Get arguments
    yperiodic = get(kwargs, :yperiodic, false)
    yperiodic = yperiodic && (Ny > 2)
    lx = get(kwargs, :lx, 1)
    ly = get(kwargs, :ly, 1)
    max_nnorder = get(kwargs, :max_nnorder, 1)
    max_nnorder = max_nnorder < 1 ? 1 : max_nnorder 
    # Calculate the coordinates of sites
    # and links for drawing lattices
    N = Nx*Ny
    cplterms = CplTerm[]
    links = Vector{Link}(undef, 2N-Nx-Ny)
    xs = Vector{Float64}(undef, N)
    ys = Vector{Float64}(undef, N)
    l = 0
    for n = 1:N
        xs[n] = ((x = div(n-1, Ny) + 1) - 1) * lx
        ys[n] = ((y = mod(n-1, Ny) + 1) - 1) * ly
        if x < Nx
            links[l+=1] = Link([xs[n],xs[n]+lx], [ys[n],ys[n]])
        end
        if y < Ny
            links[l+=1] = Link([xs[n],xs[n]], [ys[n],ys[n]+ly])
        end
    end
    # Calculate coupling terms
    cplterms = CplTerm[]
    for n1 = 1:N
        for n2 = n1+1:N
            dx = abs(xs[n1] - xs[n2])
            dy = abs(ys[n1] - ys[n2])
            dy = yperiodic ? min(dy, Ny-dy) : dy
            cp_range = sqrt(dx^2 + dy^2)
            push!(cplterms, CplTerm([n1, n2], "cpr=$cp_range"))
        end
    end
    # Sort the coupling terms according to coupling range
    get_cpr(lct) = parse(Float64, lct.type[5:end])
    sort!(cplterms, by=get_cpr)
    # Extract the first `max_nnorder` nearest-neighbor 
    # coupling terms
    cpr_list = Float64[]
    n_cplt = 0
    for cplt in cplterms
        cpr = get_cpr(cplt)
        if length(cpr_list) == max_nnorder && cpr != cpr_list[end]
            break
        end
        n_cplt += 1
        if length(cpr_list) == 0 || cpr != cpr_list[end]
            push!(cpr_list, cpr)
        end
    end
    return LatticeModel2D(xs, ys, cplterms[1:n_cplt], links)
end
