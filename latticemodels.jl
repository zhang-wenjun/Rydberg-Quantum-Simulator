# A better way to realize 2D lattice for DMRG 
# calculation.
#
# Created by: Wenjun Zhang
#             University of Chicago
#             zhangwenjun1998@gmail.com
#             June 7, 2021
#
# Last modified: 20:15, June 7, 2021
# * Add `links` to `LatticeModel2D`

import Plots: plot, plot!, annotate!

###########################
# LocalTerm               # 
###########################

"""
A LocalTerm is a struct which represents
a single local term in a geometrical lattice
defining a physical model such as a quantum 
Hamiltonian.

LocalTerm has the following data fields:
* s::Vector{Int} -- numbers of the sites involved
* type::String   -- optional description of term type
"""
struct LocalTerm
    s::Vector{Int}
    type::String
end

"""
    LocalTerm(s1::Int,s2::Int)

    LocalTerm(s1::Int,s2::Int,
                type::String="")

    LocalTerm(s::Vector{Int})

    LocalTerm(s::Vector{Int},
                type::String="")

Construct a LocalTerm struct by
specifying just the numbers of sites
involved, or additional details such
as an optional type string.
"""
function LocalTerm(s1::Int, s2::Int,
                   type::String="")
    return LocalTerm([s1, s2], type)
end

function LocalTerm(s::Vector{Int})
    return LocalTerm(s, "")
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
* Nx::Int -- total number of sites in ``x`` direction
* Ny::Int -- total number of sites in ``y`` direction
* xs::Vector{Float64} -- ``x``-coordinates of sites
* ys::Vector{Float64} -- ``y``-coordinates of sites
* locterms::Vector{LocalTerm} -- array of local terms
"""
struct LatticeModel2D
    Nx::Int
    Ny::Int
    xs::Vector{Float64}
    ys::Vector{Float64}
    locterms::Vector{LocalTerm}
    links::Vector{Link}
end

# TODO: add an initializer function
# TODO: add a plot function to show local 
#       observables by using colormap

##
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
##

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
          marker_z = loc_obs,
          markershape = :circle,
          markercolor = :match,
          markersize = 6,
          seriescolor = cmap,
          linestyle = :dash,
          linecolor = :grey,
          label = "",
          aspect_ratio = :equal)
    return p
end
##

################################
# Functions to create lattices # 
################################

"""
    square_lattice_wjz(Nx::Int,
                   Ny::Int;
                   kwargs...)::LatticeModel2D

Return a LatticeModel2D corresponding to the 
two-dimensional square lattice of dimensions 
(Nx,Ny). By default the lattice has open 
boundaries, but can be made periodic in the 
y direction by specifying the keyword argument 
`yperiodic=true`.
"""
function square_lattice_wjz(Nx::Int,
                        Ny::Int;
                        kwargs...)::LatticeModel2D
    yperiodic = get(kwargs, :yperiodic, false)
    yperiodic = yperiodic && (Ny > 2)
    a = get(kwargs, :lattconst, 1)

    N = Nx*Ny
    Nbond = 2N-Ny + (yperiodic ? 0 : -Nx)
    locterms = Vector{LocalTerm}(undef, Nbond)
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
            locterms[b+=1] = LocalTerm(n, n+Ny, "x")
            links[l+=1] = Link([xs[n],xs[n]+a], [ys[n],ys[n]])
        end
        if Ny > 1
            if y < Ny
                # Create a $y$-bond
                locterms[b+=1] = LocalTerm(n, n+1, "y");
                links[l+=1] = Link([xs[n],xs[n]], [ys[n],ys[n]+a])
            end
            if yperiodic && y==1
                # Create a $y$-bond due to the periodic 
                # boundary condition
                locterms[b+=1] = LocalTerm(n, n+Ny-1, "y")
            end
        end
    end
    return LatticeModel2D(Nx, Ny, xs, ys, locterms, links)
end

##
# TODO: add Wilson loops
"""
    toric_code(Nx::Int,
               Ny::Int;
               kwargs...)::LatticeModel2D

Return a LatticeModel2D corresponding to the 
toric-code (dual-square) lattice of dimensions 
(Nx,Ny), thus containing ``2 N_x N_y`` sites. 
By default the lattice has open boundaries, 
but can be made periodic in the y direction 
by specifying the keyword argument `yperiodic=true`.
"""
function toric_code(Nx::Int,
                    Ny::Int;
                    kwargs...)::LatticeModel2D
    yperiodic = get(kwargs, :yperiodic, false)
    yperiodic = yperiodic && (Ny > 2)
    a = get(kwargs, :lattconst, 1)

    N = 2*Nx*Ny
    # Nplaq = Nx*Ny
    # Nstar = Nx*Ny
    locterms = Vector{LocalTerm}(undef, N)
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

            locterms[n] = LocalTerm(s, "plaq")
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

            locterms[n] = LocalTerm(s, "star")
        end
    end
    return LatticeModel2D(Nx, Ny, xs, ys, locterms, links)
end


##
#TODO: modify the triangular_lattice to suit the LatticeModel2D
# """
#     triangular_Lattice(Nx::Int,
#                        Ny::Int;
#                        kwargs...)::Lattice

# Return a Lattice (array of LocalTerm
# objects) corresponding to the two-dimensional
# triangular lattice of dimensions (Nx,Ny).
# By default the lattice has open boundaries,
# but can be made periodic in the y direction
# by specifying the keyword argument 
# `yperiodic=true`.
# """
# function triangular_Lattice(Nx::Int,
#                             Ny::Int;
#                             kwargs...)::Lattice
#   yperiodic = get(kwargs,:yperiodic,false)
#   yperiodic = yperiodic && (Ny > 2)
#   N = Nx*Ny
#   Nbond = 3N-2Ny + (yperiodic ? 0 : -2Nx+1)
#   locterms = loctermsice(undef,Nbond)
#   b = 0
#   for n=1:N
#     x = div(n-1,Ny)+1
#     y = mod(n-1,Ny)+1

#     # x-direction bonds
#     if x < Nx
#       locterms[b+=1] = LocalTerm(n,n+Ny)
#     end

#     # 2d bonds
#     if Ny > 1
#       # vertical / y-periodic diagonal bond
#       if (n+1 <= N) && ((y < Ny) || yperiodic)
#         locterms[b+=1] = LocalTerm(n,n+1);
#       end
#       # periodic vertical bond
#       if yperiodic && y==1
#         locterms[b+=1] = LocalTerm(n,n+Ny-1)
#       end
#       # diagonal bonds
#       if x < Nx && y < Ny
#         locterms[b+=1] = LocalTerm(n,n+Ny+1)
#       end
#     end
#   end
#   return locterms
# end

