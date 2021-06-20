using ITensors

##
# Construct MPO with built-in `AutoMPO`
# Each `AutoMPO` instance contains a vector of built-in 
# `MPOTerm` struct accessed by `AutoMPO.data`
# Each `MPOTerm` instance contains a coefficient `MPOTerm.coef`
# and a vector of `OpTerm` struct, which is a vector of `SiteOp`
# Each `SiteOp` instance contains a string `name` specifying 
# the operator, and a `site` specifying the site where the operator acts.
N = 100
ampo = AutoMPO()
for j = 1:N-1
  # This is equivalent to
  # `add!(ampo, (0.5, "S+", j, "S-", j+1) )`
  ampo .+= 0.5,"S+",j,"S-",j+1
  ampo .+= 0.5,"S-",j,"S+",j+1
  ampo .+= "Sz",j,"Sz",j+1
end
H = MPO(ampo, sites)

# The function for converting AutoMPO to MPO
# For details see autompo.jl
# 
# function MPO(ampo::AutoMPO,
#     sites::Vector{<:Index};
#     kwargs...)::MPO
#     length(data(ampo)) == 0 && error("AutoMPO has no terms")

#     sorteachterm!(ampo,sites)
#     sortmergeterms!(ampo)

#     if hasqns(sites[1])
#         return qn_svdMPO(ampo,sites;kwargs...)
#     end
#     return svdMPO(ampo,sites;kwargs...)
# end
