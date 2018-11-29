

function evaluate(V::NBodyFunction{N},
                  desc::NBSiteDescriptor,
                  Rs::AbstractVector{JVec{T}},
                  J::SVector{K, Integer},
                  Spi::Integer,Spj::Vector{Integer},Species::Vector{Integer}) where {N, T, K}
   # check species
   skip_simplex_species(Spi,Spj,Species,J) && return zero(T)
   skip_simplex_species_order!(desc,Spi,Spj,Species,J) && return zero(T)
   return evaluate(V,desc,Rs,J)
   # get the physical descriptor: bond-lengths (+ bond-angles)
   # rθ = ricoords(desc, Rs, J)
   # # # order the variables
   # # rθ = sort_ricoords!(desc, Rs, J, Spi, Spj)
   # # check whether to skip this N-body term?
   # skip_simplex(desc, rθ) && return zero(T)
   # # compute the cutoff (and skip this site if the cutoff is zero)
   # fc = fcut(desc, rθ)
   # fc == 0 && return zero(T)
   # # compute the invariants (this also applies the transform)
   # II = invariants(desc, rθ)
   # # evaluate the inner potential function (e.g. polynomial)
   # return evaluate_I(V, II) * fc
end

evaluate(V::NBodyFunction,
         Rs::AbstractVector{JVec{T}},
         J::SVector{K, Integer},
         Spi,Spj,Species) where {T,K} = evaluate(V,descriptor(V),Rs,J,Spi,Spj,Species)

function evaluate_d!(dVsite,
                     V::NBodyFunction{N},
                     desc::NBSiteDescriptor,
                     Rs::AbstractVector{JVec{T}},
                     J,
                     Spi::Integer,Spj::Vector{Integer},Species::Vector{Integer}) where {N,T}
   # check species
   skip_simplex_species(Spi,Spj,Species,J) && return dVsite
   skip_simplex_species_order!(desc,Spi,Spj,Species,J) && return dVsite
   evaluate_d!(dVsite, V, desc, Rs, J)
end


function evaluate_many!(Es,
                        B::AbstractVector{TB},
                        desc::NBSiteDescriptor,
                        Rs, J, Spi,Spj,Species)  where {TB <: NBodyFunctionM{N}} where {N}
   ind = find(skip_simplex_species_many(Spi,Spj,Species,J))
   # Es[ind] = evaluate_many!(Es[ind],B[ind],desc,Rs,J)
   for i in ind
      if !(skip_simplex_species_order!(desc,Spi,Spj,Species[i],J))
         @show B[i]
         @show typeof(Rs)
         @show typeof(J)
         Es[i] += evaluate(B[i],desc,Rs,J)
      end
   end
   return Es
end

evaluate_many!(out, B, Rs, J, Spi, Spj, Species) =
      evaluate_many!(out, B, descriptor(B[1]), Rs, J, Spi, Spj, Species)



function evaluate_many_d!(dVsite::AbstractVector,
                          B::AbstractVector{TB},
                          desc::NBSiteDescriptor,
                          Rs,
                          J, Spi,Spj,Species)  where {TB <: NBodyFunctionM{N}} where {N}
   ind = find(skip_simplex_species_many(Spi,Spj,Species,J))
   # dVsite[ind] = evaluate_many_d!(dVsite[ind],B[ind],desc,Rs,J)
   for i in ind
      if !(skip_simplex_species_order!(desc,Spi,Spj,Species[i],J))
         dVsite[i] = evaluate_d!(dVsite[i], B[i], desc, Rs, J)
      end
   end
   return dVsite
end

evaluate_many_d!(out, B, Rs, J, Spi, Spj, Species) =
      evaluate_many_d!(out, B, descriptor(B[1]), Rs, J, Spi,Spj,Species)


include("descriptors_multi_3B.jl")
