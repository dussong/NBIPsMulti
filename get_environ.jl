using JuLIP, NeighbourLists

import JuLIP.site_energies

function site_energies(pp, at::AbstractAtoms,rcut)
   Es = zeros(length(at))
   Z = atomic_numbers(at)
   for (i, j, r, R) in sites(at, rcut)
      for k in 1:length(j)
         atnb = sort([Z[i],Z[j][k]])
         tup_at = (chemical_symbol(atnb[1]),chemical_symbol(atnb[2]))
         Es[i] += pp[tup_at](r[k]) 
      end
   end
   return Es
end




Vcucu(x) = x
Vcuzr(x) = 2*x
Vzrzr(x) = 3*x
V2s = Dict(
            (:Cu, :Cu) => Vcucu,
            (:Cu, :Zr) => Vcuzr,
            (:Zr, :Zr) => Vcuzr
         )

atu = bulk(:Cu, cubic=true)*2
rcut = 2.7
site_energies(V2s, atu,rcut)

site(atu, 1)

at = bulk(:Cu, cubic=true)
Es = zeros(length(at))
Z = atomic_numbers(at)

V2s.(:Cu, :Cu)

JuLIP.pairs(at, 2.7)

atu = bulk(:Cu, cubic=true)
positions(atu) |> mat
Z = atomic_numbers(atu)
Z[[2,4]] = atomic_number(:Zr)
atu.Z[:] = Z[:]

at = atu * 3
chemical_symbols(at)

rnn(:Cu)



chemical_symbol(40)



sort(:Cu,:Zr)

V2s[(:Cu, :Cu)]

@show(sites(atu,3.5).nlist)

ctr = 0
S = Val.(chemical_symbols(at))
for (i, j, r, R) in sites(at, 3.5)
   symi = S[i]
   symj = S[j]
   @show i
   @show symi, symj
   ctr += 1
   if ctr == 10
      break
   end
end
