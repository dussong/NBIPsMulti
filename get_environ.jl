using JuLIP, NeighbourLists

import JuLIP: site_energies, energy, partial_energy


struct PairPotentialM
   pot::Dict
end
# PairPotentialM : Dict of PairPotentials

# TODO: site_energies, energy, forces, partial_energy, partial_forces

function site_energies(ppM::PairPotentialM, at::AbstractAtoms,rcut)
   pp = ppM.pot
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

energy(ppM::PairPotentialM,at::AbstractAtoms,rcut) =
                                             sum(site_energies(ppM, at,rcut))

function partial_energy(pot_func,species,at::AbstractAtoms,rcut)
   E = 0.
   Z = atomic_numbers(at)
   for (i, j, r, R) in sites(at, rcut)
      for k in 1:length(j)
         atnb = sort([Z[i],Z[j][k]])
         tup_at = (chemical_symbol(atnb[1]),chemical_symbol(atnb[2]))
         if tup_at == species
            E += pot_func(r[k])
         end
      end
   end
   return E
end

function partial_forces(V,species, at::AbstractAtoms,rcut)
   dE = zerovecs(length(at))
   for (i, j, r, R) in sites(at, rcut)
      dV = @D V(r, R)
      for k in 1:length(j)
         atnb = sort([Z[i],Z[j][k]])
         tup_at = (chemical_symbol(atnb[1]),chemical_symbol(atnb[2]))
         if tup_at == species
            dE[j[k]] += dV[k]
            dE[i] += dV[k]
         end
      end
   end
   return dE
end

function forces(V::SitePotential, at::AbstractAtoms)
   frc = zerovecs(length(at))
   for (i, j, r, R) in sites(at, cutoff(V))
      dV = @D V(r, R)
      for a = 1:length(j)
         frc[j[a]] -= dV[a]
      end
      frc[i] += sum(dV)
   end
   return frc
end

Vcucu(x) = x
Vcuzr(x) = 2*x
Vzrzr(x) = 3*x
V2s = PairPotentialM(Dict(
            (:Cu, :Cu) => Vcucu,
            (:Cu, :Zr) => Vcuzr,
            (:Zr, :Zr) => Vcuzr
         ))
atu = bulk(:Cu, cubic=true)*2
rcut = 2.7
site_energies(V2s, atu,rcut)
energy(V2s,atu,rcut)


partial_energy(Vcucu,(:Cu,:Cu),atu,rcut)
partial_forces(Vcucu,(:Cu,:Cu),atu,rcut)
# site(atu,1.)

at = bulk(:Cu, cubic=true)
Es = zeros(length(at))
Z = atomic_numbers(at)

V2s[(:Cu, :Cu)]

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
