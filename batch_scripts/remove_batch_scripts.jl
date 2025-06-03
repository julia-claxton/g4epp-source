using Statistics, LinearAlgebra
using Glob

include("/Users/luna/Research/G4EPP_2.0/Frontend_Functions.jl")
energies_to_remove, pitch_angles_to_remove = _get_beam_locations()

for i in eachindex(energies_to_remove)
    path = "$(@__DIR__)/$(energies_to_remove[i])keV_$(pitch_angles_to_remove[i])deg.sh"
    if isfile(path) == false; continue; end
    rm(path)
end