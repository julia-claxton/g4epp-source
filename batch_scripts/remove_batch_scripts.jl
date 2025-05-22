using Statistics, LinearAlgebra
using Glob


energies_to_remove = []
pitch_angles_to_remove = []
@assert length(energies_to_remove) == length(pitch_angles_to_remove)

[rm("$(@__DIR__)/$(energies_to_remove[i])keV_$(pitch_angles_to_remove[i])deg.sh") for i in eachindex(energies_to_remove)]