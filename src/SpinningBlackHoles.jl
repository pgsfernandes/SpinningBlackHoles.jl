module SpinningBlackHoles

using ClassicalOrthogonalPolynomials, DelimitedFiles, Distributions, Roots, Cubature

include("Field.jl")

include("Load.jl")

include("Other.jl")

include("Spectral.jl")

include("Quantities.jl")

include("KerrNewman.jl")

include("Petrov.jl")

include("MetricFunctions.jl")

include("Ergosphere.jl")

include("LightRing.jl")

include("ISCO.jl")

include("CurvatureScalars.jl")

LoadSystem()

export Field, LoadSystem, PrintData, interpolate1D, interpolate1D!, interpolate, interpolate!, GetMass, GetJ, GetTh, GetAh, GetωχKerr, GetωχKerrN, quantities_kerr, GetLe, GetLp, Sphericity, vH, get_quantities, fKerrN, gKerrN, hKerrN, WKerrN, AtKerrN, AφKerrN, Print_Petrov, gtt, gtphi, gphiphi, grr, gthetatheta, Ergosphere, CircumferencialRadius, LightRing, ISCO, Ricci_Horizon, GB_Horizon, RicciScalar, GBScalar, X, Mx, Y, My, Nx, Ny

end