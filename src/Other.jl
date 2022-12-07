#Prints data array func to file file_name
function PrintData(file_name::String,func)
    open(file_name, "w") do io
        writedlm(io, func)
    end
end

function x2r(x::Float64,rh::Float64)
    return 2*rh/(1-x)
end

function r2x(r::Float64,rh::Float64)
    return 1-2*rh/r
end

function rBLKerr2x(rbl::Float64,rh::Float64,χ::Float64=0.0,q::Float64=0.0)
    M=2*rh/sqrt(1-χ^2-q^2)
    riso = 1/2*(-M+rbl+sqrt((M-rbl)^2-4*rh^2))
    return r2x(riso,rh)
end