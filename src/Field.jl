mutable struct Field
    a::Matrix{Float64} #spectral coefficients
    type::Int8 #type for the angular expansion: 1 for even cosine, 2 for odd cosine, 3 for even sine, 4 for odd sine
    Field() = new(zeros(Nx,Ny),1)
end

function (F::Field)(x::Float64,y::Float64, Nx::Integer=Nx, Ny::Integer=Ny; dx::Integer=0,dy::Integer=0)
    s=0
    type = F.type
    a = F.a
    aux = type==3 ? 1 : 0

    for jj in 1:Nx
        for kk in 1:Ny
            s+=a[jj,kk] * dT(jj-1,x,dx) * dTrig(kk-1+aux,y,dy,type)
        end
    end
    return s
end

function (F::Field)(j::Integer,k::Integer,Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My, Nx::Integer=Nx, Ny::Integer=Ny; dx::Integer=0,dy::Integer=0)
    s::Float64=0.0
    for jj in 1:Nx
        for kk in 1:Ny
            s+=(F.a)[jj,kk] * Mx[dx+1,jj,j] * My[dy+1,kk,k,F.type]
        end
    end
    return s
end