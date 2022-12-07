function Ergosphere(f::Field, g::Field, h::Field, W::Field, rh::Float64=1.0, file_name::String = "ergosphere.dat", Npoints::Integer=100; q::Float64=0.0)
    sol=Array{Float64}(undef,Npoints+2)
    points=sort(append!(collect(pi/(2*Npoints) .* ((0:Npoints-1) .+ 1/2)),[0.0,pi/2]))

    V=Array{Float64}(undef,0,2)

    M=GetMass(f,g,h,W,rh)
    j=GetJ(f,g,h,W,rh)
    χ = j/M^2
    rhK=M/2*sqrt(1-χ^2-q^2)

    #=
    function xErgoKerr(y::Float64,M::Float64,χ::Float64,rh::Float64)
        return r2x(M/2 * (sqrt(1-χ^2*cos(y)^2) + χ*sin(y)),rh)
    end
    =#

    function xErgoKerrN(y::Float64,M::Float64,χ::Float64,rh::Float64,q::Float64)
        return rBLKerr2x(M*(1+sqrt(1-q^2-χ^2*cos(y)^2)),rh,χ,q)
    end

    println("Printing ergosphere to file...")

    for i in 1:Npoints+2
        y=points[i]
        
        if i==1
            sol[i] = fzero(x->gtt(x,y,f,g,h,W), -1.0)
        else
            sol[i] = fzero(x->gtt(x,y,f,g,h,W), xErgoKerrN(y,M,χ,rhK,q))
            #sol[i] = find_zero((x->gtt(x,y,f,g,h,W),x->gtt(x,y,f,g,h,W,dx=1)), rEKerr(y,M,χ,rhK), Roots.Newton())
        end

        r=2/(1-sol[i])
        x=r*sin(y)
        z=r*cos(y)

        V=vcat(V,[x z])

    end

    PrintData(file_name,V)
end