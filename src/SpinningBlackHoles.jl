module SpinningBlackHoles

using ClassicalOrthogonalPolynomials, DelimitedFiles, Distributions, Roots, Cubature

mutable struct Field
    a::Matrix{Float64} #spectral coefficients
    type::Int8 #type for the angular expansion: 1 for even cosine, 2 for odd cosine, 3 for even sine, 4 for odd sine
    Field() = new(zeros(Nx,Ny),1)
end

function (F::Field)(x::Float64,y::Float64, Nx::Integer=Nx, Ny::Integer=Ny; dx::Integer=0,dy::Integer=0)
    s=0
    for jj in 1:Nx
        for kk in 1:Ny
            s+=(F.a)[jj,kk] * dT(jj-1,x,dx) * dTrig(kk-1,y,dy,F.type)
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

function LoadSystem()
	
	if isdefined(Main,:Nx)
		@eval const Nx = Main.Nx
	else
		print("Enter x resolution: ")
		@eval const Nx = parse(Int, chomp(readline()))
		println()
	end

	if isdefined(Main,:Ny)
		@eval const Ny = Main.Ny
	else
		print("Enter y resolution: ")
		@eval const Ny = parse(Int, chomp(readline()))
		println()
	end

    @eval begin
        const X = xpoints(Nx)
        const Mx = GetChebyshev(Nx,X)
        const Y = ypoints(Ny)
        const My = GetTrig(Ny,Y)
    end
	println()
	println("Spinning Black Holes package loaded!")
	print("Resolution: Nx=",Nx,", Ny=",Ny)
	println()
	return X, Mx, Y, My
end

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

#APPROXIMATING A 1D FUNCTION f VIA CHEBYSHEV INTERPOLATION. RETURNS THE SPECTRAL COEFFICIENTS.
function interpolate1D(f::Function,Nx::Integer,Ny::Integer)
    c=zeros(Float64, (Nx,Ny))
    x = Array{Float64, 1}(undef, Nx)
    cheb=Array{Float64, 2}(undef, Nx, Nx)
    for i in 1:Nx
        x[i] = -cos((2*i-1.0)/(2*Nx)*pi)
    end
    for k in 0:Nx-1
        for i in 0:Nx-1
            cheb[k+1,i+1] = chebyshevt(k,x[i+1])
        end
    end
    for j in 0:(Nx-1)
        for x_it in 1:Nx
            c[j+1,1] += 2.0/Nx * f(x[x_it]) * cheb[j+1,x_it]
        end
    end
    c[1,1]/=2.0
    return c
end
function interpolate1D!(F::Field,f::Function,Nx::Integer,Ny::Integer)
    F.a=zeros(Float64, (Nx,Ny))
    x = Array{Float64, 1}(undef, Nx)
    cheb=Array{Float64, 2}(undef, Nx, Nx)
    for i in 1:Nx
        x[i] = -cos((2*i-1.0)/(2*Nx)*pi)
    end
    for k in 0:Nx-1
        for i in 0:Nx-1
            cheb[k+1,i+1] = chebyshevt(k,x[i+1])
        end
    end
    for j in 0:(Nx-1)
        for x_it in 1:Nx
            (F.a)[j+1,1] += 2.0/Nx * f(x[x_it]) * cheb[j+1,x_it]
        end
    end
    (F.a)[1,1]/=2.0
end

#APPROXIMATING A 2D FUNCTION f VIA CHEBYSHEV INTERPOLATION. RETURNS THE SPECTRAL COEFFICIENTS.
function interpolate(f::Function,Nx::Integer,Ny::Integer,type::Integer)
    c=zeros(Float64, (Nx,Ny))
    x = Array{Float64, 1}(undef, Nx)
    y = Array{Float64, 1}(undef, Ny)
    cheb=Array{Float64, 2}(undef, Nx, Nx)
    trig=Array{Float64, 2}(undef, Ny, Ny)
    
    for i in 1:Nx
        x[i] = -cos((2*i-1.0)/(2*Nx)*pi)
    end
    for i in 0:Ny-1
        y[i+1] = pi*(i+0.5)/(2*Ny)
    end

    for k in 0:Nx-1
        for i in 0:Nx-1
            cheb[k+1,i+1] = chebyshevt(k,x[i+1])
        end
    end

    for k in 0:Ny-1
        for i in 0:Ny-1
            trig[k+1,i+1] = dTrig(k,y[i+1],0,type)
        end
    end
    
    for j in 0:(Nx-1)
        for k in 0:(Ny-1)
            for x_it in 1:Nx
                for y_it in 1:Ny
                    if k==0 && (type == 2 || type == 4)
                        c[j+1,k+1] += 8.0/(Nx*Ny) * f(x[x_it],y[y_it]) * cheb[j+1,x_it] * trig[k+1,y_it]
                    else
                        c[j+1,k+1] += 4.0/(Nx*Ny) * f(x[x_it],y[y_it]) * cheb[j+1,x_it] * trig[k+1,y_it]
                    end
                end
            end
        end
    end
    
    for j in 1:Nx
        c[j,1] /= 2.0
    end
    for k in 1:Ny
        c[1,k] /= 2.0
    end
    
    return c
end

function interpolate!(F::Field,f::Function,Nx::Integer,Ny::Integer)
    F.a=zeros(Float64, (Nx,Ny))
    type=F.type
    x = Array{Float64, 1}(undef, Nx)
    y = Array{Float64, 1}(undef, Ny)
    cheb=Array{Float64, 2}(undef, Nx, Nx)
    trig=Array{Float64, 2}(undef, Ny, Ny)
    
    for i in 1:Nx
        x[i] = -cos((2*i-1.0)/(2*Nx)*pi)
    end
    for i in 0:Ny-1
        y[i+1] = pi*(i+0.5)/(2*Ny)
    end

    for k in 0:Nx-1
        for i in 0:Nx-1
            cheb[k+1,i+1] = chebyshevt(k,x[i+1])
        end
    end

    for k in 0:Ny-1
        for i in 0:Ny-1
            trig[k+1,i+1] = dTrig(k,y[i+1],0,type)
        end
    end
    
    for j in 0:(Nx-1)
        for k in 0:(Ny-1)
            for x_it in 1:Nx
                for y_it in 1:Ny
                    if k==0 && (type == 2 || type == 4)
                        (F.a)[j+1,k+1] += 8.0/(Nx*Ny) * f(x[x_it],y[y_it]) * cheb[j+1,x_it] * trig[k+1,y_it]
                    else
                        (F.a)[j+1,k+1] += 4.0/(Nx*Ny) * f(x[x_it],y[y_it]) * cheb[j+1,x_it] * trig[k+1,y_it]
                    end
                end
            end
        end
    end
    
    for j in 1:Nx
        (F.a)[j,1] /= 2.0
    end
    for k in 1:Ny
        (F.a)[1,k] /= 2.0
    end
end

#Converts coefficients of 1D approximation such that they can be used in a 2D problem. In our context, takes a static solution and puts it in a form such that it can be used as an approximation for a slowly rotating axisymmetric solution.
function To2D(v::AbstractArray,Ny::Integer)
    return hcat(v,zeros(length(v),Ny-1))
end

#Returns Chebyshev nodes as our grid in the x direction
function xpoints(Nx::Integer)
    #ALL POINTS. BOUNDARYS x=-1 AND x=1 ARE LOCATED AT POSITIONS 1 AND Nx RESPECTIVELY
    return sort(append!(cos.((2 .* (1:Nx-2) .- 1) ./ (2*(Nx-2))*pi),[-1.0,1.0]))
end

#Equally spaced points in the range [0,pi/2]. Our grid in y direction
function ypoints(Ny::Integer)
    #INTERIOR POINTS ARE LOCATED IN THE POSITIONS RANGING FROM 2 TO Ny+1, y=0 IS LOCATED AT POSITION 1, AND y=pi/2 IS LOCATED AT POSITION Ny+2
    return sort(append!(collect(pi/(2*Ny) .* ((0:Ny-1) .+ 1/2)),[0.0,pi/2]))
end

#COMPUTES THE CHEBYSHEV POLYNOMIAL OF ORDER n (OR ITS dx DERIVATIVE) AT POINT x
function dT(n::Integer,x::Float64,dx::Integer=0)
    if n<=0
        if dx==1 || dx==2
            return 0
        else
            return chebyshevt(-n,x)
        end
    else
        if dx==1
            return n*chebyshevu(n-1,x)
        elseif dx==2
            if x==-1
                return (n^4-n^2)/3 * (-1)^n
            elseif x==1
                return (n^4-n^2)/3
            else
                return n*(n*chebyshevt(n,x)-x*chebyshevu(n-1,x))/(-1+x^2)
            end
        else
            return chebyshevt(n,x)
        end
    end 
end

#COMPUTES THE EVEN COSINE OF ORDER n (OR ITS dy DERIVATIVE) AT POINT y
function dTrig(n::Integer,y::Float64,dy::Integer=0,type::Integer=1)
    if dy == 0
        if type==1
            return cos(2*n*y)
        elseif type==2
            return cos((2*n+1)*y)
        elseif type==3
            return sin(2*n*y)
        else
            return sin((2*n+1)*y)
        end
    elseif dy==1
        if type==1
            return -2*n*sin(2*n*y)
        elseif type==2
            return -(2*n+1)*sin((2*n+1)*y)
        elseif type==3
            return 2*n*cos(2*n*y)
        else
            return (2*n+1)*cos((2*n+1)*y)
        end
    else
        if type==1
            return -4*n*n*cos(2*n*y)
        elseif type==2
            return -(2*n+1)^2*cos((2*n+1)*y)
        elseif type==3
            return -4*n*n*sin(2*n*y)
        else
            return -(2*n+1)^2*sin((2*n+1)*y)
        end
    end
end

#Evaluates all Chebyshev polynomials, first and second derivatives in all points of our x grid. This allows us to not repeat computations when evaluating functions on the grid.
function GetChebyshev(Nx::Integer,x::Vector{Float64})
    C = zeros((3,Nx,Nx))
    for i in 1:3
        for j in 1:Nx
            for k in 1:Nx
                #derivative, poly order, position
                C[i,j,k] = dT(j-1,x[k],i-1)
            end
        end
    end
    return C
end

#Evaluates all even cosines, first and second derivatives in all points of our y grid. This allows us to not repeat computations when evaluating functions on the grid.
function GetTrig(Ny::Integer,y::Vector{Float64})
    C = zeros((3,Ny,Ny+2,4))
    for t in 1:4
        for i in 0:2
            for j in 0:Ny-1
                for k in 1:Ny+2
                    C[i+1,j+1,k,t] = dTrig(j,y[k],i,t)
                end
            end
        end
    end
    return C
end

#Computes ADM mass of the solution
function GetMass(f::Field,g::Field,h::Field,W::Field,rh::Float64)
    #return rh*(1+f(Nx,1,dx=1))
    return -rh*gtt(1.0,pi/2,f,g,h,W,rh,dx=1)
end

#Computes angular momentum of the solution. Assumes ds^2 contains (d\phi - W/r dt)^2
function GetJ(f::Field,g::Field,h::Field,W::Field,rh::Float64, Nx::Integer = Nx, Ny::Integer = Ny)
    #return -rh^2*W(Nx,Ny+2,dx=1)
    return rh*gtphi(1.0,pi/2,f,g,h,W,rh,dx=1)
end

#Computes event horizon temperature of the solution
function GetTh(f::Field,g::Field,h::Field,rh::Float64)
    #return 1/(2*pi*rh)*f(1,1)/sqrt(m(1,1))
    return 1/(2*pi*rh)*f(1,1)/sqrt(g(1,1)*h(1,1))
end

#Computes area of the event horizon of the solution
function GetAh(f::Field,g::Field,h::Field,rh::Float64)
    #Th is constant on the horizon
    #integral, err = quadgk(y -> sin(y)*g(-1.0,y)*sqrt(h(-1.0,y))/f(-1.0,y), 0, pi, rtol=1e-14)
    #integral, err = quadgk(y -> sin(y)*sqrt(g(-1.0,y)), 0, pi, rtol=1e-14)
    integral, err = pquadrature(y -> sin(y)*sqrt(g(-1.0,y)), 0, pi, abstol=1e-14)
    return rh/GetTh(f,g,h,rh) * integral
end

#THEORETICAL VALUES FOR M, J, T_H, A_H FOR THE KERR BLACK HOLE
function GetωχKerr(WBC::Float64,spin=true,branch::Integer=1)
    ω=0.0
    χ=0.0
    if spin
        χ=WBC
        if WBC!=0
            ω = (-1+WBC^2+sqrt(1-WBC^2))/(4*WBC)
        end
    else
        ω = WBC
        aux = sqrt(Complex(-3+528*ω^2+768*ω^4))
        if branch == 2
            χ=real((3+(4*ω+(-8*ω*(9+8*ω^2)+3*aux)^(1/3))^2)/(3*(-8*ω*(9+8*ω^2)+3*aux)^(1/3)))
        else
            χ=real(1/6*(16*ω+(im*(im+sqrt(3))*(3+16*ω^2))/(-8*ω*(9+8*ω^2)+3*aux)^(1/3)-(1+im*sqrt(3))*(-8*ω*(9+8*ω^2)+3*aux)^(1/3)))
        end
    end
    return [ω, χ]
end

function GetωχKerrN(WBC::Float64,q::Float64,spin=true,branch::Integer=1)
    ω=0.0
    χ=0.0
    if spin
        χ=WBC
        if WBC!=0
            ω = χ*sqrt(1-q^2-χ^2)/(4-2*q^2+4*sqrt(1-q^2-χ^2))
        end
    else
        ω = WBC
        aux = sqrt(Complex(-3+528*ω^2+768*ω^4))
        if branch == 2
            χ=real((3+(4*ω+(-8*ω*(9+8*ω^2)+3*aux)^(1/3))^2)/(3*(-8*ω*(9+8*ω^2)+3*aux)^(1/3)))
        else
            χ=real(1/6*(16*ω+(im*(im+sqrt(3))*(3+16*ω^2))/(-8*ω*(9+8*ω^2)+3*aux)^(1/3)-(1+im*sqrt(3))*(-8*ω*(9+8*ω^2)+3*aux)^(1/3)))
        end
        χ=find_zero(x->ω - x*sqrt(1-q^2-x^2)/(4-2*q^2+4*sqrt(1-q^2-x^2)), χ)
    end
    return [ω, χ]
end

function quantities_kerr(WBC::Float64,rh::Float64,spin=true,branch::Integer=1)
    ω, χ = GetωχKerr(WBC,spin,branch)

    M = 2*rh/sqrt(1-χ^2)
    J = χ*M^2
    Th = 1/(4*pi*M * (1+M/(2*rh)))
    Ah = 8*pi*M*(M+2*rh)
    return [M J χ Th Ah ω/rh]
end

function GetLe(f::Field,g::Field,h::Field,rh::Float64)
    #Returns length of the circumference along the equator
    return 2*pi*rh*sqrt(g(1,Ny+2)/f(1,Ny+2))
end
function GetLp(f::Field,g::Field,h::Field,rh::Float64)
    #Returns length of the circumference along the poles
    integral, err = pquadrature(y -> sqrt(g(-1.0,y)*h(-1.0,y)/f(-1.0,y)), 0, pi, abstol=1e-14)
    return 2*rh * integral
end
function Sphericity(f::Field,g::Field,h::Field,rh::Float64)
    return GetLe(f,g,h,rh)/GetLp(f,g,h,rh)
end
function vH(f::Field,g::Field,h::Field,rh::Float64,OmH::Float64)
    #Returns linear velocity of the horizon
    return GetLe(f,g,h,rh)/(2*pi) * OmH
end

#Returns all the above quantitites along with the dimensionless spin j/M^2
function get_quantities(f::Field,g::Field,h::Field,W::Field,rh::Float64)
    M=GetMass(f,g,h,W,rh)
    j=GetJ(f,g,h,W,rh)
    chi = j/M^2
    Th=GetTh(f,g,h,rh)
    Ah=GetAh(f,g,h,rh)
    Ωh=W(1,1)/rh
    return[M j chi Th Ah Ωh]
end

#KERR-NEWMAN SOLUTION
function auxKN1(x::Float64,y::Float64,rh::Float64,M::Float64,Q::Float64)
    return -((-2*M^2*(-1+x)^2-1/4*rh^2*(-3+x)^2*(1+x)^2+M*rh*(-1+x)*(5+(-2+x)*x)+Q^2*(-1+x)^2*cos(y)^2+(M^2-4*rh^2)*(-1+x)^2*sin(y)^2)/(4*rh^2))
end
function auxKN2(x::Float64,y::Float64,rh::Float64,M::Float64,Q::Float64)
    return ((8*M^2*(-1+x)^2-4*Q^2*(-1+x)^2+rh^2*(-3+x)^2*(1+x)^2-4*M*rh*(-1+x)*(5+(-2+x)*x))^2+4*rh^2*(-M^2+Q^2+4*rh^2)*(-3+x)^2*(-1+x)^2*(1+x)^2*sin(y)^2)/(256*rh^4)
end

function fKerrN(x::Float64,y::Float64,rh::Float64,ω::Float64,χ::Float64,q::Float64=0.0,spin::Bool=true,branch::Integer=1)
    M = 2*rh/sqrt(1-χ^2-q^2)
    Q=q*M

    return 1/4*(-3+x)^2*auxKN1(x,y,rh,M,Q)/auxKN2(x,y,rh,M,Q)
end
function gKerrN(x::Float64,y::Float64,rh::Float64,ω::Float64,χ::Float64,q::Float64=0.0,spin::Bool=true,branch::Integer=1)
    return 1/4*(-3+x)^2
end
function hKerrN(x::Float64,y::Float64,rh::Float64,ω::Float64,χ::Float64,q::Float64=0.0,spin::Bool=true,branch::Integer=1)
    M = 2*rh/sqrt(1-χ^2-q^2)
    Q=q*M

    return auxKN1(x,y,rh,M,Q)^2/auxKN2(x,y,rh,M,Q)
end
function WKerrN(x::Float64,y::Float64,rh::Float64,ω::Float64,χ::Float64,q::Float64=0.0,spin::Bool=true,branch::Integer=1)
    M = 2*rh/sqrt(1-χ^2-q^2)
    Q=q*M

    return -((sqrt(M^2-Q^2-4*rh^2)*(-1+x)*(-2*M^2*(-1+x)+Q^2*(-1+x)+M*rh*(5+(-2+x)*x)))/(4*rh^2))/auxKN2(x,y,rh,M,Q)/rh
end
function AtKerrN(x::Float64,y::Float64,rh::Float64,ω::Float64,χ::Float64,q::Float64=0.0,spin::Bool=true,branch::Integer=1)
    M = 2*rh/sqrt(1-χ^2-q^2)
    Q=q*M

    #return -((Q*(M-(rh*(5+(-2+x)*x))/(2*(-1+x))))/((M-(rh*(5+(-2+x)*x))/(2*(-1+x)))^2+(M^2-Q^2-4*rh^2)*cos(y)^2))
    return -((Q*(M+2*rh))/(Q^2-2*M*(M+2*rh)))-(Q*(M-(rh*(5+(-2+x)*x))/(2*(-1+x)))*(1-(16*(M^2-Q^2-4*rh^2)*(-1+x)^3*(2*M^2*(-1+x)-Q^2*(-1+x)-M*rh*(5+(-2+x)*x))*sin(y)^2)/((8*M^2*(-1+x)^2-4*Q^2*(-1+x)^2+rh^2*(-3+x)^2*(1+x)^2-4*M*rh*(-1+x)*(5+(-2+x)*x))^2+4*rh^2*(-M^2+Q^2+4*rh^2)*(-3+x)^2*(-1+x)^2*(1+x)^2*sin(y)^2)))/((M-(rh*(5+(-2+x)*x))/(2*(-1+x)))^2+(M^2-Q^2-4*rh^2)*cos(y)^2)
end
function AφKerrN(x::Float64,y::Float64,rh::Float64,ω::Float64,χ::Float64,q::Float64=0.0,spin::Bool=true,branch::Integer=1)
    M = 2*rh/sqrt(1-χ^2-q^2)
    Q=q*M

    #return (Q*sqrt(M^2-Q^2-4*rh^2)*(M-(rh*(5+(-2+x)*x))/(2*(-1+x)))*sin(y)^2)/((M-(rh*(5+(-2+x)*x))/(2*(-1+x)))^2+(M^2-Q^2-4*rh^2)*cos(y)^2)
    return (Q*sqrt(M^2-Q^2-4*rh^2)*(M-(rh*(5+(-2+x)*x))/(2*(-1+x))))/((M-(rh*(5+(-2+x)*x))/(2*(-1+x)))^2+(M^2-Q^2-4*rh^2)*cos(y)^2)
end

# PETROV TYPE
function psi_0(i::Integer, j::Integer, f::Field,g::Field,h::Field,W::Field,rh::Float64,x::Vector{Float64}=X,y::Vector{Float64}=Y, Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    return ((-1 + x[i])^2*csc(y[j])*(3*(1 + x[i])^2*f(i,j)^2*h(i,j)*sin(y[j])*(im*g(i,j,dy=1) + (-1 + x[i])*g(i,j,dx=1))^2 + 2*(-1 + x[i])*(1 + x[i])*f(i,j)*g(i,j)^(3/2)*h(i,j)*sin(y[j])^2*(im*g(i,j,dy=1) + (-1 + x[i])*g(i,j,dx=1))*(2*W(i,j) + im*W(i,j,dy=1) + (-1 + x[i])*W(i,j,dx=1)) - 4*(-1 + x[i])^2*g(i,j)^3*h(i,j)*sin(y[j])^3*(2*W(i,j) + im*W(i,j,dy=1) + (-1 + x[i])*W(i,j,dx=1))^2 - 4*(1 + x[i])*f(i,j)*g(i,j)^2*(f(i,j)*(-3*(-1 + x[i])*h(i,j)*sin(y[j]) + ((1 + x[i])*cos(y[j]) + (2*im)*x[i]*sin(y[j]))*(h(i,j,dy=1) - im*(-1 + x[i])*h(i,j,dx=1))) + sin(y[j])*((1 + x[i])*(im*f(i,j,dy=1) + (-1 + x[i])*f(i,j,dx=1))*(im*h(i,j,dy=1) + (-1 + x[i])*h(i,j,dx=1)) - h(i,j)*((4*im)*x[i]*f(i,j,dy=1) - (1 + x[i])*f(i,j,dy=2) + (-1 + x[i])*((1 + 5*x[i])*f(i,j,dx=1) + (1 + x[i])*((2*im)*f(i,j,dx=1,dy=1) + (-1 + x[i])*f(i,j,dx=2)))))) + 2*(1 + x[i])*f(i,j)*g(i,j)*sin(y[j])*(2*(1 + x[i])*h(i,j)*(f(i,j,dy=1) - im*(-1 + x[i])*f(i,j,dx=1))*(g(i,j,dy=1) - im*(-1 + x[i])*g(i,j,dx=1)) - f(i,j)*(-((1 + x[i])*(im*g(i,j,dy=1) + (-1 + x[i])*g(i,j,dx=1))*(im*h(i,j,dy=1) + (-1 + x[i])*h(i,j,dx=1))) + h(i,j)*((4*im)*x[i]*g(i,j,dy=1) - (1 + x[i])*g(i,j,dy=2) + (-1 + x[i])*((1 + 5*x[i])*g(i,j,dx=1) + (1 + x[i])*((2*im)*g(i,j,dx=1,dy=1) + (-1 + x[i])*g(i,j,dx=2)))))) + 4*(-1 + x[i])*g(i,j)^(5/2)*sin(y[j])*(-((1 + x[i])*h(i,j)*sin(y[j])*(im*f(i,j,dy=1) + (-1 + x[i])*f(i,j,dx=1))*(2*W(i,j) + im*W(i,j,dy=1) + (-1 + x[i])*W(i,j,dx=1))) + f(i,j)*(-((1 + x[i])*sin(y[j])*(im*h(i,j,dy=1) + (-1 + x[i])*h(i,j,dx=1))*(2*W(i,j) + im*W(i,j,dy=1) + (-1 + x[i])*W(i,j,dx=1))) + h(i,j)*(((6*im)*(1 + x[i])*cos(y[j]) + 4*sin(y[j]))*W(i,j) + (-3*(1 + x[i])*cos(y[j]) + (2*im)*(2 + x[i])*sin(y[j]))*W(i,j,dy=1) - sin(y[j])*W(i,j,dy=2) - x[i]*sin(y[j])*W(i,j,dy=2) - (3*im)*cos(y[j])*W(i,j,dx=1) + (3*im)*x[i]^2*cos(y[j])*W(i,j,dx=1) - 5*sin(y[j])*W(i,j,dx=1) + 2*x[i]*sin(y[j])*W(i,j,dx=1) + 3*x[i]^2*sin(y[j])*W(i,j,dx=1) - (2*im)*sin(y[j])*W(i,j,dx=1,dy=1) + (2*im)*x[i]^2*sin(y[j])*W(i,j,dx=1,dy=1) + sin(y[j])*W(i,j,dx=2) - x[i]*sin(y[j])*W(i,j,dx=2) - x[i]^2*sin(y[j])*W(i,j,dx=2) + x[i]^3*sin(y[j])*W(i,j,dx=2))))))/(64*rh^2*(1 + x[i])^2*f(i,j)*g(i,j)^3*h(i,j)^2)
end

function psi_2(i::Integer, j::Integer, f::Field,g::Field,h::Field,W::Field,rh::Float64,x::Vector{Float64}=X,y::Vector{Float64}=Y, Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    return ((3*(-1 + x[i])^2*(1 + x[i])^2*sqrt(f(i,j)^9*g(i,j))*h(i,j)^2*(g(i,j,dy=1)^2 + (-1 + x[i])^2*g(i,j,dx=1)^2))/4 + f(i,j)^(5/2)*g(i,j)^3*h(i,j)^2*sin(y[j])*(-((-1 + x[i])^4*sqrt(g(i,j))*sin(y[j])*(4*W(i,j)^2 + W(i,j,dy=1)^2 + 4*(-1 + x[i])*W(i,j)*W(i,j,dx=1) + (-1 + x[i])^2*W(i,j,dx=1)^2)) + (3*im)*(1 - x[i])*(-1 + x[i])^2*(1 + x[i])*(2*W(i,j)*f(i,j,dy=1) - (-1 + x[i])*(W(i,j,dy=1)*f(i,j,dx=1) - f(i,j,dy=1)*W(i,j,dx=1)))) - f(i,j)^(7/2)*g(i,j)^(3/2)*h(i,j)^2*((-1 + x[i]^2)^2*(f(i,j,dy=1)*g(i,j,dy=1) + (-1 + x[i])^2*f(i,j,dx=1)*g(i,j,dx=1)) - (3*im)*(-1 + x[i])^3*g(i,j)^(3/2)*(2*(1 + x[i])*cos(y[j])*W(i,j) + 2*x[i]*sin(y[j])*W(i,j,dy=1) + (-1 + x[i]^2)*cos(y[j])*W(i,j,dx=1)) - ((3*im)/2)*(-1 + x[i])^3*(1 + x[i])*sqrt(g(i,j))*sin(y[j])*(2*W(i,j)*g(i,j,dy=1) - (-1 + x[i])*(W(i,j,dy=1)*g(i,j,dx=1) - g(i,j,dy=1)*W(i,j,dx=1))) - 2*(1 + x[i])*g(i,j)*((-1 + x[i])^4*f(i,j,dx=1) + ((-1 + x[i])^2*(1 + x[i])*(-2*cot(y[j])*f(i,j,dy=1) + f(i,j,dy=2) + (-1 + x[i])*(3*f(i,j,dx=1) + (-1 + x[i])*f(i,j,dx=2))))/2)) - (1 + x[i])*f(i,j)^(9/2)*g(i,j)^(3/2)*(h(i,j)^2*((-1 + x[i])^4*g(i,j,dx=1) + ((-1 + x[i])^2*(1 + x[i])*(-2*cot(y[j])*g(i,j,dy=1) + g(i,j,dy=2) + (-1 + x[i])*(3*g(i,j,dx=1) + (-1 + x[i])*g(i,j,dx=2))))/2) + (-1 + x[i])^2*g(i,j)*(-3*(-1 + x[i])*h(i,j)^2 - (1 + x[i])*(h(i,j,dy=1)^2 + (-1 + x[i])^2*h(i,j,dx=1)^2) + (1 + x[i])*h(i,j)*(h(i,j,dy=2) + (-1 + x[i])*(h(i,j,dx=1) + (-1 + x[i])*h(i,j,dx=2))))))/(48*rh^2*(1 + x[i])^2*(f(i,j)*g(i,j))^(7/2)*h(i,j)^3)
end

function psi_4(i::Integer, j::Integer, f::Field,g::Field,h::Field,W::Field,rh::Float64,x::Vector{Float64}=X,y::Vector{Float64}=Y, Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    return ((-1 + x[i])^2*csc(y[j])*(3*(1 + x[i])^2*f(i,j)^2*h(i,j)*sin(y[j])*((-im)*g(i,j,dy=1) + (-1 + x[i])*g(i,j,dx=1))^2 - 2*(-1 + x[i])*(1 + x[i])*f(i,j)*g(i,j)^(3/2)*h(i,j)*sin(y[j])^2*((-im)*g(i,j,dy=1) + (-1 + x[i])*g(i,j,dx=1))*(2*W(i,j) - im*W(i,j,dy=1) + (-1 + x[i])*W(i,j,dx=1)) - 4*(-1 + x[i])^2*g(i,j)^3*h(i,j)*sin(y[j])^3*(2*W(i,j) - im*W(i,j,dy=1) + (-1 + x[i])*W(i,j,dx=1))^2 - 4*(1 + x[i])*f(i,j)*g(i,j)^2*(-3*(-1 + x[i])*f(i,j)*h(i,j)*sin(y[j]) + f(i,j)*((1 + x[i])*cos(y[j]) - (2*im)*x[i]*sin(y[j]))*(h(i,j,dy=1) + im*(-1 + x[i])*h(i,j,dx=1)) + (1 + x[i])*sin(y[j])*((-im)*f(i,j,dy=1) + (-1 + x[i])*f(i,j,dx=1))*((-im)*h(i,j,dy=1) + (-1 + x[i])*h(i,j,dx=1)) - h(i,j)*sin(y[j])*((-4*im)*x[i]*f(i,j,dy=1) - (1 + x[i])*f(i,j,dy=2) + (-1 + x[i])*((1 + 5*x[i])*f(i,j,dx=1) + (1 + x[i])*((-2*im)*f(i,j,dx=1,dy=1) + (-1 + x[i])*f(i,j,dx=2))))) + 2*(1 + x[i])*f(i,j)*g(i,j)*sin(y[j])*(2*(1 + x[i])*h(i,j)*(f(i,j,dy=1) + im*(-1 + x[i])*f(i,j,dx=1))*(g(i,j,dy=1) + im*(-1 + x[i])*g(i,j,dx=1)) + f(i,j)*((1 + x[i])*((-im)*g(i,j,dy=1) + (-1 + x[i])*g(i,j,dx=1))*((-im)*h(i,j,dy=1) + (-1 + x[i])*h(i,j,dx=1)) - h(i,j)*((-4*im)*x[i]*g(i,j,dy=1) - (1 + x[i])*g(i,j,dy=2) + (-1 + x[i])*((1 + 5*x[i])*g(i,j,dx=1) + (1 + x[i])*((-2*im)*g(i,j,dx=1,dy=1) + (-1 + x[i])*g(i,j,dx=2)))))) + 4*(-1 + x[i])*g(i,j)^(5/2)*sin(y[j])*((1 + x[i])*h(i,j)*sin(y[j])*((-im)*f(i,j,dy=1) + (-1 + x[i])*f(i,j,dx=1))*(2*W(i,j) - im*W(i,j,dy=1) + (-1 + x[i])*W(i,j,dx=1)) + f(i,j)*((1 + x[i])*sin(y[j])*((-im)*h(i,j,dy=1) + (-1 + x[i])*h(i,j,dx=1))*(2*W(i,j) - im*W(i,j,dy=1) + (-1 + x[i])*W(i,j,dx=1)) - h(i,j)*(2*((-3*im)*(1 + x[i])*cos(y[j]) + 2*sin(y[j]))*W(i,j) - (3*(1 + x[i])*cos(y[j]) + (2*im)*(2 + x[i])*sin(y[j]))*W(i,j,dy=1) - sin(y[j])*W(i,j,dy=2) - x[i]*sin(y[j])*W(i,j,dy=2) + (3*im)*cos(y[j])*W(i,j,dx=1) - (3*im)*x[i]^2*cos(y[j])*W(i,j,dx=1) - 5*sin(y[j])*W(i,j,dx=1) + 2*x[i]*sin(y[j])*W(i,j,dx=1) + 3*x[i]^2*sin(y[j])*W(i,j,dx=1) + (2*im)*sin(y[j])*W(i,j,dx=1,dy=1) - (2*im)*x[i]^2*sin(y[j])*W(i,j,dx=1,dy=1) + sin(y[j])*W(i,j,dx=2) - x[i]*sin(y[j])*W(i,j,dx=2) - x[i]^2*sin(y[j])*W(i,j,dx=2) + x[i]^3*sin(y[j])*W(i,j,dx=2))))))/(64*rh^2*(1 + x[i])^2*f(i,j)*g(i,j)^3*h(i,j)^2)
end

function I_Petrov(i::Integer, j::Integer, f::Field,g::Field,h::Field,W::Field,rh::Float64,x::Vector{Float64}=X,y::Vector{Float64}=Y, Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    return 3*psi_2(i,j,f,g,h,W,rh)^2 + psi_0(i,j,f,g,h,W,rh)*psi_4(i,j,f,g,h,W,rh)
end
function J_Petrov(i::Integer, j::Integer, f::Field,g::Field,h::Field,W::Field,rh::Float64,x::Vector{Float64}=X,y::Vector{Float64}=Y, Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    return -psi_2(i,j,f,g,h,W,rh)^3 + psi_0(i,j,f,g,h,W,rh)*psi_2(i,j,f,g,h,W,rh)*psi_4(i,j,f,g,h,W,rh)
end
function L_Petrov(i::Integer, j::Integer, f::Field,g::Field,h::Field,W::Field,rh::Float64,x::Vector{Float64}=X,y::Vector{Float64}=Y, Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    return psi_4(i,j,f,g,h,W,rh)*psi_2(i,j,f,g,h,W,rh)
end
function N_Petrov(i::Integer, j::Integer, f::Field,g::Field,h::Field,W::Field,rh::Float64,x::Vector{Float64}=X,y::Vector{Float64}=Y, Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    return 12*L_Petrov(i,j,f,g,h,W,rh)^2 - psi_4(i,j,f,g,h,W,rh)^2*I_Petrov(i,j,f,g,h,W,rh)
end
function D_Petrov(i::Integer, j::Integer, f::Field,g::Field,h::Field,W::Field,rh::Float64,x::Vector{Float64}=X,y::Vector{Float64}=Y, Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    return 27*J_Petrov(i,j,f,g,h,W,rh)^2-I_Petrov(i,j,f,g,h,W,rh)^3
end
function S_Petrov(i::Integer, j::Integer, f::Field,g::Field,h::Field,W::Field,rh::Float64,x::Vector{Float64}=X,y::Vector{Float64}=Y, Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    return 27*J_Petrov(i,j,f,g,h,W,rh)^2/I_Petrov(i,j,f,g,h,W,rh)^3
end
function Print_Petrov(f::Field,g::Field,h::Field,W::Field,rh::Float64,x::Vector{Float64}=X,y::Vector{Float64}=Y, Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My,Nx::Integer=Nx,Ny::Integer=Ny;file_name::String="Petrov_Type.dat",detailed::Bool=false)
    println("Printing Petrov type data...")
    open(file_name, "w") do io
        for i in 2:Nx-1
            for j in 2:Ny+1
                if detailed
                    write(io,string("x=",X[i],"     ::     y=",Y[j],"\n"))
                    write(io,string("|1-S|=", abs(1-S_Petrov(i,j,f,g,h,W,rh)),"\n"))
                    write(io,string("D = ",D_Petrov(i,j,f,g,h,W,rh),"\n"))
                    write(io,string("I = ",I_Petrov(i,j,f,g,h,W,rh),"\n"))
                    write(io,string("J = ",J_Petrov(i,j,f,g,h,W,rh),"\n"))
                    write(io,string("N = ",N_Petrov(i,j,f,g,h,W,rh),"\n"))
                    write(io,"\n")
                else
                    write(io,string(X[i],"  ",Y[j],"  ",abs(1-S_Petrov(i,j,f,g,h,W,rh)),"\n"))
                end
            end
        end
    end
end

function gtt(x::Float64, y::Float64, f::Field, g::Field, h::Field, W::Field,rh::Float64=1.0;dx::Integer=0,dr::Bool=false)
    if dx==0
        return -1/4*f(x,y)*(1+x)^2 + g(x,y)*sin(y)^2*W(x,y)^2*(1-x)^2/(4*f(x,y))
    elseif dx==1
        if !dr
            return 1/4*(-2*(1+x)*f(x,y)-(1+x)^2*f(x,y,dx=1)-((-1+x)^2*g(x,y)*sin(y)^2*W(x,y)^2*f(x,y,dx=1))/f(x,y)^2+(1/f(x,y))*(-1+x)*sin(y)^2*W(x,y)*((-1+x)*W(x,y)*g(x,y,dx=1)+2*g(x,y)*(W(x,y)+(-1+x)*W(x,y,dx=1))))
        else
            return (1-x)^2/(2*rh) * gtt(x,y,f,g,h,W,rh,dx=1,dr=false)
        end
    else
        if !dr
            return (-2*f(x,y) - 4*(1 + x)*f(x,y,dx=1) + (2*(-1 + x)*sin(y)^2*W(x,y)*(2*f(x,y) - (-1 + x)*f(x,y,dx=1))*(W(x,y)*g(x,y,dx=1) + 2*g(x,y)*W(x,y,dx=1)))/f(x,y)^2 - (1 + x)^2*f(x,y,dx=2) + (g(x,y)*sin(y)^2*W(x,y)^2*(2*(f(x,y) - (-1 + x)*f(x,y,dx=1))^2 - (-1 + x)^2*f(x,y)*f(x,y,dx=2)))/f(x,y)^3 + ((-1 + x)^2*sin(y)^2*(2*g(x,y)*W(x,y,dx=1)^2 + W(x,y)^2*g(x,y,dx=2) + 2*W(x,y)*(2*g(x,y,dx=1)*W(x,y,dx=1) + g(x,y)*W(x,y,dx=2))))/f(x,y))/4
        else
            return (-1+x)^3/(4*rh^2) * ( 2*gtt(x,y,f,g,h,W,rh,dx=1,dr=false) + (-1+x)*gtt(x,y,f,g,h,W,rh,dx=2,dr=false) )
        end
    end
end

function gtphi(x::Float64, y::Float64, f::Field, g::Field, h::Field, W::Field,rh::Float64=1.0;dx::Integer=0,dr::Bool=false)
    if dx==0
        return -rh*g(x,y)*sin(y)^2*W(x,y)/f(x,y)
    elseif dx==1
        if !dr
            return rh*sin(y)^2/f(x,y)^2 * ( -f(x,y)*W(x,y)*g(x,y,dx=1) + g(x,y)*( W(x,y)*f(x,y,dx=1) - f(x,y)*W(x,y,dx=1) ) )
        else
            return (1-x)^2/(2*rh) * gtphi(x,y,f,g,h,W,rh,dx=1,dr=false)
        end
    else
        if !dr
            return -((rh*sin(y)^2*(f(x,y)*(2*g(x,y,dx=1)*(-(W(x,y)*f(x,y,dx=1)) + f(x,y)*W(x,y,dx=1)) + f(x,y)*W(x,y)*g(x,y,dx=2)) + g(x,y)*(W(x,y)*(2*f(x,y,dx=1)^2 - f(x,y)*f(x,y,dx=2)) + f(x,y)*(-2*f(x,y,dx=1)*W(x,y,dx=1) + f(x,y)*W(x,y,dx=2)))))/f(x,y)^3)
        else
            return (-1+x)^3/(4*rh^2) * ( 2*gtphi(x,y,f,g,h,W,rh,dx=1,dr=false) + (-1+x)*gtphi(x,y,f,g,h,W,rh,dx=2,dr=false) )
        end
    end
end

function gphiphi(x::Float64, y::Float64, f::Field, g::Field, h::Field, W::Field,rh::Float64=1.0;dx::Integer=0,dr::Bool=false)
    if dx==0
        return 4*rh^2*g(x,y)*sin(y)^2/((-1+x)^2*f(x,y))
    elseif dx==1
        if !dr
            return -4*rh^2*sin(y)^2/((-1+x)^3*f(x,y)^2) * ( (-1+x)*g(x,y)*f(x,y,dx=1) + f(x,y)*( 2*g(x,y) - (-1+x)*g(x,y,dx=1) ) )
        else
            return (1-x)^2/(2*rh) * gphiphi(x,y,f,g,h,W,rh,dx=1,dr=false)
        end
    else
        if !dr
            return (4*rh^2*sin(y)^2*(2*(-1 + x)^2*g(x,y)*f(x,y,dx=1)^2 + (-1 + x)*f(x,y)*(2*f(x,y,dx=1)*(2*g(x,y) - (-1 + x)*g(x,y,dx=1)) - (-1 + x)*g(x,y)*f(x,y,dx=2)) + f(x,y)^2*(6*g(x,y) + (-1 + x)*(-4*g(x,y,dx=1) + (-1 + x)*g(x,y,dx=2)))))/((-1 + x)^4*f(x,y)^3)
        else
            return (-1+x)^3/(4*rh^2) * ( 2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=false) + (-1+x)*gphiphi(x,y,f,g,h,W,rh,dx=2,dr=false) )
        end
    end
end

function grr(x::Float64, y::Float64, f::Field, g::Field, h::Field, W::Field,rh::Float64=1.0;dx::Integer=0,dr::Bool=false)
    if dx==0
        return (g(x,y)*h(x,y))/f(x,y)
    elseif dx==1
        if !dr
            return (f(x,y)*h(x,y)*g(x,y,dx=1) + g(x,y)*(-(h(x,y)*f(x,y,dx=1)) + f(x,y)*h(x,y,dx=1)))/f(x,y)^2
        else
            return (1-x)^2/(2*rh) * grr(x,y,f,g,h,W,rh,dx=1,dr=false)
        end
    else
        if !dr
            return (f(x,y)*(2*g(x,y,dx=1)*(-(h(x,y)*f(x,y,dx=1)) + f(x,y)*h(x,y,dx=1)) + f(x,y)*h(x,y)*g(x,y,dx=2)) + g(x,y)*(h(x,y)*(2*f(x,y,dx=1)^2 - f(x,y)*f(x,y,dx=2)) + f(x,y)*(-2*f(x,y,dx=1)*h(x,y,dx=1) + f(x,y)*h(x,y,dx=2))))/f(x,y)^3
        else
            return (-1+x)^3/(4*rh^2) * ( 2*grr(x,y,f,g,h,W,rh,dx=1,dr=false) + (-1+x)*grr(x,y,f,g,h,W,rh,dx=2,dr=false) )
        end
    end
end

function gthetatheta(x::Float64, y::Float64, f::Field, g::Field, h::Field, W::Field,rh::Float64=1.0;dx::Integer=0,dr::Bool=false)
    if dx==0
        return (4*rh^2*g(x,y)*h(x,y))/((-1 + x)^2*f(x,y))
    elseif dx==1
        if !dr
            return (4*rh^2*(-((-1 + x)*g(x,y)*h(x,y)*f(x,y,dx=1)) + f(x,y)*((-1 + x)*h(x,y)*g(x,y,dx=1) + g(x,y)*(-2*h(x,y) + (-1 + x)*h(x,y,dx=1)))))/((-1 + x)^3*f(x,y)^2)
        else
            return (1-x)^2/(2*rh) * gthetatheta(x,y,f,g,h,W,rh,dx=1,dr=false)
        end
    else
        if !dr
            return (4*rh^2*(6*f(x,y)^2*g(x,y)*h(x,y) - 4*(-1 + x)*f(x,y)*(f(x,y)*h(x,y)*g(x,y,dx=1) + g(x,y)*(-(h(x,y)*f(x,y,dx=1)) + f(x,y)*h(x,y,dx=1))) - (-1 + x)^2*(f(x,y)*(2*g(x,y,dx=1)*(h(x,y)*f(x,y,dx=1) - f(x,y)*h(x,y,dx=1)) - f(x,y)*h(x,y)*g(x,y,dx=2)) + g(x,y)*(h(x,y)*(-2*f(x,y,dx=1)^2 + f(x,y)*f(x,y,dx=2)) + f(x,y)*(2*f(x,y,dx=1)*h(x,y,dx=1) - f(x,y)*h(x,y,dx=2))))))/((-1 + x)^4*f(x,y)^3)
        else
            return (-1+x)^3/(4*rh^2) * ( 2*gthetatheta(x,y,f,g,h,W,rh,dx=1,dr=false) + (-1+x)*gthetatheta(x,y,f,g,h,W,rh,dx=2,dr=false) )
        end
    end
end

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

function CircumferencialRadius(x::Float64, f::Field, g::Field, h::Field, W::Field, rh::Float64)
    return sqrt(gphiphi(x,pi/2,f,g,h,W,rh))
end

function LightRing(f::Field, g::Field, h::Field, W::Field,rh::Float64,q::Float64=0.0)
    y=pi/2
    M=GetMass(f,g,h,W,rh)
    j=GetJ(f,g,h,W,rh)
    chi = j/M^2
    rhK=M/2*sqrt(1-chi^2-q^2)
    #guess_minus = rBLKerr2x(2*M*(1+cos(2/3*acos(+chi))),rh,chi)
    #guess_plus = rBLKerr2x(2*M*(1+cos(2/3*acos(-chi))),rh,chi)

    rbl_minus=M*find_zero(r->2*q^2+(-3+r)*r-2*chi*sqrt(r-q^2),2*(1+cos(2/3*acos(chi))))
    rbl_plus=M*find_zero(r->2*q^2+(-3+r)*r+2*chi*sqrt(r-q^2),2*(1+cos(2/3*acos(-chi))))

    guess_minus = rBLKerr2x(rbl_minus,rh,chi,q)
    guess_plus = rBLKerr2x(rbl_plus,rh,chi,q)

    factor_minus = x->( -gtphi(x,y,f,g,h,W,rh) + sqrt(gtphi(x,y,f,g,h,W,rh)^2 - gtt(x,y,f,g,h,W,rh)*gphiphi(x,y,f,g,h,W,rh)) )/gtt(x,y,f,g,h,W,rh)
    factor_plus = x->( -gtphi(x,y,f,g,h,W,rh) - sqrt(gtphi(x,y,f,g,h,W,rh)^2 - gtt(x,y,f,g,h,W,rh)*gphiphi(x,y,f,g,h,W,rh)) )/gtt(x,y,f,g,h,W,rh)

    eq_minus = x->gphiphi(x,y,f,g,h,W,rh,dx=1) + 2*gtphi(x,y,f,g,h,W,rh,dx=1)*factor_minus(x) + gtt(x,y,f,g,h,W,rh,dx=1)*factor_minus(x)^2
    eq_plus = x->gphiphi(x,y,f,g,h,W,rh,dx=1) + 2*gtphi(x,y,f,g,h,W,rh,dx=1)*factor_plus(x) + gtt(x,y,f,g,h,W,rh,dx=1)*factor_plus(x)^2

    sol_minus = find_zero(x->eq_minus(x),guess_minus)
    sol_plus = find_zero(x->eq_plus(x),guess_plus)

    #FOR KERR WE HAVE M*w+- = +-1/√(48 cos^4 (1/3 acos(-+ χ)) + χ^2)
    w_minus = (-gtphi(sol_minus,y,f,g,h,W,rh,dx=1,dr=true)-sqrt(gtphi(sol_minus,y,f,g,h,W,rh,dx=1,dr=true)^2-gtt(sol_minus,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(sol_minus,y,f,g,h,W,rh,dx=1,dr=true)))/gphiphi(sol_minus,y,f,g,h,W,rh,dx=1,dr=true)
    w_plus = (-gtphi(sol_plus,y,f,g,h,W,rh,dx=1,dr=true)+sqrt(gtphi(sol_plus,y,f,g,h,W,rh,dx=1,dr=true)^2-gtt(sol_plus,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(sol_plus,y,f,g,h,W,rh,dx=1,dr=true)))/gphiphi(sol_plus,y,f,g,h,W,rh,dx=1,dr=true)

    #wKerr_plus = 1/(M*sqrt(48*cos(1/3*acos(-chi))^4+chi^2))
    #wKerr_minus = -1/(M*sqrt(48*cos(1/3*acos(chi))^4+chi^2))

    wKerr_plus= (rbl_plus^2*sqrt(M*(rbl_plus-M*q^2)) + M^3*q^2*chi - M^2*rbl_plus*chi) / (rbl_plus^4+M^3*(M*q^2-rbl_plus)*chi^2)
    wKerr_minus= (-rbl_minus^2*sqrt(M*(rbl_minus-M*q^2)) + M^3*q^2*chi - M^2*rbl_minus*chi) / (rbl_minus^4+M^3*(M*q^2-rbl_minus)*chi^2)

    RKerr_minus = 2*rhK*sqrt(gKerrN(guess_minus,y,rhK,0.0,chi,q)/((-1+guess_minus)^2*fKerrN(guess_minus,y,rhK,0.0,chi,q)))
    RKerr_plus = 2*rhK*sqrt(gKerrN(guess_plus,y,rhK,0.0,chi,q)/((-1+guess_plus)^2*fKerrN(guess_plus,y,rhK,0.0,chi,q)))

    R_minus=CircumferencialRadius(sol_minus,f,g,h,W,rh)
    R_plus=CircumferencialRadius(sol_plus,f,g,h,W,rh)

    println()
    println("LIGHT RING DETAILS")
    println("Prograde (co-rotating) Circular Photon Orbit Located at x = ", sol_plus, ". r/rh = ", 2/(1-sol_plus), ". Circumferencial Radius/M = ", R_plus/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R_plus/RKerr_plus-1, ". ω*M = ", w_plus*M, ". ω/ωKerrN - 1 = ", w_plus/wKerr_plus-1)
    println("Retrogade (counter-rotating) Circular Photon Orbit Located at x = ", sol_minus, ". r/rh = ", 2/(1-sol_minus), ". Circumferencial Radius/M = ", R_minus/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R_minus/RKerr_minus-1, ". ω*M = ", w_minus*M, ". ω/ωKerrN - 1 = ", w_minus/wKerr_minus-1)
    println()

    return [R_plus/RKerr_plus-1, w_plus/wKerr_plus-1, R_minus/RKerr_minus-1, w_minus/wKerr_minus-1]
end

function ISCO_eq_plus(x::Float64,f::Field, g::Field, h::Field, W::Field,rh::Float64,static=false)
    #ISCO EQUATION TO OBTAIN LOCAITON OF CO-ROTATING ORBIT
    y=pi/2
    if static
        return (2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gphiphi(x,y,f,g,h,W,rh) - 2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh) - gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh) + gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))/(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(x,y,f,g,h,W,rh) - gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh))
    else
        return (-2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^4*gphiphi(x,y,f,g,h,W,rh)^2*gtphi(x,y,f,g,h,W,rh) - 8*gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^5*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)^2 + 8*gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^4*gtphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)*(2*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(x,y,f,g,h,W,rh) + gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh)) + gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + 2*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*(-(gtphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)) + gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh))*sqrt((gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2 - gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true))*(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh) - gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh))^2*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))^2) + 2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*(-(gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh)) + gtphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtt(x,y,f,g,h,W,rh))*sqrt((gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2 - gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true))*(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh) - gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh))^2*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))^2) + gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtphi(x,y,f,g,h,W,rh)*(-2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtt(x,y,f,g,h,W,rh)^2 + (gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh) - 4*gtphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh) + gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtt(x,y,f,g,h,W,rh))*(-gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))) + 2*gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^3*gphiphi(x,y,f,g,h,W,rh)^2*gtphi(x,y,f,g,h,W,rh) - 3*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)^2 + 2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtphi(x,y,f,g,h,W,rh)*(2*gtphi(x,y,f,g,h,W,rh)^2 - 5*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh)*(gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh) + gtphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtt(x,y,f,g,h,W,rh))*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))*(gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtphi(x,y,f,g,h,W,rh) + gtphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh) + 2*gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))) + gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^3*(gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtphi(x,y,f,g,h,W,rh)*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*(-8*gtphi(x,y,f,g,h,W,rh)^3 + 4*gphiphi(x,y,f,g,h,W,rh)*gtphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))) - 2*gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^3*(2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh)*(4*gtphi(x,y,f,g,h,W,rh)^2 - 3*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gphiphi(x,y,f,g,h,W,rh)*(4*gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + gtt(x,y,f,g,h,W,rh)*(gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + gtt(x,y,f,g,h,W,rh)*(gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtt(x,y,f,g,h,W,rh) - gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)))) + gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*(-2*gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtt(x,y,f,g,h,W,rh)*sqrt((gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2 - gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true))*(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh) - gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh))^2*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))^2) + 2*gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*sqrt((gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2 - gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true))*(-(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh)) + gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh))^2*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))^2) + gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtt(x,y,f,g,h,W,rh)^2*(-gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + 2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^3*gphiphi(x,y,f,g,h,W,rh)*(4*gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^2*(-4*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtt(x,y,f,g,h,W,rh)*(-4*gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + (-gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))*(2*gtphi(x,y,f,g,h,W,rh)*(gtphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh) + gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh)) + gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))) + gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*(2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^3*gtt(x,y,f,g,h,W,rh)^3 - gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*(-gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))*(-2*gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh)^2 + gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh) - 6*gtphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh) + gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtt(x,y,f,g,h,W,rh)^2)))) / ((-(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh)) + gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh))*(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gphiphi(x,y,f,g,h,W,rh)^2 + gtt(x,y,f,g,h,W,rh)*(4*gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gphiphi(x,y,f,g,h,W,rh) - 4*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh) + gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtt(x,y,f,g,h,W,rh)) - 2*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*(2*gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtphi(x,y,f,g,h,W,rh) - 2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))))
    end
end

function ISCO_eq_minus(x::Float64,f::Field, g::Field, h::Field, W::Field,rh::Float64,static=false)
    #ISCO EQUATION TO OBTAIN LOCAITON OF COUNTER-ROTATING ORBIT
    y=pi/2
    if static
        return (2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gphiphi(x,y,f,g,h,W,rh) - 2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh) - gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh) + gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))/(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(x,y,f,g,h,W,rh) - gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh))
    else
        return (-2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^4*gphiphi(x,y,f,g,h,W,rh)^2*gtphi(x,y,f,g,h,W,rh) - 8*gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^5*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)^2 + 8*gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^4*gtphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)*(2*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(x,y,f,g,h,W,rh) + gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh)) + 2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*(gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh) - gtphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtt(x,y,f,g,h,W,rh))*sqrt((gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2 - gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true))*(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh) - gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh))^2*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))^2) - gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*(2*(-(gtphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)) + gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh))*sqrt((gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2 - gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true))*(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh) - gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh))^2*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))^2) + gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)*(-gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))) + gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtphi(x,y,f,g,h,W,rh)*(-2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtt(x,y,f,g,h,W,rh)^2 + (gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh) - 4*gtphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh) + gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtt(x,y,f,g,h,W,rh))*(-gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))) + 2*gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^3*gphiphi(x,y,f,g,h,W,rh)^2*gtphi(x,y,f,g,h,W,rh) - 3*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)^2 + 2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtphi(x,y,f,g,h,W,rh)*(2*gtphi(x,y,f,g,h,W,rh)^2 - 5*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh)*(gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh) + gtphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtt(x,y,f,g,h,W,rh))*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))*(gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtphi(x,y,f,g,h,W,rh) + gtphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh) + 2*gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))) + gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^3*(gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtphi(x,y,f,g,h,W,rh)*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*(-8*gtphi(x,y,f,g,h,W,rh)^3 + 4*gphiphi(x,y,f,g,h,W,rh)*gtphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))) - 2*gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^3*(2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh)*(4*gtphi(x,y,f,g,h,W,rh)^2 - 3*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gphiphi(x,y,f,g,h,W,rh)*(4*gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + gtt(x,y,f,g,h,W,rh)*(gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + gtt(x,y,f,g,h,W,rh)*(gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtt(x,y,f,g,h,W,rh) - gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)))) + gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*(2*gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtt(x,y,f,g,h,W,rh)*sqrt((gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2 - gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true))*(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh) - gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh))^2*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))^2) - 2*gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*sqrt((gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2 - gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true))*(-(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh)) + gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh))^2*(gtphi(x,y,f,g,h,W,rh)^2 - gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))^2) + gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtt(x,y,f,g,h,W,rh)^2*(-gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + 2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^3*gphiphi(x,y,f,g,h,W,rh)*(4*gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^2*(-4*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtt(x,y,f,g,h,W,rh)*(-4*gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh)) + (-gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))*(2*gtphi(x,y,f,g,h,W,rh)*(gtphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh) + gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh)) + gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))) + gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*(2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^3*gtt(x,y,f,g,h,W,rh)^3 - gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*(-gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))*(-2*gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh)^2 + gtt(x,y,f,g,h,W,rh,dx=2,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh) - 6*gtphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh) + gphiphi(x,y,f,g,h,W,rh,dx=2,dr=true)*gtt(x,y,f,g,h,W,rh)^2))))/((-(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh)) + gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtt(x,y,f,g,h,W,rh))*(gtt(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gphiphi(x,y,f,g,h,W,rh)^2 + gtt(x,y,f,g,h,W,rh)*(4*gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gphiphi(x,y,f,g,h,W,rh) - 4*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh) + gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)^2*gtt(x,y,f,g,h,W,rh)) - 2*gtt(x,y,f,g,h,W,rh,dx=1,dr=true)*(2*gtphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtphi(x,y,f,g,h,W,rh) - 2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gtphi(x,y,f,g,h,W,rh)^2 + gphiphi(x,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(x,y,f,g,h,W,rh)*gtt(x,y,f,g,h,W,rh))))
    end
end

function ISCO(f::Field, g::Field, h::Field, W::Field,rh::Float64,q::Float64=0.0)
    y=pi/2
    M=GetMass(f,g,h,W,rh)
    j=GetJ(f,g,h,W,rh)
    chi = j/M^2
    rhK=M/2*sqrt(1-chi^2-q^2)

    static=false
    if chi < 1e-10
        static=true
    end

    #SEE https://duetosymmetry.com/tool/kerr-isco-calculator/
    z1=1+(1-chi^2)^(1/3) * ( (1+chi)^(1/3) + (1-chi)^(1/3) )
    z2=sqrt(3*chi^2+z1^2)
    rbl_minus = (3+z2+sqrt((3-z1)*(3+z1+2*z2)))
    rbl_plus= (3+z2-sqrt((3-z1)*(3+z1+2*z2)))

    rbl_minus=M*find_zero(r->4*(q^2-r)*(sqrt(r-q^2)+chi)^2+r*(q^2+(r-2)*r+chi^2),rbl_minus)
    rbl_plus=M*find_zero(r->4*(q^2-r)*(sqrt(r-q^2)-chi)^2+r*(q^2+(r-2)*r+chi^2),rbl_plus)

    guess_plus=rBLKerr2x(rbl_plus,rhK,chi,q)
    guess_minus=rBLKerr2x(rbl_minus,rhK,chi,q)

    sol_plus = fzero(x->ISCO_eq_plus(x,f,g,h,W,rh,static),guess_plus)
    sol_minus = fzero(x->ISCO_eq_minus(x,f,g,h,W,rh,static),guess_minus)

    #FOR KERR WE HAVE M*w+- = +- 1/((rBL+-/M)^(3/2) +- χ). SCHWARZSCHILD = 1/(6*√6)
    w_plus = (-gtphi(sol_plus,y,f,g,h,W,rh,dx=1,dr=true)+sqrt(gtphi(sol_plus,y,f,g,h,W,rh,dx=1,dr=true)^2-gtt(sol_plus,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(sol_plus,y,f,g,h,W,rh,dx=1,dr=true)))/gphiphi(sol_plus,y,f,g,h,W,rh,dx=1,dr=true)
    w_minus = (-gtphi(sol_minus,y,f,g,h,W,rh,dx=1,dr=true)-sqrt(gtphi(sol_minus,y,f,g,h,W,rh,dx=1,dr=true)^2-gtt(sol_minus,y,f,g,h,W,rh,dx=1,dr=true)*gphiphi(sol_minus,y,f,g,h,W,rh,dx=1,dr=true)))/gphiphi(sol_minus,y,f,g,h,W,rh,dx=1,dr=true)

    #wKerr_plus = 1/((rbl_plus/M)^(3/2)+chi)/M
    #wKerr_minus = -1/((rbl_minus/M)^(3/2)-chi)/M

    wKerr_plus= (rbl_plus^2*sqrt(M*(rbl_plus-M*q^2)) + M^3*q^2*chi - M^2*rbl_plus*chi) / (rbl_plus^4+M^3*(M*q^2-rbl_plus)*chi^2)
    wKerr_minus= (-rbl_minus^2*sqrt(M*(rbl_minus-M*q^2)) + M^3*q^2*chi - M^2*rbl_minus*chi) / (rbl_minus^4+M^3*(M*q^2-rbl_minus)*chi^2)

    RKerr_plus = 2*rhK*sqrt(gKerrN(guess_plus,y,rhK,0.0,chi,q)/((-1+guess_plus)^2*fKerrN(guess_plus,y,rhK,0.0,chi,q)))
    RKerr_minus = 2*rhK*sqrt(gKerrN(guess_minus,y,rhK,0.0,chi,q)/((-1+guess_minus)^2*fKerrN(guess_minus,y,rhK,0.0,chi,q)))

    R_plus = CircumferencialRadius(sol_plus,f,g,h,W,rh)
    R_minus = CircumferencialRadius(sol_minus,f,g,h,W,rh)

    println()
    println("ISCO DETAILS")
    println("Prograde (co-rotating) Circular Orbit of Massive Particles Located at x = ", sol_plus, ". r/rh = ", 2/(1-sol_plus), ". Circumferencial Radius/M = ", R_plus/M, ". Coordinate shift w.r.t. KerrN (x-xKerrN) = ", sol_plus-guess_plus, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R_plus/RKerr_plus-1,  ". ω*M = ", w_plus*M, ". ω/ωKerrN - 1 = ", w_plus/wKerr_plus-1)
    println("Retrogade (counter-rotating) Circular Orbit of Massive Particles Located at x = ", sol_minus, ". r/rh = ", 2/(1-sol_minus), ". Circumferencial Radius/M = ", R_minus/M, ". Coordinate shift w.r.t. KerrN (x-xKerrN) = ", sol_minus-guess_minus, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R_minus/RKerr_minus-1, ". ω*M = ", w_minus*M, ". ω/ωKerrN - 1 = ", w_minus/wKerr_minus-1)
    println()

    return [R_plus/RKerr_plus-1, w_plus/wKerr_plus-1, R_minus/RKerr_minus-1, w_minus/wKerr_minus-1]
end

function Ricci_Horizon(y::Float64,rh::Float64,f::Field, g::Field, h::Field, W::Field)
    return (g(-1.0,y)^2*h(-1.0,y)^2*(4*g(-1.0,y)*sin(y)^2*(3*W(-1.0,y) - 2*W(-1.0,y,dx=2))^2 - 3*f(-1.0,y,dy=1)^2) + f(-1.0,y)*g(-1.0,y)*h(-1.0,y)^2*(f(-1.0,y,dy=1)*(2*cot(y)*g(-1.0,y) + g(-1.0,y,dy=1)) + 2*g(-1.0,y)*f(-1.0,y,dy=2)) + f(-1.0,y)^2*(3*h(-1.0,y)^2*g(-1.0,y,dy=1)^2 - 4*g(-1.0,y)*h(-1.0,y)^2*(6*g(-1.0,y,dx=2) + cot(y)*g(-1.0,y,dy=1) + g(-1.0,y,dy=2)) + g(-1.0,y)^2*(2*h(-1.0,y,dy=1)^2 + h(-1.0,y)*(3*h(-1.0,y) - 8*h(-1.0,y,dx=2) - 2*h(-1.0,y,dy=2)))))/(2*rh^2*f(-1.0,y)*g(-1.0,y)^3*h(-1.0,y)^3)
end

function GB_Horizon(y::Float64,rh::Float64,f::Field, g::Field, h::Field, W::Field)
    return (12*g(-1.0,y)^3*h(-1.0,y)^2*f(-1.0,y,dy=1)^2*(-12*g(-1.0,y)*sin(y)^2*(3*W(-1.0,y) - 2*W(-1.0,y,dx=2))^2 - f(-1.0,y,dy=1)^2) - 3*f(-1.0,y)*g(-1.0,y)^2*h(-1.0,y)*(16*f(-1.0,y,dx=2)*g(-1.0,y)*h(-1.0,y)*(4*g(-1.0,y)*sin(y)^2*(3*W(-1.0,y) - 2*W(-1.0,y,dx=2))^2 - f(-1.0,y,dy=1)^2) - 3*h(-1.0,y)*f(-1.0,y,dy=1)^3*g(-1.0,y,dy=1) - 2*g(-1.0,y)*h(-1.0,y)*f(-1.0,y,dy=1)*(3*cot(y)*f(-1.0,y,dy=1)^2 + 14*sin(y)^2*(3*W(-1.0,y) - 2*W(-1.0,y,dx=2))^2*g(-1.0,y,dy=1) + 3*f(-1.0,y,dy=1)*f(-1.0,y,dy=2)) - 8*g(-1.0,y)^2*sin(y)*(3*W(-1.0,y) - 2*W(-1.0,y,dx=2))*(4*sin(y)*(-3*W(-1.0,y) + 2*W(-1.0,y,dx=2))*f(-1.0,y,dy=1)*h(-1.0,y,dy=1) + h(-1.0,y)*(f(-1.0,y,dy=1)*(45*W(-1.0,y)*cos(y) - 30*cos(y)*W(-1.0,y,dx=2) - 8*sin(y)*W(-1.0,y,dx=2,dy=1)) + 3*sin(y)*(3*W(-1.0,y) - 2*W(-1.0,y,dx=2))*f(-1.0,y,dy=2)))) + f(-1.0,y)^2*g(-1.0,y)*(384*f(-1.0,y,dx=2)^2*g(-1.0,y)^2*h(-1.0,y)^2 + 18*h(-1.0,y)^2*f(-1.0,y,dy=1)^2*g(-1.0,y,dy=1)^2 - 96*f(-1.0,y,dx=2)*g(-1.0,y)^2*h(-1.0,y)*(-(f(-1.0,y,dy=1)*h(-1.0,y,dy=1)) + 2*h(-1.0,y)*(cot(y)*f(-1.0,y,dy=1) + f(-1.0,y,dy=2))) - 3*g(-1.0,y)*h(-1.0,y)*(32*g(-1.0,y,dx=2)*h(-1.0,y)*f(-1.0,y,dy=1)^2 - 5*f(-1.0,y,dy=1)^2*g(-1.0,y,dy=1)*h(-1.0,y,dy=1) + 2*h(-1.0,y)*(-4*sin(y)^2*(3*W(-1.0,y) - 2*W(-1.0,y,dx=2))^2*g(-1.0,y,dy=1)^2 + 4*f(-1.0,y,dy=1)*g(-1.0,y,dy=1)*f(-1.0,y,dy=2) + f(-1.0,y,dy=1)^2*(cot(y)*g(-1.0,y,dy=1) + 2*g(-1.0,y,dy=2)))) + 6*g(-1.0,y)^2*(32*g(-1.0,y,dx=2)*h(-1.0,y)^2*sin(y)^2*(3*W(-1.0,y) - 2*W(-1.0,y,dx=2))^2 + 2*f(-1.0,y,dy=1)^2*h(-1.0,y,dy=1)^2 - h(-1.0,y)^2*(17*f(-1.0,y,dy=1)^2 + f(-1.0,y,dy=1)*(-32*f(-1.0,y,dx=2,dy=1) + 4*cot(y)*f(-1.0,y,dy=2)) + 4*sin(y)*(3*W(-1.0,y) - 2*W(-1.0,y,dx=2))*(g(-1.0,y,dy=1)*(27*W(-1.0,y)*cos(y) - 18*cos(y)*W(-1.0,y,dx=2) - 4*sin(y)*W(-1.0,y,dx=2,dy=1)) + 2*sin(y)*(3*W(-1.0,y) - 2*W(-1.0,y,dx=2))*g(-1.0,y,dy=2))) + h(-1.0,y)*(-24*h(-1.0,y,dx=2)*f(-1.0,y,dy=1)^2 + 10*sin(y)^2*(3*W(-1.0,y) - 2*W(-1.0,y,dx=2))^2*g(-1.0,y,dy=1)*h(-1.0,y,dy=1) - 2*f(-1.0,y,dy=1)*h(-1.0,y,dy=1)*f(-1.0,y,dy=2) + f(-1.0,y,dy=1)^2*(cot(y)*h(-1.0,y,dy=1) - h(-1.0,y,dy=2)))) + 4*g(-1.0,y)^3*(-12*sin(y)^2*(3*W(-1.0,y) - 2*W(-1.0,y,dx=2))^2*h(-1.0,y,dy=1)^2 - 6*h(-1.0,y)^2*(9*W(-1.0,y)^2*(1 + 11*cos(2*y)) + (4 + 44*cos(2*y))*W(-1.0,y,dx=2)^2 - 36*W(-1.0,y)*sin(2*y)*W(-1.0,y,dx=2,dy=1) + 12*sin(y)^2*W(-1.0,y,dx=2,dy=1)^2 - 12*W(-1.0,y,dx=2)*(W(-1.0,y) + 11*W(-1.0,y)*cos(2*y) - 2*sin(2*y)*W(-1.0,y,dx=2,dy=1))) + 6*h(-1.0,y)*sin(y)*(3*W(-1.0,y) - 2*W(-1.0,y,dx=2))*(8*h(-1.0,y,dx=2)*sin(y)*(3*W(-1.0,y) - 2*W(-1.0,y,dx=2)) + h(-1.0,y,dy=1)*(27*W(-1.0,y)*cos(y) - 18*cos(y)*W(-1.0,y,dx=2) - 8*sin(y)*W(-1.0,y,dx=2,dy=1)) + sin(y)*(3*W(-1.0,y) - 2*W(-1.0,y,dx=2))*h(-1.0,y,dy=2)))) + 3*f(-1.0,y)^4*(-40*g(-1.0,y,dx=2)*h(-1.0,y)^2*g(-1.0,y,dy=1)^2 + 2*g(-1.0,y)*h(-1.0,y)*(64*g(-1.0,y,dx=2)^2*h(-1.0,y) + g(-1.0,y,dy=1)*(4*h(-1.0,y,dx=2)*g(-1.0,y,dy=1) + h(-1.0,y)*(9*g(-1.0,y,dy=1) + 16*g(-1.0,y,dx=2,dy=1))) - 24*g(-1.0,y,dx=2)*g(-1.0,y,dy=1)*h(-1.0,y,dy=1)) + 2*g(-1.0,y)^3*(46*h(-1.0,y)^2 - 8*h(-1.0,y,dy=1)^2 + h(-1.0,y)*(48*h(-1.0,y,dx=2) + 19*cot(y)*h(-1.0,y,dy=1) + 4*h(-1.0,y,dy=2))) + g(-1.0,y)^2*(-(h(-1.0,y)*(-((3*g(-1.0,y,dy=1) + 32*g(-1.0,y,dx=2,dy=1))*h(-1.0,y,dy=1)) + 22*h(-1.0,y)*(cot(y)*g(-1.0,y,dy=1) + g(-1.0,y,dy=2)) + 16*h(-1.0,y,dx=2)*(2*cot(y)*g(-1.0,y,dy=1) + g(-1.0,y,dy=2)))) + 16*g(-1.0,y,dx=2)*(10*h(-1.0,y)^2 - 2*h(-1.0,y,dy=1)^2 + h(-1.0,y)*(8*h(-1.0,y,dx=2) + 2*cot(y)*h(-1.0,y,dy=1) + h(-1.0,y,dy=2))))) - 3*f(-1.0,y)^3*(5*h(-1.0,y)^2*f(-1.0,y,dy=1)*g(-1.0,y,dy=1)^3 + g(-1.0,y)*h(-1.0,y)*g(-1.0,y,dy=1)*(-96*g(-1.0,y,dx=2)*h(-1.0,y)*f(-1.0,y,dy=1) + 5*f(-1.0,y,dy=1)*g(-1.0,y,dy=1)*h(-1.0,y,dy=1) + h(-1.0,y)*(-2*g(-1.0,y,dy=1)*f(-1.0,y,dy=2) + 4*f(-1.0,y,dy=1)*(cot(y)*g(-1.0,y,dy=1) - g(-1.0,y,dy=2)))) + 8*f(-1.0,y,dx=2)*g(-1.0,y)*(5*h(-1.0,y)^2*g(-1.0,y,dy=1)^2 + g(-1.0,y)*h(-1.0,y)*(32*g(-1.0,y,dx=2)*h(-1.0,y) + g(-1.0,y,dy=1)*h(-1.0,y,dy=1) - 6*h(-1.0,y)*(cot(y)*g(-1.0,y,dy=1) + g(-1.0,y,dy=2))) + 2*g(-1.0,y)^2*(16*h(-1.0,y)^2 - 2*h(-1.0,y,dy=1)^2 + h(-1.0,y)*(8*h(-1.0,y,dx=2) + 5*cot(y)*h(-1.0,y,dy=1) + h(-1.0,y,dy=2)))) - g(-1.0,y)^3*(-8*cot(y)*f(-1.0,y,dy=1)*h(-1.0,y,dy=1)^2 + 30*h(-1.0,y)^2*(cot(y)*f(-1.0,y,dy=1) + f(-1.0,y,dy=2)) + h(-1.0,y)*(16*h(-1.0,y,dx=2)*(5*cot(y)*f(-1.0,y,dy=1) + f(-1.0,y,dy=2)) + 4*h(-1.0,y,dy=1)*(-8*f(-1.0,y,dx=2,dy=1) + cot(y)*f(-1.0,y,dy=2)) + f(-1.0,y,dy=1)*(-3*h(-1.0,y,dy=1) + 4*cot(y)*h(-1.0,y,dy=2)))) - 2*g(-1.0,y)^2*(-2*f(-1.0,y,dy=1)*g(-1.0,y,dy=1)*h(-1.0,y,dy=1)^2 + 4*g(-1.0,y,dx=2)*h(-1.0,y)*(3*f(-1.0,y,dy=1)*h(-1.0,y,dy=1) + 2*h(-1.0,y)*(cot(y)*f(-1.0,y,dy=1) + f(-1.0,y,dy=2))) + 2*h(-1.0,y)^2*(g(-1.0,y,dy=1)*(-8*f(-1.0,y,dx=2,dy=1) + cot(y)*f(-1.0,y,dy=2)) + f(-1.0,y,dy=1)*(3*g(-1.0,y,dy=1) - 16*g(-1.0,y,dx=2,dy=1) + cot(y)*g(-1.0,y,dy=2))) + h(-1.0,y)*(20*h(-1.0,y,dx=2)*f(-1.0,y,dy=1)*g(-1.0,y,dy=1) + g(-1.0,y,dy=1)*h(-1.0,y,dy=1)*f(-1.0,y,dy=2) + f(-1.0,y,dy=1)*(h(-1.0,y,dy=1)*g(-1.0,y,dy=2) + g(-1.0,y,dy=1)*(-2*cot(y)*h(-1.0,y,dy=1) + h(-1.0,y,dy=2)))))))/(6*rh^4*f(-1.0,y)^2*g(-1.0,y)^5*h(-1.0,y)^4)
end

function RicciScalar(x::Float64, y::Float64,rh::Float64,f::Field, g::Field, h::Field, W::Field)
    return (3*(-1 + x)^2*(1 + x)^2*f(x,y)^2*h(x,y)^2*(g(x,y,dy=1)^2 + (-1 + x)^2*g(x,y,dx=1)^2) + (-1 + x)^4*g(x,y)^3*h(x,y)^2*sin(y)^2*(4*W(x,y)^2 + W(x,y,dy=1)^2 + 4*(-1 + x)*W(x,y)*W(x,y,dx=1) + (-1 + x)^2*W(x,y,dx=1)^2) + 2*(1 + x)*f(x,y)*g(x,y)*h(x,y)^2*(((-1 + x)^2*(1 + x)*(f(x,y,dy=1)*g(x,y,dy=1) + (-1 + x)^2*f(x,y,dx=1)*g(x,y,dx=1)))/2 - f(x,y)*((-1 + x)^4*g(x,y,dx=1) + 2*(-1 + x)^2*(1 + x)*(cot(y)*g(x,y,dy=1) + g(x,y,dy=2) + (-1 + x)^2*g(x,y,dx=2)))) - 2*(1 + x)*g(x,y)^2*((3*(-1 + x)^2*(1 + x)*h(x,y)^2*(f(x,y,dy=1)^2 + (-1 + x)^2*f(x,y,dx=1)^2))/2 + f(x,y)*h(x,y)^2*((-1 + x)^4*f(x,y,dx=1) - (-1 + x)^2*(1 + x)*(cot(y)*f(x,y,dy=1) + f(x,y,dy=2) + (-1 + x)^2*f(x,y,dx=2))) - (-1 + x)^2*(1 + x)*f(x,y)^2*(h(x,y,dy=1)^2 + (-1 + x)^2*h(x,y,dx=1)^2 - h(x,y)*(h(x,y,dy=2) + (-1 + x)*(h(x,y,dx=1) + (-1 + x)*h(x,y,dx=2))))))/(8*rh^2*(1 + x)^2*f(x,y)*g(x,y)^3*h(x,y)^3)
end

function GBScalar(x::Float64, y::Float64,rh::Float64,f::Field, g::Field, h::Field, W::Field)
    return ((-1 + x)^2*(-4*rh^2*(-1 + x)^2*(1 + x)^2*g(x,y)^3*h(x,y)^2*((1 + x)^2*(f(x,y,dy=1)^2 + (-1 + x)^2*f(x,y,dx=1)^2)^2 + (1 - x)*g(x,y)*sin(y)^2*(4*(1 - x)*W(x,y)^2*(3*f(x,y,dy=1)^2 - (-1 + x)^2*f(x,y,dx=1)^2) + 8*(-1 + x)^3*f(x,y,dy=1)*W(x,y,dy=1)*f(x,y,dx=1)*W(x,y,dx=1) + 4*(-1 + x)^2*W(x,y)*(4*f(x,y,dy=1)*W(x,y,dy=1)*f(x,y,dx=1) - 3*f(x,y,dy=1)^2*W(x,y,dx=1) + (-1 + x)^2*f(x,y,dx=1)^2*W(x,y,dx=1)) - (1 - x)*f(x,y,dy=1)^2*(W(x,y,dy=1)^2 - 3*(-1 + x)^2*W(x,y,dx=1)^2) + (1 - x)^3*f(x,y,dx=1)^2*(3*W(x,y,dy=1)^2 - (-1 + x)^2*W(x,y,dx=1)^2))) - 2*rh^2*(1 - x)*(1 + x)^3*f(x,y)^4*(5*(1 - x)^3*h(x,y)^2*g(x,y,dx=1)*(g(x,y,dy=1)^2 + (-1 + x)^2*g(x,y,dx=1)^2) + (1 - x)*g(x,y)*h(x,y)*(6*(-1 + x)^2*g(x,y,dy=1)*h(x,y,dy=1)*g(x,y,dx=1) - 4*(1 - x)*h(x,y)*(2*g(x,y,dy=1)^2 - (-1 + x)^2*g(x,y,dx=1)^2) - (-1 + x)^2*g(x,y,dy=1)^2*h(x,y,dx=1) + 5*(-1 + x)^4*g(x,y,dx=1)^2*h(x,y,dx=1) + 2*h(x,y)*(-3*(-1 + x)*g(x,y,dy=1)^2 - 2*(-1 + x)^2*g(x,y,dy=1)*g(x,y,dx=1,dy=1) + 2*(1 - x)^3*g(x,y,dx=1)*(4*g(x,y,dx=1) + (-1 + x)*g(x,y,dx=2)))) + 2*g(x,y)^2*(2*(1 - x)^3*g(x,y,dx=1)*(h(x,y,dy=1)^2 + (-1 + x)^2*h(x,y,dx=1)^2) - 2*(-1 + x)^2*h(x,y)^2*(-(cot(y)*g(x,y,dy=1)) - g(x,y,dy=2) + 2*(-1 + x)*(3*g(x,y,dx=1) + (-1 + x)*g(x,y,dx=2))) - (1 - x)*h(x,y)*((-1 + x)^2*h(x,y,dy=2)*g(x,y,dx=1) - (-1 + x)^2*g(x,y,dy=2)*h(x,y,dx=1) + 3*(-1 + x)^3*g(x,y,dx=1)*h(x,y,dx=1) - (1 - x)*g(x,y,dy=1)*(3*h(x,y,dy=1) - 2*(-1 + x)*cot(y)*h(x,y,dx=1)) + 2*(1 - x)*(g(x,y,dy=1)*h(x,y,dy=1) - (-1 + x)^2*g(x,y,dx=1)*h(x,y,dx=1)) + 2*(-1 + x)^2*h(x,y,dy=1)*(cot(y)*g(x,y,dx=1) + g(x,y,dx=1,dy=1)) - (1 - x)^3*h(x,y,dx=1)*(2*g(x,y,dx=1) + (-1 + x)*g(x,y,dx=2)) - (1 - x)^3*g(x,y,dx=1)*(2*h(x,y,dx=1) + (-1 + x)*h(x,y,dx=2)))) - 4*(-1 + x)^2*g(x,y)^3*(-2*(h(x,y,dy=1)^2 + (-1 + x)^2*h(x,y,dx=1)^2) + h(x,y)*(3*cot(y)*h(x,y,dy=1) + h(x,y,dy=2) + (-1 + x)*(4*h(x,y,dx=1) + (-1 + x)*h(x,y,dx=2))))) - 2*rh^2*(1 + x)^3*f(x,y)^3*((5*(-1 + x)^2*(1 + x)*h(x,y)^2*(f(x,y,dy=1)*g(x,y,dy=1) + (-1 + x)^2*f(x,y,dx=1)*g(x,y,dx=1))*(g(x,y,dy=1)^2 + (-1 + x)^2*g(x,y,dx=1)^2))/2 + g(x,y)*h(x,y)*(((-1 + x)^2*(1 + x)*(f(x,y,dy=1)*(5*g(x,y,dy=1)^2*h(x,y,dy=1) - (-1 + x)^2*h(x,y,dy=1)*g(x,y,dx=1)^2 + 6*(-1 + x)^2*g(x,y,dy=1)*g(x,y,dx=1)*h(x,y,dx=1)) + (-1 + x)^2*f(x,y,dx=1)*(6*g(x,y,dy=1)*h(x,y,dy=1)*g(x,y,dx=1) - g(x,y,dy=1)^2*h(x,y,dx=1) + 5*(-1 + x)^2*g(x,y,dx=1)^2*h(x,y,dx=1))))/2 + h(x,y)*(3*(-1 + x)^4*(g(x,y,dy=1)^2*f(x,y,dx=1) - 4*f(x,y,dy=1)*g(x,y,dy=1)*g(x,y,dx=1) - 3*(-1 + x)^2*f(x,y,dx=1)*g(x,y,dx=1)^2) - (-1 + x)^2*(1 + x)*(-2*f(x,y,dy=1)*(cot(y)*g(x,y,dy=1)^2 - g(x,y,dy=1)*(g(x,y,dy=2) + 3*(-1 + x)*g(x,y,dx=1)) + (-1 + x)^2*g(x,y,dx=1)*(cot(y)*g(x,y,dx=1) - g(x,y,dx=1,dy=1))) + 2*(-1 + x)^2*g(x,y,dy=1)*(3*g(x,y,dx=1)*f(x,y,dx=1,dy=1) + f(x,y,dx=1)*g(x,y,dx=1,dy=1)) + g(x,y,dy=1)^2*(f(x,y,dy=2) - (-1 + x)*(f(x,y,dx=1) + 2*(-1 + x)*f(x,y,dx=2))) - g(x,y,dx=1)*(2*(-1 + x)^2*f(x,y,dy=2)*g(x,y,dx=1) + (1 - x)^3*g(x,y,dx=1)*(2*f(x,y,dx=1) + (-1 + x)*f(x,y,dx=2)) + 2*(1 - x)^3*f(x,y,dx=1)*(4*g(x,y,dx=1) + (-1 + x)*g(x,y,dx=2)))))) + 2*g(x,y)^2*((-1 + x)^2*(1 + x)*(f(x,y,dy=1)*g(x,y,dy=1) + (-1 + x)^2*f(x,y,dx=1)*g(x,y,dx=1))*(h(x,y,dy=1)^2 + (-1 + x)^2*h(x,y,dx=1)^2) - 2*h(x,y)^2*((-1 + x)^2*((-1 + x)^2*g(x,y,dy=2)*f(x,y,dx=1) - 2*(-1 + x)^3*f(x,y,dx=1)*g(x,y,dx=1) + (-1 + x)^2*g(x,y,dy=1)*(cot(y)*f(x,y,dx=1) - f(x,y,dx=1,dy=1)) + 2*(1 - x)*f(x,y,dy=1)*(g(x,y,dy=1) + (-1 + x)*g(x,y,dx=1,dy=1)) + (1 - x)^3*g(x,y,dx=1)*(2*f(x,y,dx=1) + (-1 + x)*f(x,y,dx=2)) + (1 - x)^3*f(x,y,dx=1)*(4*g(x,y,dx=1) + (-1 + x)*g(x,y,dx=2))) - ((1 + x)*((1 - x)^3*g(x,y,dy=2)*(f(x,y,dx=1) + (-1 + x)*f(x,y,dx=2)) - (-1 + x)^2*g(x,y,dy=1)*(cot(y)*f(x,y,dy=2) + (-1 + x)*(cot(y)*f(x,y,dx=1) - 2*f(x,y,dx=1,dy=1) + (-1 + x)*cot(y)*f(x,y,dx=2))) + (1 - x)^3*f(x,y,dy=2)*(g(x,y,dx=1) + (-1 + x)*g(x,y,dx=2)) - 2*(1 - x)*(-1 + x)^3*(f(x,y,dx=1,dy=1)*g(x,y,dx=1,dy=1) + (-1 + x)*g(x,y,dx=1)*f(x,y,dx=2) + f(x,y,dx=1)*(3*g(x,y,dx=1) + (-1 + x)*g(x,y,dx=2))) - (-1 + x)^2*f(x,y,dy=1)*(-2*g(x,y,dy=1) + cot(y)*g(x,y,dy=2) + (-1 + x)*(cot(y)*g(x,y,dx=1) - 2*g(x,y,dx=1,dy=1) + (-1 + x)*cot(y)*g(x,y,dx=2)))))/2) - h(x,y)*(4*(-1 + x)^6*f(x,y,dx=1)*g(x,y,dx=1)*h(x,y,dx=1) + 2*(-1 + x)^4*f(x,y,dy=1)*(h(x,y,dy=1)*g(x,y,dx=1) + g(x,y,dy=1)*h(x,y,dx=1)) + ((-1 + x)^2*(1 + x)*((-1 + x)^2*h(x,y,dy=2)*f(x,y,dx=1)*g(x,y,dx=1) - (-1 + x)^2*g(x,y,dy=2)*f(x,y,dx=1)*h(x,y,dx=1) - (-1 + x)^2*f(x,y,dy=2)*g(x,y,dx=1)*h(x,y,dx=1) + 9*(-1 + x)^3*f(x,y,dx=1)*g(x,y,dx=1)*h(x,y,dx=1) - 2*(-1 + x)^2*g(x,y,dy=1)*h(x,y,dx=1)*(cot(y)*f(x,y,dx=1) - f(x,y,dx=1,dy=1)) + 2*(-1 + x)^2*h(x,y,dy=1)*g(x,y,dx=1)*f(x,y,dx=1,dy=1) - 2*(-1 + x)^2*f(x,y,dy=1)*h(x,y,dx=1)*(cot(y)*g(x,y,dx=1) - g(x,y,dx=1,dy=1)) + 2*(-1 + x)^2*h(x,y,dy=1)*f(x,y,dx=1)*(cot(y)*g(x,y,dx=1) + g(x,y,dx=1,dy=1)) + (-1 + x)^4*g(x,y,dx=1)*h(x,y,dx=1)*f(x,y,dx=2) + g(x,y,dy=1)*h(x,y,dy=1)*(f(x,y,dy=2) + (-1 + x)*(f(x,y,dx=1) - (-1 + x)*f(x,y,dx=2))) + (-1 + x)^4*f(x,y,dx=1)*h(x,y,dx=1)*g(x,y,dx=2) + f(x,y,dy=1)*h(x,y,dy=1)*(g(x,y,dy=2) + (-1 + x)*(g(x,y,dx=1) - (-1 + x)*g(x,y,dx=2))) + (-1 + x)^4*f(x,y,dx=1)*g(x,y,dx=1)*h(x,y,dx=2) - f(x,y,dy=1)*g(x,y,dy=1)*(2*cot(y)*h(x,y,dy=1) - h(x,y,dy=2) - (-1 + x)*(3*h(x,y,dx=1) + (-1 + x)*h(x,y,dx=2)))))/2)) - 2*(-1 + x)^2*g(x,y)^3*(4*(-1/2*((1 + x)*cot(y)*f(x,y,dy=1)) + (-1 + x)*x*f(x,y,dx=1))*(h(x,y,dy=1)^2 + (-1 + x)^2*h(x,y,dx=1)^2) + 2*(-1 + x)*h(x,y)^2*(-(cot(y)*f(x,y,dy=1)) - f(x,y,dy=2) + 2*(-1 + x)*(3*f(x,y,dx=1) + (-1 + x)*f(x,y,dx=2))) - h(x,y)*(-(h(x,y,dy=1)*((1 + x)*cot(y)*f(x,y,dy=2) - (-1 + x)*((-1 + 7*x)*cot(y)*f(x,y,dx=1) + 4*x*f(x,y,dx=1,dy=1) + (-1 + x^2)*cot(y)*f(x,y,dx=2)))) + (-1 + x)*(2*x*h(x,y,dy=2)*f(x,y,dx=1) - 2*x*f(x,y,dy=2)*h(x,y,dx=1) + (-1 + x)*(-2*h(x,y,dx=1)*((1 + x)*cot(y)*f(x,y,dx=1,dy=1) - (-1 + x)*x*f(x,y,dx=2)) + f(x,y,dx=1)*((-7 + 15*x)*h(x,y,dx=1) + 2*(-1 + x)*x*h(x,y,dx=2)))) + f(x,y,dy=1)*((-1 + 5*x)*h(x,y,dy=1) - cot(y)*((1 + x)*h(x,y,dy=2) + (-1 + x)*((-1 + 7*x)*h(x,y,dx=1) + (-1 + x^2)*h(x,y,dx=2))))))) + 2*(-1 + x)^2*(1 + x)*f(x,y)*g(x,y)^2*h(x,y)*((3*rh^2*(1 + x)^3*h(x,y)*(f(x,y,dy=1)^2 + (-1 + x)^2*f(x,y,dx=1)^2)*(f(x,y,dy=1)*g(x,y,dy=1) + (-1 + x)^2*f(x,y,dx=1)*g(x,y,dx=1)))/2 - ((1 + x)*g(x,y)*h(x,y)*(2*(1 + x)*(rh - rh*x)^2*f(x,y,dx=1)*(f(x,y,dy=1)^2 + (-1 + x)^2*f(x,y,dx=1)^2) - rh*(1 - x)*sin(y)^2*(4*rh*(1 - x)*W(x,y)^2*(7*f(x,y,dy=1)*g(x,y,dy=1) - 3*(-1 + x)^2*f(x,y,dx=1)*g(x,y,dx=1)) - rh*(1 - x)*(-1 + x)^2*f(x,y,dx=1)*(-7*W(x,y,dy=1)^2*g(x,y,dx=1) + 10*g(x,y,dy=1)*W(x,y,dy=1)*W(x,y,dx=1) + 3*(-1 + x)^2*g(x,y,dx=1)*W(x,y,dx=1)^2) + 4*(-1 + x)^2*W(x,y)*(rh*g(x,y,dy=1)*(5*W(x,y,dy=1)*f(x,y,dx=1) - 7*f(x,y,dy=1)*W(x,y,dx=1)) + rh*g(x,y,dx=1)*(5*f(x,y,dy=1)*W(x,y,dy=1) + 3*(-1 + x)^2*f(x,y,dx=1)*W(x,y,dx=1))) - rh*(1 - x)*f(x,y,dy=1)*(10*(-1 + x)^2*W(x,y,dy=1)*g(x,y,dx=1)*W(x,y,dx=1) - g(x,y,dy=1)*(-3*W(x,y,dy=1)^2 + 7*(-1 + x)^2*W(x,y,dx=1)^2))) - 6*rh^2*(1 + x)^2*(f(x,y,dy=1)^2 + (-1 + x)^2*f(x,y,dx=1)^2)*(cot(y)*f(x,y,dy=1) + f(x,y,dy=2) + (-1 + x)^2*f(x,y,dx=2))))/2 - 2*g(x,y)^2*sin(y)*(2*rh^2*(-1 + x)^2*(1 + x)*sin(y)*(2*W(x,y)*f(x,y,dy=1) - (-1 + x)*(W(x,y,dy=1)*f(x,y,dx=1) - f(x,y,dy=1)*W(x,y,dx=1)))*(2*W(x,y)*h(x,y,dy=1) - (-1 + x)*(W(x,y,dy=1)*h(x,y,dx=1) - h(x,y,dy=1)*W(x,y,dx=1))) + (h(x,y)*(-(rh^2*(-1 + x)^4*sin(y)*(-5*W(x,y,dy=1)^2*f(x,y,dx=1) + 6*f(x,y,dy=1)*W(x,y,dy=1)*W(x,y,dx=1) + (-1 + x)^2*f(x,y,dx=1)*W(x,y,dx=1)^2)) - 4*rh^2*W(x,y)^2*((-1 + x)^4*sin(y)*f(x,y,dx=1) + (-1 + x)^2*(1 + x)*(15*cos(y)*f(x,y,dy=1) + sin(y)*(3*f(x,y,dy=2) - (-1 + x)*(2*f(x,y,dx=1) + (-1 + x)*f(x,y,dx=2))))) + (1 + x)*(2*rh^2*(-1 + x)^4*W(x,y,dy=1)*(4*sin(y)*W(x,y,dx=1)*f(x,y,dx=1,dy=1) + f(x,y,dx=1)*(9*cos(y)*W(x,y,dx=1) - 2*sin(y)*W(x,y,dx=1,dy=1))) - rh^2*(-1 + x)^2*sin(y)*W(x,y,dy=1)^2*(-f(x,y,dy=2) + (-1 + x)*(4*f(x,y,dx=1) + 3*(-1 + x)*f(x,y,dx=2))) + (rh - rh*x)^2*sin(y)*W(x,y,dx=1)*(6*(-1 + x)^2*W(x,y,dy=2)*f(x,y,dx=1) - 3*(-1 + x)^2*f(x,y,dy=2)*W(x,y,dx=1) - (1 - x)^3*W(x,y,dx=1)*(2*f(x,y,dx=1) + (-1 + x)*f(x,y,dx=2)) - 2*(1 - x)^3*f(x,y,dx=1)*(3*W(x,y,dx=1) + (-1 + x)*W(x,y,dx=2))) + rh*f(x,y,dy=1)*(3*rh*(-1 + x)^2*cos(y)*W(x,y,dy=1)^2 - rh*(-1 + x)^4*W(x,y,dx=1)*(15*cos(y)*W(x,y,dx=1) + 4*sin(y)*W(x,y,dx=1,dy=1)) + 2*sin(y)*W(x,y,dy=1)*(rh*(-1 + x)^2*W(x,y,dy=2) - 3*rh*(1 - x)^3*(2*W(x,y,dx=1) + (-1 + x)*W(x,y,dx=2))))) + 4*rh*W(x,y)*((1 - x)^3*sin(y)*(3*rh*f(x,y,dy=1)*W(x,y,dy=1) + rh*(-1 + x)^2*f(x,y,dx=1)*W(x,y,dx=1)) - rh*(1 + x)*(-((-1 + x)^3*W(x,y,dy=1)*(9*cos(y)*f(x,y,dx=1) + 4*sin(y)*f(x,y,dx=1,dy=1))) + (-1 + x)^2*f(x,y,dy=1)*(3*sin(y)*W(x,y,dy=1) + (-1 + x)*(15*cos(y)*W(x,y,dx=1) + 2*sin(y)*W(x,y,dx=1,dy=1))) + (1 - x)*sin(y)*(3*(-1 + x)^2*W(x,y,dy=2)*f(x,y,dx=1) - 3*(-1 + x)^2*f(x,y,dy=2)*W(x,y,dx=1) - (1 - x)^3*W(x,y,dx=1)*(2*f(x,y,dx=1) + (-1 + x)*f(x,y,dx=2)) - (1 - x)^3*f(x,y,dx=1)*(3*W(x,y,dx=1) + (-1 + x)*W(x,y,dx=2)))))))/2)) + f(x,y)^2*g(x,y)*(3*rh^2*(-1 + x)^2*(1 + x)^4*h(x,y)^2*(2*(-1 + x)^2*f(x,y,dy=1)*g(x,y,dy=1)*f(x,y,dx=1)*g(x,y,dx=1) + f(x,y,dy=1)^2*(2*g(x,y,dy=1)^2 + (-1 + x)^2*g(x,y,dx=1)^2) + (-1 + x)^2*f(x,y,dx=1)^2*(g(x,y,dy=1)^2 + 2*(-1 + x)^2*g(x,y,dx=1)^2)) + rh^2*(1 + x)^2*g(x,y)*h(x,y)*((-1 + x^2)^2*(2*(-1 + x)^2*f(x,y,dy=1)*f(x,y,dx=1)*(h(x,y,dy=1)*g(x,y,dx=1) + g(x,y,dy=1)*h(x,y,dx=1)) + f(x,y,dy=1)^2*(5*g(x,y,dy=1)*h(x,y,dy=1) + 3*(-1 + x)^2*g(x,y,dx=1)*h(x,y,dx=1)) + (-1 + x)^2*f(x,y,dx=1)^2*(3*g(x,y,dy=1)*h(x,y,dy=1) + 5*(-1 + x)^2*g(x,y,dx=1)*h(x,y,dx=1))) - h(x,y)*(6*(-1 + x)^4*(1 + x)*(f(x,y,dy=1)^2 + (-1 + x)^2*f(x,y,dx=1)^2)*g(x,y,dx=1) - (1 - x)^3*sin(y)^2*((1 - x)*W(x,y)^2*(8*g(x,y,dy=1)^2 - 4*(-1 + x)^2*g(x,y,dx=1)^2) + 6*(-1 + x)^3*g(x,y,dy=1)*W(x,y,dy=1)*g(x,y,dx=1)*W(x,y,dx=1) + 4*(-1 + x)^2*W(x,y)*(3*g(x,y,dy=1)*W(x,y,dy=1)*g(x,y,dx=1) - 2*g(x,y,dy=1)^2*W(x,y,dx=1) + (-1 + x)^2*g(x,y,dx=1)^2*W(x,y,dx=1)) - (1 - x)*g(x,y,dy=1)^2*(W(x,y,dy=1)^2 - 2*(-1 + x)^2*W(x,y,dx=1)^2) + (1 - x)^3*g(x,y,dx=1)^2*(2*W(x,y,dy=1)^2 - (-1 + x)^2*W(x,y,dx=1)^2)) + 2*(-1 + x^2)^2*(2*f(x,y,dy=1)*(g(x,y,dy=1)*(2*f(x,y,dy=2) + (-1 + x)*f(x,y,dx=1)) + (-1 + x)^2*(2*g(x,y,dx=1)*f(x,y,dx=1,dy=1) + f(x,y,dx=1)*g(x,y,dx=1,dy=1))) + f(x,y,dx=1)*((-1 + x)^2*g(x,y,dy=2)*f(x,y,dx=1) + (-1 + x)^2*g(x,y,dy=1)*(cot(y)*f(x,y,dx=1) + 4*f(x,y,dx=1,dy=1)) - 4*(1 - x)^3*g(x,y,dx=1)*(2*f(x,y,dx=1) + (-1 + x)*f(x,y,dx=2)) - 2*(1 - x)^3*f(x,y,dx=1)*(g(x,y,dx=1) + (-1 + x)*g(x,y,dx=2))) + f(x,y,dy=1)^2*(cot(y)*g(x,y,dy=1) + 2*g(x,y,dy=2) + (-1 + x)*(3*g(x,y,dx=1) + (-1 + x)*g(x,y,dx=2)))))) + 2*(1 + x)*g(x,y)^2*(2*(1 + x)^3*(rh - rh*x)^2*(f(x,y,dy=1)^2 + (-1 + x)^2*f(x,y,dx=1)^2)*(h(x,y,dy=1)^2 + (-1 + x)^2*h(x,y,dx=1)^2) + ((1 + x)*h(x,y)*(2*(-1 + x)^2*(1 + x)*(rh - rh*x)^2*(2*f(x,y,dy=1)*h(x,y,dy=1)*f(x,y,dx=1) - 5*f(x,y,dy=1)^2*h(x,y,dx=1) - 3*(-1 + x)^2*f(x,y,dx=1)^2*h(x,y,dx=1)) + rh*(1 - x)^3*sin(y)^2*(4*rh*(1 - x)*W(x,y)^2*(5*g(x,y,dy=1)*h(x,y,dy=1) - (-1 + x)^2*g(x,y,dx=1)*h(x,y,dx=1)) - rh*(1 - x)*(-1 + x)^2*g(x,y,dx=1)*(-5*W(x,y,dy=1)^2*h(x,y,dx=1) + 6*h(x,y,dy=1)*W(x,y,dy=1)*W(x,y,dx=1) + (-1 + x)^2*h(x,y,dx=1)*W(x,y,dx=1)^2) + 4*(-1 + x)^2*W(x,y)*(rh*h(x,y,dy=1)*(3*W(x,y,dy=1)*g(x,y,dx=1) - 5*g(x,y,dy=1)*W(x,y,dx=1)) + rh*h(x,y,dx=1)*(3*g(x,y,dy=1)*W(x,y,dy=1) + (-1 + x)^2*g(x,y,dx=1)*W(x,y,dx=1))) - rh*(1 - x)*g(x,y,dy=1)*(6*(-1 + x)^2*W(x,y,dy=1)*h(x,y,dx=1)*W(x,y,dx=1) + h(x,y,dy=1)*(W(x,y,dy=1)^2 - 5*(-1 + x)^2*W(x,y,dx=1)^2))) + 2*rh^2*(-1 + x)^2*(1 + x)^2*(2*f(x,y,dy=1)*((-1 + x)^2*h(x,y,dx=1)*(cot(y)*f(x,y,dx=1) - 2*f(x,y,dx=1,dy=1)) - h(x,y,dy=1)*(f(x,y,dy=2) - (-1 + x)^2*f(x,y,dx=2))) - (-1 + x)^2*f(x,y,dx=1)*(h(x,y,dy=2)*f(x,y,dx=1) - 2*f(x,y,dy=2)*h(x,y,dx=1) - 8*f(x,y,dx=1)*h(x,y,dx=1) + 8*x*f(x,y,dx=1)*h(x,y,dx=1) + h(x,y,dy=1)*(cot(y)*f(x,y,dx=1) + 4*f(x,y,dx=1,dy=1)) + 2*h(x,y,dx=1)*f(x,y,dx=2) - 4*x*h(x,y,dx=1)*f(x,y,dx=2) + 2*x^2*h(x,y,dx=1)*f(x,y,dx=2) + f(x,y,dx=1)*h(x,y,dx=2) - 2*x*f(x,y,dx=1)*h(x,y,dx=2) + x^2*f(x,y,dx=1)*h(x,y,dx=2)) + f(x,y,dy=1)^2*(cot(y)*h(x,y,dy=1) - h(x,y,dy=2) - (-1 + x)*(4*h(x,y,dx=1) + (-1 + x)*h(x,y,dx=2))))))/2 - 2*h(x,y)^2*((rh^2*(-1 + x)^5*sin(y)^2*(4*(-1 + x)*W(x,y)^2*g(x,y,dx=1) + 4*W(x,y)*(2*g(x,y,dy=1)*W(x,y,dy=1) + (-1 + x)^2*g(x,y,dx=1)*W(x,y,dx=1)) + (-1 + x)*(-3*W(x,y,dy=1)^2*g(x,y,dx=1) + 4*g(x,y,dy=1)*W(x,y,dy=1)*W(x,y,dx=1) + (-1 + x)^2*g(x,y,dx=1)*W(x,y,dx=1)^2)))/2 + (rh - rh*x^2)^2*(-((-1 + x)*f(x,y,dy=1)^2) + 2*(1 - x)*(2*f(x,y,dy=1)^2 + (-1 + x)^2*f(x,y,dx=1)^2) + 2*(-1 + x)^2*f(x,y,dy=1)*(cot(y)*f(x,y,dx=1) - 2*f(x,y,dx=1,dy=1)) + 2*(-1 + x)^2*f(x,y,dx=1)*(f(x,y,dy=2) - (-1 + x)*(2*f(x,y,dx=1) + (-1 + x)*f(x,y,dx=2)))) - 2*rh^2*(1 + x)^3*((-1 + x)^2*f(x,y,dy=1)^2 - (-1 + x)^4*f(x,y,dx=1)^2 + (-1 + x)^4*f(x,y,dx=1,dy=1)^2 + (1 - x)^3*f(x,y,dy=2)*(f(x,y,dx=1) + (-1 + x)*f(x,y,dx=2)) + 2*(-1 + x)^4*f(x,y,dx=1)*(2*f(x,y,dx=1) + (-1 + x)*f(x,y,dx=2)) - (-1 + x)^2*f(x,y,dy=1)*(cot(y)*f(x,y,dy=2) + (-1 + x)*(cot(y)*f(x,y,dx=1) - 2*f(x,y,dx=1,dy=1) + (-1 + x)*cot(y)*f(x,y,dx=2)))) - (rh^2*(-1 + x)^2*(1 + x)*sin(y)*(-36*(-1 + x)^2*cos(y)*W(x,y)^2*g(x,y,dy=1) + 3*(-1 + x)^2*cos(y)*g(x,y,dy=1)*W(x,y,dy=1)^2 + (-1 + x)^2*sin(y)*W(x,y,dy=1)^2*g(x,y,dy=2) - 2*(-1 + x)^3*sin(y)*W(x,y,dy=1)^2*g(x,y,dx=1) + 12*(-1 + x)^4*cos(y)*W(x,y,dy=1)*g(x,y,dx=1)*W(x,y,dx=1) + 4*(-1 + x)^4*sin(y)*W(x,y,dy=2)*g(x,y,dx=1)*W(x,y,dx=1) - 2*(-1 + x)^4*sin(y)*g(x,y,dy=2)*W(x,y,dx=1)^2 + 9*(-1 + x)^5*sin(y)*g(x,y,dx=1)*W(x,y,dx=1)^2 + 6*(-1 + x)^4*sin(y)*W(x,y,dy=1)*W(x,y,dx=1)*g(x,y,dx=1,dy=1) + 12*(-1 + x)^3*W(x,y)*W(x,y,dy=1)*(2*cos(y)*g(x,y,dx=1) + sin(y)*g(x,y,dx=1,dy=1)) - 2*(-1 + x)^4*sin(y)*W(x,y,dy=1)*g(x,y,dx=1)*W(x,y,dx=1,dy=1) - (-1 + x)^4*g(x,y,dy=1)*W(x,y,dx=1)*(9*cos(y)*W(x,y,dx=1) + 2*sin(y)*W(x,y,dx=1,dy=1)) - 4*(-1 + x)^2*W(x,y)*g(x,y,dy=1)*(sin(y)*W(x,y,dy=1) + (-1 + x)*(9*cos(y)*W(x,y,dx=1) + sin(y)*W(x,y,dx=1,dy=1))) - 2*(-1 + x)^4*sin(y)*W(x,y,dy=1)^2*g(x,y,dx=2) + (-1 + x)^6*sin(y)*W(x,y,dx=1)^2*g(x,y,dx=2) - 4*(-1 + x)^2*sin(y)*W(x,y)^2*(2*g(x,y,dy=2) - (-1 + x)*(3*g(x,y,dx=1) + (-1 + x)*g(x,y,dx=2))) + 2*(-1 + x)^6*sin(y)*g(x,y,dx=1)*W(x,y,dx=1)*W(x,y,dx=2) - 4*(1 - x)*sin(y)*W(x,y)*(2*(-1 + x)^2*W(x,y,dy=2)*g(x,y,dx=1) - 2*(-1 + x)^2*g(x,y,dy=2)*W(x,y,dx=1) - (1 - x)^3*W(x,y,dx=1)*(2*g(x,y,dx=1) + (-1 + x)*g(x,y,dx=2)) - (1 - x)^3*g(x,y,dx=1)*(4*W(x,y,dx=1) + (-1 + x)*W(x,y,dx=2))) + 2*(-1 + x)^2*sin(y)*g(x,y,dy=1)*W(x,y,dy=1)*(W(x,y,dy=2) + (-1 + x)*(5*W(x,y,dx=1) + 2*(-1 + x)*W(x,y,dx=2)))))/2)) - 2*(-1 + x)^2*g(x,y)^3*(2*rh^2*(-1 + x)^2*(1 + x)^2*sin(y)^2*(h(x,y,dy=1)^2 + (-1 + x)^2*h(x,y,dx=1)^2)*(4*W(x,y)^2 + W(x,y,dy=1)^2 + 4*(-1 + x)*W(x,y)*W(x,y,dx=1) + (-1 + x)^2*W(x,y,dx=1)^2) + (1 + x)*h(x,y)*sin(y)*(rh^2*(-1 + x)^4*sin(y)*(3*W(x,y,dy=1)^2*h(x,y,dx=1) - 2*h(x,y,dy=1)*W(x,y,dy=1)*W(x,y,dx=1) + (-1 + x)^2*h(x,y,dx=1)*W(x,y,dx=1)^2) + 4*rh^2*W(x,y)^2*((-1 + x)^4*sin(y)*h(x,y,dx=1) - (-1 + x)^2*(1 + x)*(9*cos(y)*h(x,y,dy=1) + sin(y)*(h(x,y,dy=2) + (-1 + x)*(6*h(x,y,dx=1) + (-1 + x)*h(x,y,dx=2))))) + (1 + x)*(2*rh^2*(-1 + x)^4*W(x,y,dy=1)*h(x,y,dx=1)*(3*cos(y)*W(x,y,dx=1) - 2*sin(y)*W(x,y,dx=1,dy=1)) - rh^2*(-1 + x)^2*sin(y)*W(x,y,dy=1)^2*(h(x,y,dy=2) + (-1 + x)*(4*h(x,y,dx=1) + (-1 + x)*h(x,y,dx=2))) + 2*(rh - rh*x)^2*sin(y)*W(x,y,dx=1)*((-1 + x)^2*W(x,y,dy=2)*h(x,y,dx=1) - ((-1 + x)^2*h(x,y,dy=2)*W(x,y,dx=1))/2 - ((-1 + x)^3*W(x,y,dx=1)*(2*h(x,y,dx=1) + (-1 + x)*h(x,y,dx=2)))/2 + (1 - x)^3*h(x,y,dx=1)*(5*W(x,y,dx=1) + (-1 + x)*W(x,y,dx=2))) - rh^2*(-1 + x)^2*h(x,y,dy=1)*(3*cos(y)*W(x,y,dy=1)^2 + (-1 + x)^2*W(x,y,dx=1)*(9*cos(y)*W(x,y,dx=1) + 4*sin(y)*W(x,y,dx=1,dy=1)) + 2*sin(y)*W(x,y,dy=1)*(W(x,y,dy=2) + (-1 + x)*(2*W(x,y,dx=1) - (-1 + x)*W(x,y,dx=2))))) + 4*rh*W(x,y)*((1 - x)^3*sin(y)*(rh*h(x,y,dy=1)*W(x,y,dy=1) - rh*(-1 + x)^2*h(x,y,dx=1)*W(x,y,dx=1)) - rh*(1 + x)*(-3*(-1 + x)^3*cos(y)*W(x,y,dy=1)*h(x,y,dx=1) + (-1 + x)^2*h(x,y,dy=1)*(5*sin(y)*W(x,y,dy=1) + (-1 + x)*(9*cos(y)*W(x,y,dx=1) + 2*sin(y)*W(x,y,dx=1,dy=1))) + (1 - x)*sin(y)*((-1 + x)^2*W(x,y,dy=2)*h(x,y,dx=1) - (-1 + x)^2*h(x,y,dy=2)*W(x,y,dx=1) + (1 - x)^3*W(x,y,dx=1)*(2*h(x,y,dx=1) + (-1 + x)*h(x,y,dx=2)) + (1 - x)^3*h(x,y,dx=1)*(7*W(x,y,dx=1) + (-1 + x)*W(x,y,dx=2)))))) + 2*rh*h(x,y)^2*(2*rh*(-1 + x)^2*(1 + x)*(7 + 5*x + (5 + 7*x)*cos(2*y))*W(x,y)^2 + 2*rh*(-1 + x)^4*sin(y)^2*W(x,y,dy=1)^2 - rh*(-1 + x)^3*(1 + x)*sin(y)*(2*sin(y)*W(x,y,dy=1)^2 + (-1 + x)*sin(y)*W(x,y,dx=1)*(-2*W(x,y,dy=2) + (-1 + x)*W(x,y,dx=1)) - 2*(-1 + x)*W(x,y,dy=1)*(3*cos(y)*W(x,y,dx=1) - sin(y)*W(x,y,dx=1,dy=1))) + 2*rh*(-1 + x^2)^2*(3*sin(y)^2*W(x,y,dy=1)^2 - (-1 + x)*sin(y)^2*W(x,y,dy=2)*(3*W(x,y,dx=1) + (-1 + x)*W(x,y,dx=2)) - 3*(-1 + x)*sin(y)*W(x,y,dy=1)*(cos(y)*W(x,y,dx=1) - sin(y)*W(x,y,dx=1,dy=1) + (-1 + x)*cos(y)*W(x,y,dx=2)) + (-1 + x)^2*(3*W(x,y,dx=1)^2 + sin(y)^2*W(x,y,dx=1,dy=1)^2 + sin(y)*W(x,y,dx=1)*(3*cos(y)*W(x,y,dx=1,dy=1) + (-1 + x)*sin(y)*W(x,y,dx=2)))) - 4*(1 + x)*W(x,y)*((1 - x)^3*sin(y)*(3*rh*cos(y)*W(x,y,dy=1) + rh*sin(y)*(W(x,y,dy=2) - (-1 + x)*W(x,y,dx=1))) - (rh*(-1 + x)^2*(1 + x)*(6*sin(2*y)*W(x,y,dy=1) + (-1 + x)*(3*(3 + cos(2*y))*W(x,y,dx=1) + 2*sin(y)*(3*cos(y)*W(x,y,dx=1,dy=1) + (-1 + x)*sin(y)*W(x,y,dx=2)))))/2))))))/(32*rh^6*(1 + x)^4*f(x,y)^2*g(x,y)^5*h(x,y)^4)
end

LoadSystem()

export Field, LoadSystem, PrintData, interpolate1D, interpolate1D!, interpolate, interpolate!, GetMass, GetJ, GetTh, GetAh, GetωχKerr, GetωχKerrN, quantities_kerr, GetLe, GetLp, Sphericity, vH, get_quantities, fKerrN, gKerrN, hKerrN, WKerrN, AtKerrN, AφKerrN, Print_Petrov, gtt, gtphi, gphiphi, grr, gthetatheta, Ergosphere, CircumferencialRadius, LightRing, ISCO, Ricci_Horizon, GB_Horizon, RicciScalar, GBScalar, X, Mx, Y, My, Nx, Ny

end