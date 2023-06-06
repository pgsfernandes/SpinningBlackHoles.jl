function V(p::Float64,dp::Int64=0)
    μ=0.0
    λ=0.0
    if dp==0
        return 0.5*μ^2*p^2 + λ/4 * p^4
    elseif dp==1
        return μ^2*p + λ*p^3
    elseif dp==2
        return μ^2 + 3*λ*p^2
    else
        return 6*λ*p
    end
end

function ϕpert(x::Float64,y::Float64,k::Float64,rh::Float64,χ::Float64,q::Float64)
    #RETURNS A PERTURBATIVE SOLUTION FOR THE SCALAR FIELD TO ORDER g_aγγ AND χ
    sq=sqrt(1-q^2-χ^2)
    return (2*k*rh*χ*cos(y)*(sq/(rh+sqrt(1-q^2)*rh)+2/(rh*(1-x)*(-1-4/(-1+x)^2+4/((-1+x)*sq)))+((5-2*x+x^2)*(-1+q^2+χ^2)*log((4*sqrt(1-q^2)+5*sq+x^2*sq-2*x*(2*sqrt(1-q^2)+sq))^2/(4+5*sq+x^2*sq-2*x*(2+sq))^2)/(8*q^2*rh*(-1+x)))))/sq

    #return 0.0
end

function GetMs(rh::Float64,p::Field,f::Field,g::Field,h::Field)
    function to_integrate(v::Vector{Float64},p::Field,f::Field,g::Field,h::Field)
        x=v[1]
        y=v[2]
        return 4*rh^3*(1+x)/(1-x)^4 * g(x,y)^(3/2) * h(x,y)/f(x,y)*sin(y) *( V(p(x,y)) )
    end

    integral, err = hcubature(v->to_integrate(v,p,f,g,h), [-1.0, 0.0], [1.0-1e-9, pi], abstol=1e-10)
    return integral
end

#k is g_aγγ
function OneSolution(WBC::Float64, rh::Float64, q::Float64, k::Float64, spin::Bool=true, tol::Float64=1e-11; guess=nothing, branch::Int64=1, ToPrint::Bool=false, ergosphere::Bool=false, light_ring::Bool=false, isco::Bool=false, petrov::Bool=false, sphericity::Bool=false, linvel::Bool=false)
    println()
    println("k=",k)
    println("q=",q)

    if spin
        println("χ=",WBC)
    else
        println("Ωh=",WBC)
    end

    #DECLARE OUR FIELDS
    f=Field();
    g=Field();
    h=Field();
    W=Field();
    p=Field(); p.type=2;
    A=Field();
    B=Field();

    if guess==nothing

        ωKN, χKN = GetωχKerrN(WBC,q,spin,branch)
        
        interpolate!(f,(x,y)->fKerrN(x,y,rh,ωKN,χKN,q,spin,branch),Nx,Ny)
        interpolate!(g,(x,y)->gKerrN(x,y,rh,ωKN,χKN,q,spin,branch),Nx,Ny)
        interpolate!(h,(x,y)->hKerrN(x,y,rh,ωKN,χKN,q,spin,branch),Nx,Ny)
        interpolate!(W,(x,y)->WKerrN(x,y,rh,ωKN,χKN,q,spin,branch),Nx,Ny)
        M, J, χ, Th, Ah, Ωh = get_quantities(f,g,h,W,rh)

        interpolate!(p,(x,y)->ϕpert(x,y,k,rh,χ,q),Nx,Ny)
        interpolate!(A,(x,y)->AtKerrN(x,y,rh,ωKN,χKN,q,spin,branch),Nx,Ny)
        interpolate!(B,(x,y)->AφKerrN(x,y,rh,ωKN,χKN,q,spin,branch),Nx,Ny)
    else
        #USE guess AS INITIAL GUESS IF THE USER PROVIDES ONE
        f.a .= guess[0*Nx+1:1*Nx,:]
        g.a .= guess[1*Nx+1:2*Nx,:]
        h.a .= guess[2*Nx+1:3*Nx,:]
        W.a .= guess[3*Nx+1:4*Nx,:]
        p.a .= guess[4*Nx+1:5*Nx,:]
        A.a .= guess[5*Nx+1:6*Nx,:]
        B.a .= guess[6*Nx+1:7*Nx,:]
    end

    #HERE WE SOLVE THE SYSTEM OF PDEs
    sol=solve_system(f,g,h,W,p,A,B,tol,WBC,rh,q,k,spin,true,15)

    convergence=sol.f_converged || sol.x_converged

    if convergence
        #OBTAIN PHYSICAL QUANTITIES OF OUR SOLUTION
        M, J, χ, Th, Ah, Ωh = get_quantities(f,g,h,W,rh)

        Q = 2*rh*A(Nx,1,dx=1)
        μM = -2*rh*B(Nx,1,dx=1)
        Φ = A(Nx,1)
        gg = 2*M*μM/(Q*J)
        println("M*Ωh=",M*Ωh)
        println("Φ=",Φ)
        println("g=",gg)
        println("Th*8πM=",Th*8*π*M)
        println("Ah/16πM^2=",Ah/(16*π*M^2))
        println("1-χ^2-q^2=", 1.0-χ^2-q^2)
        Ms=GetMs(rh,p,f,g,h)
		smarr = abs(1-(0.5*Th*Ah+2*Ωh*J+Φ*Q-Ms)/M)
        println("smarr=",smarr)

        if ergosphere
            Ergosphere(f,g,h,W,rh,"ergosphere.dat",100,q=q)
        end

        if petrov
            Print_Petrov(f,g,h,W,rh)
        end

        if light_ring
            LightRing(f,g,h,W,rh,q)
        end

        if isco
            ISCO(f,g,h,W,rh,q)
        end

        if sphericity
            println("Sphericity = ",Sphericity(f,g,h,rh))
        end

        if linvel
            println("Linear velocity of the horizon = ",vH(f,g,h,rh,Ωh))
        end

        if ToPrint
            PrintData("func_f.dat",f.a)
            PrintData("func_g.dat",g.a)
            PrintData("func_h.dat",h.a)
            PrintData("func_W.dat",W.a)
            PrintData("func_p.dat",p.a)
            PrintData("func_A.dat",A.a)
            PrintData("func_B.dat",B.a)
            PrintData("par.dat",[χ,rh,q,k])
        end
    end
end

OneSolution(0.001,1.0,0.2,0.5,true,1e-8,ToPrint=false,ergosphere=false,petrov=false,light_ring=false,isco=false,sphericity=false,linvel=false)