#DEFINES THE GAUSS-BONNET COUPLING FUNCTION AND ITS DERIVATIVES WITH RESPECT TO ϕ (HERE DENOTED BY p)
function F(p::Float64,dp::Int64=0)
    γ=1.0
    if dp==0
        #return p
        return exp(γ*p)
    elseif dp==1
        #return 1.0
        return γ*exp(γ*p)
    elseif dp==2
        #return 0.0
        return γ^2*exp(γ*p)
    else
        #return 0.0
        return γ^3*exp(γ*p)
    end
end

function GetQs(rh::Float64,p::Field)
    #RETURNS THE SCALAR CHARGE OF THE SOLUTION
    return -2*rh*p(Nx,1,dx=1)
end

function GetMs(rh::Float64,g::Field,p::Field)
    #SHOULD WORK FOR ALL COUPLINGS EXCEPT THE EXPONENTIAL F(p) = exp(γ p), IN WHICH CASE Ms = Qs/γ
    function to_integrate(v::Vector{Float64},g::Field,p::Field)
        x=v[1]
        y=v[2]
        if x>1.0-1e-12
            return (1.0 - F(p(x,y))*F(p(x,y),2)/F(p(x,y),1)^2) * (1+x)*sqrt(g(x,y))*sin(y)*(p(x,y,dx=1)^2)
        else
            return (1.0 - F(p(x,y))*F(p(x,y),2)/F(p(x,y),1)^2) * (1+x)*sqrt(g(x,y))*sin(y)*(p(x,y,dy=1)^2/(1-x)^2 + p(x,y,dx=1)^2)
        end
    end

    integral, err = pcubature(v->to_integrate(v,g,p), [-1.0, 0.0], [1.0, pi], abstol=1e-10)
    return rh/2.0 * integral
end

function Ricci_H(rh::Float64,f::Field,g::Field,h::Field)
    #RETURNS r^2 siny * R_H, WHERE R_H IS THE INDUCED RICCI SCALAR ON THE HORIZON
    x=-1.0
    y=pi/2
    return (-2*g(x,y)^2*h(x,y)*sin(y)*f(x,y,dy=1)^2 + f(x,y)*g(x,y)^2*(f(x,y,dy=1)*(2*cos(y)*h(x,y) - sin(y)*h(x,y,dy=1)) + 2*h(x,y)*sin(y)*f(x,y,dy=2)) + f(x,y)^2*(2*cos(y)*g(x,y)*(-(h(x,y)*g(x,y,dy=1)) + g(x,y)*h(x,y,dy=1)) + sin(y)*(4*g(x,y)^2*h(x,y) + 2*h(x,y)*g(x,y,dy=1)^2 + g(x,y)*(g(x,y,dy=1)*h(x,y,dy=1) - 2*h(x,y)*g(x,y,dy=2)))))/(2*f(x,y)*g(x,y)^3*h(x,y)^2)
end

function Entropy(rh::Float64,A::Float64,f::Field,g::Field,h::Field,p::Field)
    #RETURNS THE ENTROPY OF A GB BLACK HOLE
    x=-1.0

    #Actually r^2 siny * R_H
    R_H = y->(-2*g(x,y)^2*h(x,y)*sin(y)*f(x,y,dy=1)^2 + f(x,y)*g(x,y)^2*(f(x,y,dy=1)*(2*cos(y)*h(x,y) - sin(y)*h(x,y,dy=1)) + 2*h(x,y)*sin(y)*f(x,y,dy=2)) + f(x,y)^2*(2*cos(y)*g(x,y)*(-(h(x,y)*g(x,y,dy=1)) + g(x,y)*h(x,y,dy=1)) + sin(y)*(4*g(x,y)^2*h(x,y) + 2*h(x,y)*g(x,y,dy=1)^2 + g(x,y)*(g(x,y,dy=1)*h(x,y,dy=1) - 2*h(x,y)*g(x,y,dy=2)))))/(2*f(x,y)*g(x,y)^3*h(x,y)^2)

    integral, err = pquadrature(y -> g(x,y)/f(x,y)*sqrt(h(x,y))*F(p(x,y))*R_H(y), 0.0, pi, abstol=1e-14)

    return GetAh(f,g,h,rh)/4.0 + A*rh^2 * pi/4.0 * integral
end

function ϕpert(x::Float64,A::Float64)
    #RETURNS A PERTURBATIVE SOLUTION FOR THE SCALAR FIELD TO ORDER α/rh^2 = A
    return A * (415-1047 * x+942 * x^2-358 * x^3+51 * x^4-3 * x^5)/(12 * (-3+x)^6)
end

function OneSolution(;WBC::Float64, rh::Float64, A::Float64, spin::Bool=true, tol::Float64=1e-11, guess=nothing, branch::Int64=1, ToPrint::Bool=false, ergosphere::Bool=false, light_ring::Bool=false, isco::Bool=false, petrov::Bool=false, sphericity::Bool=false, linvel::Bool=false)
    #FUNCTION TO COMPUTE ONE GB BLACK HOLE SOLUTION.
    #WBC IS THE VALUE CONCERNING THE BOUNDARY CONDITION OF THE FUNCTION W: THE DIMENSIONLESS SPIN χ IF "spin" IS TRUE, OR THE ANGULAR VELOCITY OF THE HORIZON Ωh IF "spin" IS FALSE. rh IS THE HORIZON RADIUS, A = α/rh^2, tol DEFINES THE TOLERANCE TO DECLARE CONVERGENCE (NORM DIFFERENCE IN THE SPECTRAL COEFFICIENTS BETWEEN TWO CONSECUTIVE ITERATIONS), guess IS AN INITIAL GUESS THAT CAN BE IMPORTED FROM DATA FILES (OTHERWISE A COMPARABLE KERR BH WILL BE USED), branch CONCERNS THE TWO BRANCHES OF KERR SOLUTIONS (IMPORTANT IF Ωh IS USED AS BC FOR W), ToPrint PRINTS THE SOLUTIONS TO A .DAT FILE, ergosphere PRINTS THE REGION OF THE ERGOSPHERE IF SET TO true, light_ring COMPUTES THE LIGHT RINGS OF THE SOLUTION IF SET TO true (SIMILAR FOR ISCO), petrov PRINTS A FILE WITH THE VALUE OF THE LORENTZ INVARIANT SCALARS ACROSS THE SPACETIME TO DETERMINE THE PETROV TYPE, sphericity COMPUTES THE SPHERICITY OF THE SOLUTION, AND linvel COMPUTES THE LINEAR VELOCITY OF THE HORIZON.
    println()
    println("Solver initated with:")
    println("α/rh^2=",A)

    if spin
        println("χ=",WBC)
    else
        println("Ωh=",WBC)
    end
    println()

    #DECLARE OUR FIELDS
    f=Field();
    g=Field();
    h=Field();
    W=Field();
    p=Field();

    if guess==nothing
        #IF THERE IS NO GUESS WE USE A COMPARABLE KERR BLACK HOLE AS INITIAL GUESS. WE PERTURB IT WITH THE PERTURBATION CONTROLED BY δ. THIS CAN BE USEFUL WHEN TESTING STUFF 

        δ=1e-8*0.0
        function pert(x::Float64,y::Float64,δ::Float64)
            return 1-δ*(1+x)^2
        end
        
        ωKN, χKN = GetωχKerrN(WBC,0.0,spin,branch)
        
        interpolate!(f,(x,y)->fKerrN(x,y,rh,ωKN,χKN,0.0,spin,branch)*pert(x,y,δ),Nx,Ny)
        interpolate!(g,(x,y)->gKerrN(x,y,rh,ωKN,χKN,0.0,spin,branch)*pert(x,y,δ),Nx,Ny)
        interpolate!(h,(x,y)->hKerrN(x,y,rh,ωKN,χKN,0.0,spin,branch)*pert(x,y,δ),Nx,Ny)
        interpolate!(W,(x,y)->WKerrN(x,y,rh,ωKN,χKN,0.0,spin,branch)*pert(x,y,δ),Nx,Ny)
        
        interpolate1D!(p,x->ϕpert(x,A),Nx,Ny)

    else
        #USE guess AS INITIAL GUESS IF THE USER PROVIDES ONE
        f.a .= guess[0*Nx+1:1*Nx,:]
        g.a .= guess[1*Nx+1:2*Nx,:]
        h.a .= guess[2*Nx+1:3*Nx,:]
        W.a .= guess[3*Nx+1:4*Nx,:]
        p.a .= guess[4*Nx+1:5*Nx,:]
    end

    #THE LINE BELOW COMPUTES THE PHYSICAL QUANTITIES FOR A COMPARABLE KERR BH
    #Mkerr, Jkerr, χkerr, Thkerr, Ahkerr, Ωhkerr = quantities_kerr(WBC,rh,spin,branch)

    #HERE WE SOLVE THE SYSTEM OF PDEs
    sol=solve_system(f,g,h,W,p,tol,WBC,rh,A,spin,true,15)

    convergence=sol.f_converged || sol.x_converged

    if convergence
        println()
        println("Success! Solution converged!")
        #OBTAIN PHYSICAL QUANTITIES OF OUR SOLUTION
        M, J, χ, Th, Ah, Ωh = get_quantities(f,g,h,W,rh)

        println("χ=$χ")
        Qs = GetQs(rh,p)
        α = A*rh^2
        println("α/M^2=$(α/(M^2))")
        println("Qs/M=",Qs/M)

        #RELATION FOR THE LINEAR COUPLING
        #QsThRelation = abs(1-2*pi*α*Th/Qs)
        #println("|1-2pi*α*Th/Qs|=",QsThRelation)

        #SMARR RELATION FOR EXPONENTIAL COUPLING
        S=Entropy(rh,A,f,g,h,p)
        Smarr = abs(1-(2*Th*S+2*Ωh*J-Qs)/M)

        #SMARR RELATION GENERAL CASE
        #Ms=GetMs(rh,g,p)
        #Smarr = abs(1-(2*Th*S+2*Ωh*J-Ms)/M)

        println("Smarr = ", Smarr)
        
        if ergosphere
            Ergosphere(f,g,h,W,rh,"ergosphere.dat",100)
        end

        if petrov
            Print_Petrov(f,g,h,W,rh,detailed=false)
        end

        if light_ring
            LightRing(f,g,h,W,rh)
        end

        if isco
            ISCO(f,g,h,W,rh)
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
        end
    else
        println("Did not converge...")
    end
end

#OBTAIN ONE SOLUTION WITH χ=0.6
OneSolution(WBC=0.6,rh=1.0,A=1.0,spin=true,tol=1e-10,ToPrint=false,ergosphere=false,petrov=false,light_ring=false,isco=false,sphericity=false,linvel=false)

#OBTAIN ONE SOLUTION WITH Ωh=1/15
#OneSolution(WBC=1.0/15.0,rh=1.0,A=1.0,spin=false,tol=1e-10,branch=1,ToPrint=false,ergosphere=false,petrov=false,light_ring=false,isco=false,sphericity=false,linvel=false)

#EXAMPLE ON HOW TO OBTAIN ONE SOLUTION USING A guess AS INITTIAL GUESS
#OneSolution(WBC=1.0/15.0,rh=1.0,A=1.0,spin=false,tol=1e-10,guess=vcat(readdlm("func_f.dat"),readdlm("func_g.dat"),readdlm("func_h.dat"),readdlm("func_W.dat"),readdlm("func_p.dat")),branch=1,ToPrint=false,ergosphere=false,petrov=false,light_ring=false,isco=false,sphericity=false,linvel=false)

nothing