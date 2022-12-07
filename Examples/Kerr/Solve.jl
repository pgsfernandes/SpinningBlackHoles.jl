function OneSolution(WBC::Float64, rh::Float64, spin::Bool=true, tol::Float64=1e-11; guess=nothing, branch::Int64=1, ToPrint::Bool=false, ergosphere::Bool=false, light_ring::Bool=false, isco::Bool=false, petrov::Bool=false, sphericity::Bool=false, linvel::Bool=false)
    #FUNCTION TO COMPUTE ONE KERR BLACK HOLE SOLUTION.
    #WBC IS THE VALUE CONCERNING THE BOUNDARY CONDITION OF THE FUNCTION W: THE DIMENSIONLESS SPIN χ IF "spin" IS TRUE, OR THE ANGULAR VELOCITY OF THE HORIZON Ωh IF "spin" IS FALSE. rh IS THE HORIZON RADIUS, tol DEFINES THE TOLERANCE TO DECLARE CONVERGENCE (NORM DIFFERENCE IN THE SPECTRAL COEFFICIENTS BETWEEN TWO CONSECUTIVE ITERATIONS), guess IS AN INITIAL GUESS THAT CAN BE IMPORTED FROM DATA FILES (OTHERWISE A PERTURBED COMPARABLE KERR BH WILL BE USED), branch CONCERNTS THE TWO BRANCHES OF KERR SOLUTIONS (IMPORTANT IF Ωh IS USED AS BC FOR W), ToPrint PRINTS THE SOLUTIONS TO A .DAT FILE, ergosphere PRINTS THE REGION OF THE ERGOSPHERE IF SET TO true, light_ring COMPUTES THE LIGHT RINGS OF THE SOLUTION IF SET TO true (SIMILAR FOR ISCO), petrov PRINTS A FILE WITH THE VALUE OF THE LORENTZ INVARIANT SCALARS ACROSS THE SPACETIME TO DETERMINE THE PETROV TYPE, sphericity COMPUTES THE SPHERICITY OF THE SOLUTION, AND linvel COMPUTES THE LINEAR VELOCITY OF THE HORIZON.
    println()
    println("Solver initated with:")

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

    if guess==nothing
        #IF THERE IS NO GUESS WE USE A COMPARABLE PERTURBED KERR BLACK HOLE AS INITIAL GUESS. WE PERTURB IT WITH THE PERTURBATION CONTROLED BY δ. THIS CAN BE USEFUL WHEN TESTING STUFF 

        δ=1e-3*1.0
        function pert(x::Float64,y::Float64,δ::Float64)
            return 1-δ*(1+x)^2
        end
        
        ωKN, χKN = GetωχKerrN(WBC,0.0,spin,branch)
        
        interpolate!(f,(x,y)->fKerrN(x,y,rh,ωKN,χKN,0.0,spin,branch)*pert(x,y,δ),Nx,Ny)
        interpolate!(g,(x,y)->gKerrN(x,y,rh,ωKN,χKN,0.0,spin,branch)*pert(x,y,δ),Nx,Ny)
        interpolate!(h,(x,y)->hKerrN(x,y,rh,ωKN,χKN,0.0,spin,branch)*pert(x,y,δ),Nx,Ny)
        interpolate!(W,(x,y)->WKerrN(x,y,rh,ωKN,χKN,0.0,spin,branch)*pert(x,y,δ),Nx,Ny)
    else
        #USE guess AS INITIAL GUESS IF THE USER PROVIDES ONE
        f.a .= guess[0*Nx+1:1*Nx,:]
        g.a .= guess[1*Nx+1:2*Nx,:]
        h.a .= guess[2*Nx+1:3*Nx,:]
        W.a .= guess[3*Nx+1:4*Nx,:]
    end

    #THE LINE BELOW COMPUTES THE PHYSICAL QUANTITIES FOR A COMPARABLE KERR BH
    Mkerr, Jkerr, χkerr, Thkerr, Ahkerr, Ωhkerr = quantities_kerr(WBC,rh,spin,branch)

    #HERE WE SOLVE THE SYSTEM OF PDEs
    sol=solve_system(f,g,h,W,tol,WBC,rh,spin,true,15)

    convergence=sol.f_converged || sol.x_converged

    if convergence
        println()
        println("Success! Solution converged!")
        #OBTAIN PHYSICAL QUANTITIES OF OUR SOLUTION
        M, J, χ, Th, Ah, Ωh = get_quantities(f,g,h,W,rh)

        println("χ=$χ")

        println("Error in M = ",abs(1-Mkerr/M))
        println("Error in J = ",abs(1-Jkerr/J))
        println("Error in Th = ",abs(1-Thkerr/Th))
        println("Error in Ah = ",abs(1-Ahkerr/Ah))

        #SMARR RELATION
        S=Ah/4.0
        Smarr = abs(1-(2*Th*S+2*Ωh*J)/M)

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
#@time OneSolution(0.3,1.0,true,1e-10,ToPrint=false,ergosphere=false,petrov=false,light_ring=false,isco=false,sphericity=false,linvel=false)

#OBTAIN ONE SOLUTION WITH Ωh=1/15
@time OneSolution(1.0/15.0,1.0,false,1e-10,branch=1,ToPrint=false,ergosphere=false,petrov=false,light_ring=false,isco=false,sphericity=false,linvel=false)

nothing