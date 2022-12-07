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