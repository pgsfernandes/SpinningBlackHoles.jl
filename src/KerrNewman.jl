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