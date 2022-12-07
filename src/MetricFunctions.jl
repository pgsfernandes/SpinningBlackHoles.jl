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