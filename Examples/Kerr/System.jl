#DEFINE RESOLUTION IN x
global const Nx=40
#DEFINE RESOLUTION IN y
global const Ny=8

using SpinningBlackHoles, NLsolve, DelimitedFiles, Cubature

LoadSystem()

global const NFields=4
#DEFINE BOUNDARY CONDITIONS AT THE HORIZON
function BC_Horizon!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,R::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64)
    R[idx+0] = f - 2*dfdx
    R[idx+1] = g + 2*dgdx
    R[idx+2] = dhdx
    R[idx+3] = -WBC + W
end

#DEFINE BOUNDARY CONDITIONS AT THE HORIZON - SPIN AS INPUT PARAMETER
function BC_Horizon_Spin!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,R::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64)
    R[idx+0] = f - 2*dfdx
    R[idx+1] = g + 2*dgdx
    R[idx+2] = dhdx
    R[idx+3] = W - dWdx
end

#DEFINE BOUNDARY CONDITIONS AT THE INFINITY
function BC_Infinity!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,R::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64)
    R[idx+0] = -1 + f
    R[idx+1] = -1 + g
    R[idx+2] = -1 + h
    R[idx+3] = W
end

#DEFINE BOUNDARY CONDITIONS AT THE INFINITY - SPIN AS INPUT PARAMETER
function BC_Infinity_Spin!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,R::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64)
    R[idx+0] = -1 + f
    R[idx+1] = -1 + g
    R[idx+2] = -1 + h
    R[idx+3] = WBC*(1 + dfdx)^2 + dWdx
end

#DEFINE THE FIELD EQUATIONS
function Field_Eqs!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,haxis::Float64,siny::Float64,cosy::Float64,R::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,x::Float64,y::Float64)
    R[idx+0] = -((-1 + x)^4*(((1 + x)*f*((1 + x)*dfdy*dgdy + (-1 + x)^2*(2*f + (1 + x)*dfdx)*dgdx))/(2*(-1 + x)^2) - g^2*siny^2*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) + ((1 + x)*g*(-((1 + x)*(dfdy^2 + (-1 + x)^2*dfdx^2)) + f*((1 + x)*cosy/siny*dfdy + (1 + x)*d2fdy2 + (-1 + x)^2*(dfdx + (1 + x)*d2fdx2))))/(-1 + x)^2))

    R[idx+1] = -((-1 + x)^2*(-4*(1 + x)*g*siny*(dfdy*dWdy + (-1 + x)*dfdx*(2*W + (-1 + x)*dWdx)) + f*(3*(1 + x)*siny*(dgdy*dWdy + (-1 + x)*dgdx*(2*W + (-1 + x)*dWdx)) + 2*g*(-4*x*siny*W + 3*(1 + x)*cosy*dWdy + siny*((1 + x)*d2Wdy2 + (-1 + x)*((3 + x)*dWdx + (-1 + x^2)*d2Wdx2))))))

    R[idx+2] = (-1 + x)^2*f*(-4*(-1 + x)*g^2 - (1 + x)*(dgdy^2 + (-1 + x)^2*dgdx^2) + 2*g*(2*(1 + x)*cosy/siny*dgdy + (1 + x)*d2gdy2 + (-1 + x)*((-3 + x)*dgdx + (-1 + x^2)*d2gdx2)))

    if j==2
        R[idx+3] = haxis-1
    else
        R[idx+3] = -((-1 + x)^2*(4*(1 + x)*f^2*g*h^2*((1 + x)*cosy/siny*dgdy - 2*(-1 + x)*dgdx) + (1 + x)^2*f^2*h^2*(dgdy^2 + (-1 + x)^2*dgdx^2) + 3*(-1 + x)^2*g^3*h^2*siny^2*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) - (1 + x)*g^2*(4*(-1 + x)^2*f*h^2*dfdx + (1 + x)*h^2*(dfdy^2 + (-1 + x)^2*dfdx^2) + 2*f^2*(4*(-1 + x)*h^2 - (1 + x)*(dhdy^2 + (-1 + x)^2*dhdx^2) + (1 + x)*h*(d2hdy2 + (-1 + x)*(dhdx + (-1 + x)*d2hdx2))))))
    end
end

#DEFINE BOUNDARY CONDITIONS AT THE HORIZON
function BC_Horizon_Jac!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,J::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,v::Int64,funcidx::Int64,type::Int8,ordx::Int64,ordy::Int64,Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    if funcidx==1
        J[v, idx+0] = (Mx[1, 1 + ordx, i] - 2*Mx[2, 1 + ordx, i])*My[1, 1 + ordy, j, type]
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
    end
    if funcidx==2
        J[v, idx+0] = 0
        J[v, idx+1] = (Mx[1, 1 + ordx, i] + 2*Mx[2, 1 + ordx, i])*My[1, 1 + ordy, j, type]
        J[v, idx+2] = 0
        J[v, idx+3] = 0
    end
    if funcidx==3
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+3] = 0
    end
    if funcidx==4
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
    end
end

function BC_Horizon_Jac_Spin!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,J::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,v::Int64,funcidx::Int64,type::Int8,ordx::Int64,ordy::Int64,Mx::Array{Float64, 3}=Mx,Cy::Array{Float64, 4}=My)
    if funcidx==1
        J[v, idx+0] = (Mx[1, 1 + ordx, i] - 2*Mx[2, 1 + ordx, i])*My[1, 1 + ordy, j, type]
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
    end
    if funcidx==2
        J[v, idx+0] = 0
        J[v, idx+1] = (Mx[1, 1 + ordx, i] + 2*Mx[2, 1 + ordx, i])*My[1, 1 + ordy, j, type]
        J[v, idx+2] = 0
        J[v, idx+3] = 0
    end
    if funcidx==3
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+3] = 0
    end
    if funcidx==4
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = (Mx[1, 1 + ordx, i] - Mx[2, 1 + ordx, i])*My[1, 1 + ordy, j, type]
    end
end

#DEFINE BOUNDARY CONDITIONS AT THE INFINITY
function BC_Infinity_Jac!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,J::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,v::Int64,funcidx::Int64,type::Int8,ordx::Int64,ordy::Int64,Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    if funcidx==1
        J[v, idx+0] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
    end
    if funcidx==2
        J[v, idx+0] = 0
        J[v, idx+1] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+2] = 0
        J[v, idx+3] = 0
    end
    if funcidx==3
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+3] = 0
    end
    if funcidx==4
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
    end
end

function BC_Infinity_Jac_Spin!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,J::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,v::Int64,funcidx::Int64,type::Int8,ordx::Int64,ordy::Int64,Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    if funcidx==1
        J[v, idx+0] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 2*WBC*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(1 + dfdx)
    end
    if funcidx==2
        J[v, idx+0] = 0
        J[v, idx+1] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+2] = 0
        J[v, idx+3] = 0
    end
    if funcidx==3
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+3] = 0
    end
    if funcidx==4
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
    end
end

#DEFINE THE FIELD EQUATIONS
function Field_Eqs_Jac!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,haxis::Float64,siny::Float64,cosy::Float64,J::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,v::Int64,funcidx::Int64,type::Int8,ordx::Int64,ordy::Int64,x::Float64,y::Float64,Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    if funcidx==1
        J[v, idx+0] = -((-1 + x)^4*(1 + x)^2*f*g*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type]) - (-1 + x)^2*(1 + x)^2*f*g*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type] - (-1 + x)^4*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(((1 + x)*g*((1 + x)*cosy/siny*f - 2*(1 + x)*dfdy))/(-1 + x)^2 + ((1 + x)^2*f*dgdy)/(2*(-1 + x)^2)) - (-1 + x)^4*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(((1 + x)*g*((-1 + x)^2*f - 2*(-1 + x)^2*(1 + x)*dfdx))/(-1 + x)^2 + ((1 + x)^2*f*dgdx)/2) - (-1 + x)^4*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*((1 + x)*f*dgdx + ((1 + x)*((1 + x)*dfdy*dgdy + (-1 + x)^2*(2*f + (1 + x)*dfdx)*dgdx))/(2*(-1 + x)^2) + ((1 + x)*g*((1 + x)*cosy/siny*dfdy + (1 + x)*d2fdy2 + (-1 + x)^2*(dfdx + (1 + x)*d2fdx2)))/(-1 + x)^2)

        J[v, idx+1] = 4*(-1 + x)^2*(1 + x)*g*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*siny*dWdy + 4*(-1 + x)^3*(1 + x)*g*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny*(2*W + (-1 + x)*dWdx) - (-1 + x)^2*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(3*(1 + x)*siny*(dgdy*dWdy + (-1 + x)*dgdx*(2*W + (-1 + x)*dWdx)) + 2*g*(-4*x*siny*W + 3*(1 + x)*cosy*dWdy + siny*((1 + x)*d2Wdy2 + (-1 + x)*((3 + x)*dWdx + (-1 + x^2)*d2Wdx2))))

        J[v, idx+2] = (-1 + x)^2*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-4*(-1 + x)*g^2 - (1 + x)*(dgdy^2 + (-1 + x)^2*dgdx^2) + 2*g*(2*(1 + x)*cosy/siny*dgdy + (1 + x)*d2gdy2 + (-1 + x)*((-3 + x)*dgdx + (-1 + x^2)*d2gdx2)))

        if j==2
            J[v, idx+3] = 0
        else
            J[v, idx+3] = 2*(-1 + x)^2*(1 + x)^2*g^2*h^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*dfdy + (-1 + x)^2*(1 + x)*g^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(4*(-1 + x)^2*f*h^2 + 2*(-1 + x)^2*(1 + x)*h^2*dfdx) - (-1 + x)^2*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(8*(1 + x)*f*g*h^2*((1 + x)*cosy/siny*dgdy - 2*(-1 + x)*dgdx) + 2*(1 + x)^2*f*h^2*(dgdy^2 + (-1 + x)^2*dgdx^2) - (1 + x)*g^2*(4*(-1 + x)^2*h^2*dfdx + 4*f*(4*(-1 + x)*h^2 - (1 + x)*(dhdy^2 + (-1 + x)^2*dhdx^2) + (1 + x)*h*(d2hdy2 + (-1 + x)*(dhdx + (-1 + x)*d2hdx2)))))

        end
    end
    if funcidx==2
        J[v, idx+0] = -1/2*((-1 + x)^2*(1 + x)^2*f*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*dfdy) - ((-1 + x)^4*(1 + x)*f*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(2*f + (1 + x)*dfdx))/2 - (-1 + x)^4*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-2*g*siny^2*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) + ((1 + x)*(-((1 + x)*(dfdy^2 + (-1 + x)^2*dfdx^2)) + f*((1 + x)*cosy/siny*dfdy + (1 + x)*d2fdy2 + (-1 + x)^2*(dfdx + (1 + x)*d2fdx2))))/(-1 + x)^2)

        J[v, idx+1] = -3*(-1 + x)^2*(1 + x)*f*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*siny*dWdy - 3*(-1 + x)^3*(1 + x)*f*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny*(2*W + (-1 + x)*dWdx) - (-1 + x)^2*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-4*(1 + x)*siny*(dfdy*dWdy + (-1 + x)*dfdx*(2*W + (-1 + x)*dWdx)) + 2*f*(-4*x*siny*W + 3*(1 + x)*cosy*dWdy + siny*((1 + x)*d2Wdy2 + (-1 + x)*((3 + x)*dWdx + (-1 + x^2)*d2Wdx2))))

        J[v, idx+2] = 2*(-1 + x)^3*(-1 + x^2)*f*g*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type] + 2*(-1 + x)^2*(1 + x)*f*g*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type] + (-1 + x)^2*f*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(4*(1 + x)*cosy/siny*g - 2*(1 + x)*dgdy) + (-1 + x)^2*f*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(2*(-3 + x)*(-1 + x)*g - 2*(-1 + x)^2*(1 + x)*dgdx) + (-1 + x)^2*f*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-8*(-1 + x)*g + 2*(2*(1 + x)*cosy/siny*dgdy + (1 + x)*d2gdy2 + (-1 + x)*((-3 + x)*dgdx + (-1 + x^2)*d2gdx2)))

        if j==2
            J[v, idx+3] = 0
        else
            J[v, idx+3] = -((-1 + x)^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(4*(1 + x)^2*cosy/siny*f^2*g*h^2 + 2*(1 + x)^2*f^2*h^2*dgdy)) - (-1 + x)^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-8*(-1 + x)*(1 + x)*f^2*g*h^2 + 2*(-1 + x)^2*(1 + x)^2*f^2*h^2*dgdx) - (-1 + x)^2*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(4*(1 + x)*f^2*h^2*((1 + x)*cosy/siny*dgdy - 2*(-1 + x)*dgdx) + 9*(-1 + x)^2*g^2*h^2*siny^2*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) - 2*(1 + x)*g*(4*(-1 + x)^2*f*h^2*dfdx + (1 + x)*h^2*(dfdy^2 + (-1 + x)^2*dfdx^2) + 2*f^2*(4*(-1 + x)*h^2 - (1 + x)*(dhdy^2 + (-1 + x)^2*dhdx^2) + (1 + x)*h*(d2hdy2 + (-1 + x)*(dhdx + (-1 + x)*d2hdx2)))))

        end
    end
    if funcidx==3
        J[v, idx+0] = 0

        J[v, idx+1] = 0

        J[v, idx+2] = 0

        if j==2
            J[v, idx+3] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, 1, 1]
        else
            J[v, idx+3] = 2*(-1 + x)^4*(1 + x)^2*f^2*g^2*h*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type] + 2*(-1 + x)^2*(1 + x)^2*f^2*g^2*h*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type] - 4*(-1 + x)^2*(1 + x)^2*f^2*g^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*dhdy + 2*(-1 + x)^2*(1 + x)*f^2*g^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*((-1 + x)*(1 + x)*h - 2*(-1 + x)^2*(1 + x)*dhdx) - (-1 + x)^2*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(8*(1 + x)*f^2*g*h*((1 + x)*cosy/siny*dgdy - 2*(-1 + x)*dgdx) + 2*(1 + x)^2*f^2*h*(dgdy^2 + (-1 + x)^2*dgdx^2) + 6*(-1 + x)^2*g^3*h*siny^2*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) - (1 + x)*g^2*(8*(-1 + x)^2*f*h*dfdx + 2*(1 + x)*h*(dfdy^2 + (-1 + x)^2*dfdx^2) + 2*f^2*(8*(-1 + x)*h + (1 + x)*(d2hdy2 + (-1 + x)*(dhdx + (-1 + x)*d2hdx2)))))

        end
    end
    if funcidx==4
        J[v, idx+0] = 2*(-1 + x)^4*g^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*siny^2*dWdy + (-1 + x)^4*g^2*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny^2*(8*W + 4*(-1 + x)*dWdx) + (-1 + x)^4*g^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny^2*(4*(-1 + x)*W + 2*(-1 + x)^2*dWdx)

        J[v, idx+1] = -2*(-1 + x)^3*(-1 + x^2)*f*g*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny - 2*(-1 + x)^2*(1 + x)*f*g*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type]*siny - (-1 + x)^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(-4*(1 + x)*g*siny*dfdy + f*(6*(1 + x)*cosy*g + 3*(1 + x)*siny*dgdy)) - (-1 + x)^2*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-8*(-1 + x)*(1 + x)*g*siny*dfdx + f*(-8*x*g*siny + 6*(-1 + x)*(1 + x)*siny*dgdx)) - (-1 + x)^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-4*(-1 + x)^2*(1 + x)*g*siny*dfdx + f*(2*(-1 + x)*(3 + x)*g*siny + 3*(-1 + x)^2*(1 + x)*siny*dgdx))

        J[v, idx+2] = 0

        if j==2
            J[v, idx+3] = 0
        else
            J[v, idx+3] = -6*(-1 + x)^4*g^3*h^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*siny^2*dWdy - 3*(-1 + x)^4*g^3*h^2*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny^2*(8*W + 4*(-1 + x)*dWdx) - 3*(-1 + x)^4*g^3*h^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny^2*(4*(-1 + x)*W + 2*(-1 + x)^2*dWdx)

        end
    end
end
#DEFINE THE WHOLE SYSTEM (BCs+FIELD EQS) ON OUR GRID
function residual_jac!(R::Matrix{Float64},J::Matrix{Float64},a::Matrix{Float64},f::Field,g::Field,h::Field,W::Field,WBC::Float64,rh::Float64,spin::Bool,x::Vector{Float64}=X,y::Vector{Float64}=Y)

    f.a=a[0*Nx+1:1*Nx,:]

    g.a=a[1*Nx+1:2*Nx,:]

    h.a=a[2*Nx+1:3*Nx,:]

    W.a=a[3*Nx+1:4*Nx,:]

    type=[f.type,g.type,h.type,W.type]
    ff=Matrix{Float64}(undef,Nx,Ny); dfdx=Matrix{Float64}(undef,Nx,Ny); d2fdx2=Matrix{Float64}(undef,Nx,Ny); dfdy=Matrix{Float64}(undef,Nx,Ny); d2fdy2=Matrix{Float64}(undef,Nx,Ny); d2fdxdy=Matrix{Float64}(undef,Nx,Ny);
    gg=Matrix{Float64}(undef,Nx,Ny); dgdx=Matrix{Float64}(undef,Nx,Ny); d2gdx2=Matrix{Float64}(undef,Nx,Ny); dgdy=Matrix{Float64}(undef,Nx,Ny); d2gdy2=Matrix{Float64}(undef,Nx,Ny); d2gdxdy=Matrix{Float64}(undef,Nx,Ny);
    hh=Matrix{Float64}(undef,Nx,Ny); dhdx=Matrix{Float64}(undef,Nx,Ny); d2hdx2=Matrix{Float64}(undef,Nx,Ny); dhdy=Matrix{Float64}(undef,Nx,Ny); d2hdy2=Matrix{Float64}(undef,Nx,Ny); d2hdxdy=Matrix{Float64}(undef,Nx,Ny);
    WW=Matrix{Float64}(undef,Nx,Ny); dWdx=Matrix{Float64}(undef,Nx,Ny); d2Wdx2=Matrix{Float64}(undef,Nx,Ny); dWdy=Matrix{Float64}(undef,Nx,Ny); d2Wdy2=Matrix{Float64}(undef,Nx,Ny); d2Wdxdy=Matrix{Float64}(undef,Nx,Ny);
    haxis=Array{Float64}(undef,Nx);

    for i in 1:Nx
        for j in 1:Ny
            ff[i,j]=f(i,j+1)
            dfdx[i,j]=f(i,j+1,dx=1)
            d2fdx2[i,j]=f(i,j+1,dx=2)
            dfdy[i,j]=f(i,j+1,dy=1)
            d2fdy2[i,j]=f(i,j+1,dy=2)
            d2fdxdy[i,j]=f(i,j+1,dx=1,dy=1)
            gg[i,j]=g(i,j+1)
            dgdx[i,j]=g(i,j+1,dx=1)
            d2gdx2[i,j]=g(i,j+1,dx=2)
            dgdy[i,j]=g(i,j+1,dy=1)
            d2gdy2[i,j]=g(i,j+1,dy=2)
            d2gdxdy[i,j]=g(i,j+1,dx=1,dy=1)
            hh[i,j]=h(i,j+1)
            dhdx[i,j]=h(i,j+1,dx=1)
            d2hdx2[i,j]=h(i,j+1,dx=2)
            dhdy[i,j]=h(i,j+1,dy=1)
            d2hdy2[i,j]=h(i,j+1,dy=2)
            d2hdxdy[i,j]=h(i,j+1,dx=1,dy=1)
            WW[i,j]=W(i,j+1)
            dWdx[i,j]=W(i,j+1,dx=1)
            d2Wdx2[i,j]=W(i,j+1,dx=2)
            dWdy[i,j]=W(i,j+1,dy=1)
            d2Wdy2[i,j]=W(i,j+1,dy=2)
            d2Wdxdy[i,j]=W(i,j+1,dx=1,dy=1)
        end
        haxis[i]=h(i,1)
    end

    siny=Array{Float64}(undef,Ny);
    cosy=Array{Float64}(undef,Ny);
    for j in 1:Ny
        siny[j]=sin(y[j+1])
        cosy[j]=cos(y[j+1])
    end

    if !(R == nothing)
        #INDEX FOR THE RESIDUAL R
        idx=1
        #LOOP ON ALL X POINTS INCLUDING -1 AND 1
        for i in 1:Nx
        #LOOP ON INTERNAL Y POINTS. NOTE THAT BC ON 0 AND PI/2 ARE AUTOMATICALLY IMPOSED BY OUR CHOICE OF BASIS FUNCTIONS (EVEN COSINES)
            for j in 1:Ny
                if i==1
                    #DEFINE BCS AT THE HORIZON
                    if !spin
                        BC_Horizon!(i,j,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],R,idx,WBC,rh)
                    else
                        BC_Horizon_Spin!(i,j,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],R,idx,WBC,rh)
                    end
                elseif i==Nx
                    #DEFINE BCS AT INFINITY
                    if !spin
                        BC_Infinity!(i,j,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],R,idx,WBC,rh)
                    else
                        BC_Infinity_Spin!(i,j,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],R,idx,WBC,rh)
                    end
                else
                    #DEFINE FIELD EQUATIONS EVERYWHERE ELSE ON THE GRID
                    Field_Eqs!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],haxis[i],siny[j],cosy[j],R,idx,WBC,rh,X[i],Y[j])
                end
                idx+=NFields
            end
        end
    end

    if !(J == nothing)
        funcidx=1
        ordx=0
        ordy=0
        for v in 1:Nx*Ny*NFields
            #INDEX FOR THE JACOBIAN J COLUMNS
            idx=1
            #LOOP ON ALL X POINTS INCLUDING -1 AND 1
            for i in 1:Nx
            #LOOP ON INTERNAL Y POINTS. NOTE THAT BC ON 0 AND PI/2 ARE AUTOMATICALLY IMPOSED BY OUR CHOICE OF BASIS FUNCTIONS (EVEN COSINES)
                for j in 1:Ny
                    if i==1
                        #DEFINE BCS AT THE HORIZON
                        if !spin
                            BC_Horizon_Jac!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],J,idx,WBC,rh,v,funcidx,type[funcidx],ordx,ordy)
                        else
                            BC_Horizon_Jac_Spin!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],J,idx,WBC,rh,v,funcidx,type[funcidx],ordx,ordy)
                        end
                    elseif i==Nx
                        #DEFINE BCS AT INFINITY
                        if !spin
                            BC_Infinity_Jac!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],J,idx,WBC,rh,v,funcidx,type[funcidx],ordx,ordy)
                        else
                            BC_Infinity_Jac_Spin!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],J,idx,WBC,rh,v,funcidx,type[funcidx],ordx,ordy)
                        end
                    else
                        #DEFINE FIELD EQUATIONS EVERYWHERE ELSE ON THE GRID
                        Field_Eqs_Jac!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],haxis[i],siny[j],cosy[j],J,idx,WBC,rh,v,funcidx,type[funcidx],ordx,ordy,X[i],Y[j])
                    end
                    idx+=NFields
                end
            end
            ordx+=1
            if mod(v,Nx)==0
                funcidx+=1
                ordx=0
            end
            if mod(v,Nx*NFields)==0
                funcidx=1
                ordx=0
                ordy+=1
            end
        end
        J .= transpose(J)
    end
end

#SOLVE THE NON-LINEAR SYSTEM TO OBTAIN THE SOLUTIONS
function solve_system(f::Field,g::Field,h::Field,W::Field,tol::Float64,WBC::Float64,rh::Float64,spin::Bool=false,show_trace::Bool=true,iterations::Int64=30,method = :newton)
    return nlsolve(only_fj!((R,J,a)->residual_jac!(R,J,a,f,g,h,W,WBC,rh,spin)), vcat(f.a,g.a,h.a,W.a), ftol=0.0,xtol=tol, show_trace=show_trace, iterations = iterations, method = method)
end
nothing