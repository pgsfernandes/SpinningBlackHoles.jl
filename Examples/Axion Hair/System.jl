#DEFINE RESOLUTION IN x
global const Nx=40
#DEFINE RESOLUTION IN y
global const Ny=8
using SpinningBlackHoles, NLsolve, DelimitedFiles, Cubature, Roots

LoadSystem()

global const NFields=7
#DEFINE BOUNDARY CONDITIONS AT THE HORIZON
function BC_Horizon!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,p::Float64,dpdx::Float64,d2pdx2::Float64,dpdy::Float64,d2pdy2::Float64,d2pdxdy::Float64,A::Float64,dAdx::Float64,d2Adx2::Float64,dAdy::Float64,d2Ady2::Float64,d2Adxdy::Float64,B::Float64,dBdx::Float64,d2Bdx2::Float64,dBdy::Float64,d2Bdy2::Float64,d2Bdxdy::Float64,R::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,q::Float64,k::Float64)
    R[idx+0] = f - 2*dfdx
    R[idx+1] = g + 2*dgdx
    R[idx+2] = dhdx
    R[idx+3] = -WBC + W
    R[idx+4] = dpdx
    R[idx+5] = A
    R[idx+6] = dBdx
end

#DEFINE BOUNDARY CONDITIONS AT THE HORIZON - SPIN AS INPUT PARAMETER
function BC_Horizon_Spin!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,p::Float64,dpdx::Float64,d2pdx2::Float64,dpdy::Float64,d2pdy2::Float64,d2pdxdy::Float64,A::Float64,dAdx::Float64,d2Adx2::Float64,dAdy::Float64,d2Ady2::Float64,d2Adxdy::Float64,B::Float64,dBdx::Float64,d2Bdx2::Float64,dBdy::Float64,d2Bdy2::Float64,d2Bdxdy::Float64,R::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,q::Float64,k::Float64)
    R[idx+0] = f - 2*dfdx
    R[idx+1] = g + 2*dgdx
    R[idx+2] = dhdx
    R[idx+3] = W - dWdx
    R[idx+4] = dpdx
    R[idx+5] = A
    R[idx+6] = dBdx
end

#DEFINE BOUNDARY CONDITIONS AT THE INFINITY
function BC_Infinity!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,p::Float64,dpdx::Float64,d2pdx2::Float64,dpdy::Float64,d2pdy2::Float64,d2pdxdy::Float64,A::Float64,dAdx::Float64,d2Adx2::Float64,dAdy::Float64,d2Ady2::Float64,d2Adxdy::Float64,B::Float64,dBdx::Float64,d2Bdx2::Float64,dBdy::Float64,d2Bdy2::Float64,d2Bdxdy::Float64,R::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,q::Float64,k::Float64)
    R[idx+0] = -1 + f
    R[idx+1] = -1 + g
    R[idx+2] = -1 + h
    R[idx+3] = W
    R[idx+4] = p
    R[idx+5] = dAdx - (q*(1 + dfdx))/2
    R[idx+6] = B
end

#DEFINE BOUNDARY CONDITIONS AT THE INFINITY - SPIN AS INPUT PARAMETER
function BC_Infinity_Spin!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,p::Float64,dpdx::Float64,d2pdx2::Float64,dpdy::Float64,d2pdy2::Float64,d2pdxdy::Float64,A::Float64,dAdx::Float64,d2Adx2::Float64,dAdy::Float64,d2Ady2::Float64,d2Adxdy::Float64,B::Float64,dBdx::Float64,d2Bdx2::Float64,dBdy::Float64,d2Bdy2::Float64,d2Bdxdy::Float64,R::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,q::Float64,k::Float64)
    R[idx+0] = -1 + f
    R[idx+1] = -1 + g
    R[idx+2] = -1 + h
    R[idx+3] = WBC*(1 + dfdx)^2 + dWdx
    R[idx+4] = p
    R[idx+5] = dAdx - (q*(1 + dfdx))/2
    R[idx+6] = B
end

#DEFINE THE FIELD EQUATIONS
function Field_Eqs!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,p::Float64,dpdx::Float64,d2pdx2::Float64,dpdy::Float64,d2pdy2::Float64,d2pdxdy::Float64,A::Float64,dAdx::Float64,d2Adx2::Float64,dAdy::Float64,d2Ady2::Float64,d2Adxdy::Float64,B::Float64,dBdx::Float64,d2Bdx2::Float64,dBdy::Float64,d2Bdy2::Float64,d2Bdxdy::Float64,haxis::Float64,siny::Float64,cosy::Float64,d0V::Float64,d1V::Float64,d2V::Float64,d3V::Float64,R::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,q::Float64,k::Float64,x::Float64,y::Float64)
    R[idx+0] = (-4*rh^2*(1 + x)*g*(dAdy*dfdy + (-1 + x)^2*dAdx*dfdx))/(-1 + x)^2 - k*rh*(1 + x)^2*f^2*g^(1/2)*siny*(dpdy*dBdx - dBdy*dpdx) + (2*f*((1 + x)*(rh - rh*x)^2*(dAdy*dgdy + (-1 + x)^2*dAdx*dgdx) + (rh*g*(-4*rh*(-1 + x)^4*dAdx + (1 + x)*(4*rh*(-1 + x)^2*cosy/siny*dAdy - (-1 + x)^4*siny^2*dBdy*dWdy + 4*rh*(-1 + x)^2*d2Ady2 - 8*rh*(-1 + x)^3*dAdx - 2*(-1 + x)^5*siny^2*W*dBdx - (-1 + x)^6*siny^2*dBdx*dWdx - 4*rh*(1 - x)^3*(2*dAdx + (-1 + x)*d2Adx2))))/2))/(-1 + x)^4 + (rh*B*(4*k*(1 + x)^2*cosy*f^2*g^(1/2)*dpdx + 2*(1 + x)*g*siny^2*(dfdy*dWdy + (-1 + x)*dfdx*(2*W + (-1 + x)*dWdx)) - f*siny*((1 + x)*siny*(dgdy*dWdy + (-1 + x)*dgdx*(2*W + (-1 + x)*dWdx)) + 2*g*(4*siny*W + 3*(1 + x)*cosy*dWdy + siny*((1 + x)*d2Wdy2 + (-1 + x)*((5 + 3*x)*dWdx + (-1 + x^2)*d2Wdx2))))))/2

    R[idx+1] = -4*rh^2*(-1 + x)^2*(1 + x)*g^2*W*(dAdy*dfdy + (-1 + x)^2*dAdx*dfdx) - (-1 + x)^2*(1 + x)^2*f^2*g*(rh*(1 + x)*(dBdy*dfdy + (-1 + x)^2*dBdx*dfdx) - k*g^(1/2)*siny*(rh*dpdy*(4*rh*1/siny^2*dAdx - (-1 + x)^2*W*dBdx) + rh*(-4*rh*1/siny^2*dAdy + (-1 + x)^2*W*dBdy)*dpdx)) + 2*f*g*((1 + x)*(rh - rh*x)^2*W*(dAdy*dgdy + (-1 + x)^2*dAdx*dgdx) + g*(-(rh*(-1 + x)^5*(1 + x)*siny^2*W^2*dBdx) + 2*(1 + x)*(rh - rh*x)^2*(dAdy*dWdy + (-1 + x)^2*dAdx*dWdx) + (rh*W*(-4*rh*(-1 + x)^4*dAdx + (1 + x)*(4*rh*(-1 + x)^2*cosy/siny*dAdy - (-1 + x)^4*siny^2*dBdy*dWdy + 4*rh*(-1 + x)^2*d2Ady2 - (-1 + x)^6*siny^2*dBdx*dWdx - 4*rh*(1 - x)^3*(2*dAdx + (-1 + x)*d2Adx2))))/2)) + rh*(1 + x)^2*f^3*(((-1 + x)^2*(1 + x)*(dBdy*dgdy + (-1 + x)^2*dBdx*dgdx))/2 - g*((-1 + x)^4*dBdx + (-1 + x)^2*(1 + x)*(3*cosy/siny*dBdy + d2Bdy2 + (-1 + x)*(2*dBdx + (-1 + x)*d2Bdx2)))) + (-1 + x)^2*B*(rh*(1 + x)^3*f^3*(2*g + cosy/siny*dgdy) + rh*(-1 + x)^2*(1 + x)*g^2*siny^2*W*(dfdy*dWdy + (-1 + x)*dfdx*(2*W + (-1 + x)*dWdx)) - (1 + x)^2*f^2*g*(2*rh*(1 + x)*cosy/siny*dfdy - k*rh*g^(1/2)*(2*(1 - x)*W*(siny*dpdy - (-1 + x)*cosy*dpdx) + (-1 + x)^2*siny*(dWdy*dpdx - dpdy*dWdx))) - (rh*(-1 + x)^2*f*g*siny*((1 + x)*siny*W*(dgdy*dWdy + (-1 + x)*dgdx*(2*W + (-1 + x)*dWdx)) + 2*g*(4*(2 + x)*siny*W^2 - (-1 + x)^2*siny*W*dWdx + (1 + x)*siny*(dWdy^2 + (-1 + x)^2*dWdx^2) + (1 + x)*W*(3*cosy*dWdy + siny*(d2Wdy2 + (-1 + x)*(8*dWdx + (-1 + x)*d2Wdx2))))))/2)

    R[idx+2] = -d1V - (k*(-1 + x)^4*f*(4*rh*siny*(-(dBdy*dAdx) + dAdy*dBdx) + (-1 + x)*B^2*siny*2*siny*cosy*(2*W + (-1 + x)*dWdx) + B*(2*(-1 + x)*siny^3*W*dBdy - 8*rh*cosy*dAdx - (-1 + x)^2*siny^3*(dWdy*dBdx - dBdy*dWdx))))/(16*rh^4*(1 + x)*g^(3/2)*h) + (f*(((-1 + x)^2*(1 + x)*(dgdy*dpdy + (-1 + x)^2*dgdx*dpdx))/2 + g*((-1 + x)^4*dpdx + (-1 + x)^2*(1 + x)*(cosy/siny*dpdy + d2pdy2 + (-1 + x)^2*d2pdx2))))/(4*rh^2*(1 + x)*g^2*h)

    R[idx+3] = 2*(1 + x)*(rh - rh*x)^2*f*(2*(-1 + x)^2*(1 + x)*B^2*cosy^2*f^2 + (-1 + x)^2*(1 + x)*B*f^2*2*siny*cosy*dBdy + ((-1 + x)^2*(1 + x)*f^2*siny^2*(dBdy^2 + (-1 + x)^2*dBdx^2))/2 - rh^2*(-1 + x)^2*f*dgdx - (rh^2*(1 + x)*(dfdy*dgdy + (-1 + x)^2*dfdx*dgdx))/2) - 2*rh^4*g^2*(16*rh^2*(1 + x)^2*f*h*d0V - 4*(-1 + x)^4*siny^2*W^2 - 4*(-1 + x)^5*siny^2*W*dWdx - (-1 + x)^4*siny^2*(dWdy^2 + (-1 + x)^2*dWdx^2)) + rh^2*g*((-1 + x)^6*B^2*f*siny^4*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) - 8*rh*(-1 + x)^4*B*f*siny^2*(dAdy*dWdy + (-1 + x)*dAdx*(2*W + (-1 + x)*dWdx)) + 4*rh^2*(((-1 + x^2)^2*(dfdy^2 + (-1 + x)^2*dfdx^2))/2 + f*(4*(-1 + x)^2*dAdy^2 + 4*(-1 + x)^4*dAdx^2 - ((-1 + x)^4*(1 + x)*dfdx)/2 - ((-1 + x^2)^2*(cosy/siny*dfdy + d2fdy2 + (-1 + x)^2*d2fdx2))/2)))

    R[idx+4] = (8*rh^3*(1 + x)*g*siny*(dfdy*dWdy + (-1 + x)*dfdx*(2*W + (-1 + x)*dWdx)))/(-1 + x)^2 + 4*rh*(1 + x)*f^2*(-2*B^2*cosy*siny^2*dWdy + (4*rh*siny*(dAdy*dBdy + (-1 + x)^2*dAdx*dBdx))/(-1 + x)^2 + B*((8*rh*cosy*dAdy)/(-1 + x)^2 - siny^3*(dBdy*dWdy + (-1 + x)*dBdx*(2*W + (-1 + x)*dWdx)))) - (2*rh^3*f*(3*(1 + x)*siny*(dgdy*dWdy + (-1 + x)*dgdx*(2*W + (-1 + x)*dWdx)) + 2*g*(-4*x*siny*W + 3*(1 + x)*cosy*dWdy + siny*((1 + x)*d2Wdy2 + (-1 + x)*((3 + x)*dWdx + (-1 + x^2)*d2Wdx2)))))/(-1 + x)^2

    R[idx+5] = -2*(-1 + x)*f*g^2 + (32*rh^2*(1 + x)*g^3*h*d0V)/(-1 + x)^2 - ((1 + x)*f*(dgdy^2 + (-1 + x)^2*dgdx^2))/2 + f*g*(2*(1 + x)*cosy/siny*dgdy + (1 + x)*d2gdy2 + (-1 + x)*((-3 + x)*dgdx + (-1 + x^2)*d2gdx2))

    if j==2
        R[idx+6] = haxis-1
    else
        R[idx+6] = -(rh^2*(-1 + x)^2*(1 + x)^2*f^2*h^2*(dgdy^2 + (-1 + x)^2*dgdx^2)) - 2*(-1 + x)^2*(1 + x)*f^2*g*h^2*(2*(-1 + x)^2*(1 + x)*B^2*cosy^2*f + (-1 + x)^2*(1 + x)*B*f*2*siny*cosy*dBdy + ((-1 + x)^2*(1 + x)*f*siny^2*(dBdy^2 + (-1 + x)^2*dBdx^2))/2 + 2*rh^2*(-1 + x)^2*dgdx - 2*rh^2*(1 + x)*(-(cosy/siny*dgdy) + (-1 + x)*dgdx)) - rh^2*g^3*h^2*(32*rh^2*(1 + x)^2*f*h*d0V + 12*(-1 + x)^4*siny^2*W^2 + 12*(-1 + x)^5*siny^2*W*dWdx + 3*(-1 + x)^4*siny^2*(dWdy^2 + (-1 + x)^2*dWdx^2)) - g^2*((-1 + x)^6*B^2*f*h^2*siny^4*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) - 8*rh*(-1 + x)^4*B*f*h^2*siny^2*(dAdy*dWdy + (-1 + x)*dAdx*(2*W + (-1 + x)*dWdx)) - rh^2*(-1 + x)^2*((1 + x)^2*h^2*(dfdy^2 + (-1 + x)^2*dfdx^2) - 16*f*h^2*(dAdy^2 - ((-1 + x)^2*(-4*dAdx^2 + (1 + x)*dfdx))/4) - 2*(1 + x)*f^2*((1 + x)*(dhdy^2 + (-1 + x)^2*dhdx^2) + 4*h^2*(1 - x - ((1 + x)*(dpdy^2 + (-1 + x)^2*dpdx^2))/2) - (1 + x)*h*(d2hdy2 + (-1 + x)*(dhdx + (-1 + x)*d2hdx2)))))
    end
end

#DEFINE BOUNDARY CONDITIONS AT THE HORIZON
function BC_Horizon_Jac!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,p::Float64,dpdx::Float64,d2pdx2::Float64,dpdy::Float64,d2pdy2::Float64,d2pdxdy::Float64,A::Float64,dAdx::Float64,d2Adx2::Float64,dAdy::Float64,d2Ady2::Float64,d2Adxdy::Float64,B::Float64,dBdx::Float64,d2Bdx2::Float64,dBdy::Float64,d2Bdy2::Float64,d2Bdxdy::Float64,J::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,q::Float64,k::Float64,v::Int64,funcidx::Int64,type::Int8,ordx::Int64,ordy::Int64,Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    if funcidx==1
        J[v, idx+0] = (Mx[1, 1 + ordx, i] - 2*Mx[2, 1 + ordx, i])*My[1, 1 + ordy, j, type]
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==2
        J[v, idx+0] = 0
        J[v, idx+1] = (Mx[1, 1 + ordx, i] + 2*Mx[2, 1 + ordx, i])*My[1, 1 + ordy, j, type]
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==3
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==4
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==5
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==6
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+6] = 0
    end
    if funcidx==7
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
    end
end

function BC_Horizon_Jac_Spin!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,p::Float64,dpdx::Float64,d2pdx2::Float64,dpdy::Float64,d2pdy2::Float64,d2pdxdy::Float64,A::Float64,dAdx::Float64,d2Adx2::Float64,dAdy::Float64,d2Ady2::Float64,d2Adxdy::Float64,B::Float64,dBdx::Float64,d2Bdx2::Float64,dBdy::Float64,d2Bdy2::Float64,d2Bdxdy::Float64,J::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,q::Float64,k::Float64,v::Int64,funcidx::Int64,type::Int8,ordx::Int64,ordy::Int64,Mx::Array{Float64, 3}=Mx,Cy::Array{Float64, 4}=My)
    if funcidx==1
        J[v, idx+0] = (Mx[1, 1 + ordx, i] - 2*Mx[2, 1 + ordx, i])*My[1, 1 + ordy, j, type]
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==2
        J[v, idx+0] = 0
        J[v, idx+1] = (Mx[1, 1 + ordx, i] + 2*Mx[2, 1 + ordx, i])*My[1, 1 + ordy, j, type]
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==3
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==4
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = (Mx[1, 1 + ordx, i] - Mx[2, 1 + ordx, i])*My[1, 1 + ordy, j, type]
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==5
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==6
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+6] = 0
    end
    if funcidx==7
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
    end
end

#DEFINE BOUNDARY CONDITIONS AT THE INFINITY
function BC_Infinity_Jac!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,p::Float64,dpdx::Float64,d2pdx2::Float64,dpdy::Float64,d2pdy2::Float64,d2pdxdy::Float64,A::Float64,dAdx::Float64,d2Adx2::Float64,dAdy::Float64,d2Ady2::Float64,d2Adxdy::Float64,B::Float64,dBdx::Float64,d2Bdx2::Float64,dBdy::Float64,d2Bdy2::Float64,d2Bdxdy::Float64,J::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,q::Float64,k::Float64,v::Int64,funcidx::Int64,type::Int8,ordx::Int64,ordy::Int64,Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    if funcidx==1
        J[v, idx+0] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = -1/2*(q*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type])
        J[v, idx+6] = 0
    end
    if funcidx==2
        J[v, idx+0] = 0
        J[v, idx+1] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==3
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==4
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==5
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==6
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+6] = 0
    end
    if funcidx==7
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
    end
end

function BC_Infinity_Jac_Spin!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,p::Float64,dpdx::Float64,d2pdx2::Float64,dpdy::Float64,d2pdy2::Float64,d2pdxdy::Float64,A::Float64,dAdx::Float64,d2Adx2::Float64,dAdy::Float64,d2Ady2::Float64,d2Adxdy::Float64,B::Float64,dBdx::Float64,d2Bdx2::Float64,dBdy::Float64,d2Bdy2::Float64,d2Bdxdy::Float64,J::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,q::Float64,k::Float64,v::Int64,funcidx::Int64,type::Int8,ordx::Int64,ordy::Int64,Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    if funcidx==1
        J[v, idx+0] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 2*WBC*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(1 + dfdx)
        J[v, idx+4] = 0
        J[v, idx+5] = -1/2*(q*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type])
        J[v, idx+6] = 0
    end
    if funcidx==2
        J[v, idx+0] = 0
        J[v, idx+1] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==3
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==4
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==5
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+5] = 0
        J[v, idx+6] = 0
    end
    if funcidx==6
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]
        J[v, idx+6] = 0
    end
    if funcidx==7
        J[v, idx+0] = 0
        J[v, idx+1] = 0
        J[v, idx+2] = 0
        J[v, idx+3] = 0
        J[v, idx+4] = 0
        J[v, idx+5] = 0
        J[v, idx+6] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]
    end
end

#DEFINE THE FIELD EQUATIONS
function Field_Eqs_Jac!(i::Int64,j::Int64,f::Float64,dfdx::Float64,d2fdx2::Float64,dfdy::Float64,d2fdy2::Float64,d2fdxdy::Float64,g::Float64,dgdx::Float64,d2gdx2::Float64,dgdy::Float64,d2gdy2::Float64,d2gdxdy::Float64,h::Float64,dhdx::Float64,d2hdx2::Float64,dhdy::Float64,d2hdy2::Float64,d2hdxdy::Float64,W::Float64,dWdx::Float64,d2Wdx2::Float64,dWdy::Float64,d2Wdy2::Float64,d2Wdxdy::Float64,p::Float64,dpdx::Float64,d2pdx2::Float64,dpdy::Float64,d2pdy2::Float64,d2pdxdy::Float64,A::Float64,dAdx::Float64,d2Adx2::Float64,dAdy::Float64,d2Ady2::Float64,d2Adxdy::Float64,B::Float64,dBdx::Float64,d2Bdx2::Float64,dBdy::Float64,d2Bdy2::Float64,d2Bdxdy::Float64,haxis::Float64,siny::Float64,cosy::Float64,d0V::Float64,d1V::Float64,d2V::Float64,d3V::Float64,J::Matrix{Float64},idx::Int64,WBC::Float64,rh::Float64,q::Float64,k::Float64,v::Int64,funcidx::Int64,type::Int8,ordx::Int64,ordy::Int64,x::Float64,y::Float64,Mx::Array{Float64, 3}=Mx,My::Array{Float64, 4}=My)
    if funcidx==1
        J[v, idx+0] = Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*((-4*rh^2*(1 + x)*g*dAdy)/(-1 + x)^2 + rh*(1 + x)*B*g*siny^2*dWdy) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-4*rh^2*(1 + x)*g*dAdx + rh*(-1 + x)*(1 + x)*B*g*siny^2*(2*W + (-1 + x)*dWdx)) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-2*k*rh*(1 + x)^2*f*g^(1/2)*siny*(dpdy*dBdx - dBdy*dpdx) + (2*((1 + x)*(rh - rh*x)^2*(dAdy*dgdy + (-1 + x)^2*dAdx*dgdx) + (rh*g*(-4*rh*(-1 + x)^4*dAdx + (1 + x)*(4*rh*(-1 + x)^2*cosy/siny*dAdy - (-1 + x)^4*siny^2*dBdy*dWdy + 4*rh*(-1 + x)^2*d2Ady2 - 8*rh*(-1 + x)^3*dAdx - 2*(-1 + x)^5*siny^2*W*dBdx - (-1 + x)^6*siny^2*dBdx*dWdx - 4*rh*(1 - x)^3*(2*dAdx + (-1 + x)*d2Adx2))))/2))/(-1 + x)^4 + (rh*B*(8*k*(1 + x)^2*cosy*f*g^(1/2)*dpdx - siny*((1 + x)*siny*(dgdy*dWdy + (-1 + x)*dgdx*(2*W + (-1 + x)*dWdx)) + 2*g*(4*siny*W + 3*(1 + x)*cosy*dWdy + siny*((1 + x)*d2Wdy2 + (-1 + x)*((5 + 3*x)*dWdx + (-1 + x^2)*d2Wdx2))))))/2)

        J[v, idx+1] = Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(-4*rh^2*(-1 + x)^2*(1 + x)*g^2*W*dAdy - rh*(-1 + x)^2*(1 + x)^3*f^2*g*dBdy + (-1 + x)^2*B*(-2*rh*(1 + x)^3*cosy/siny*f^2*g + rh*(-1 + x)^2*(1 + x)*g^2*siny^2*W*dWdy)) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-4*rh^2*(-1 + x)^4*(1 + x)*g^2*W*dAdx - rh*(-1 + x)^4*(1 + x)^3*f^2*g*dBdx + rh*(-1 + x)^5*(1 + x)*B*g^2*siny^2*W*(2*W + (-1 + x)*dWdx)) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-2*(-1 + x)^2*(1 + x)^2*f*g*(rh*(1 + x)*(dBdy*dfdy + (-1 + x)^2*dBdx*dfdx) - k*g^(1/2)*siny*(rh*dpdy*(4*rh*1/siny^2*dAdx - (-1 + x)^2*W*dBdx) + rh*(-4*rh*1/siny^2*dAdy + (-1 + x)^2*W*dBdy)*dpdx)) + 2*g*((1 + x)*(rh - rh*x)^2*W*(dAdy*dgdy + (-1 + x)^2*dAdx*dgdx) + g*(-(rh*(-1 + x)^5*(1 + x)*siny^2*W^2*dBdx) + 2*(1 + x)*(rh - rh*x)^2*(dAdy*dWdy + (-1 + x)^2*dAdx*dWdx) + (rh*W*(-4*rh*(-1 + x)^4*dAdx + (1 + x)*(4*rh*(-1 + x)^2*cosy/siny*dAdy - (-1 + x)^4*siny^2*dBdy*dWdy + 4*rh*(-1 + x)^2*d2Ady2 - (-1 + x)^6*siny^2*dBdx*dWdx - 4*rh*(1 - x)^3*(2*dAdx + (-1 + x)*d2Adx2))))/2)) + 3*rh*(1 + x)^2*f^2*(((-1 + x)^2*(1 + x)*(dBdy*dgdy + (-1 + x)^2*dBdx*dgdx))/2 - g*((-1 + x)^4*dBdx + (-1 + x)^2*(1 + x)*(3*cosy/siny*dBdy + d2Bdy2 + (-1 + x)*(2*dBdx + (-1 + x)*d2Bdx2)))) + (-1 + x)^2*B*(3*rh*(1 + x)^3*f^2*(2*g + cosy/siny*dgdy) - 2*(1 + x)^2*f*g*(2*rh*(1 + x)*cosy/siny*dfdy - k*rh*g^(1/2)*(2*(1 - x)*W*(siny*dpdy - (-1 + x)*cosy*dpdx) + (-1 + x)^2*siny*(dWdy*dpdx - dpdy*dWdx))) - (rh*(-1 + x)^2*g*siny*((1 + x)*siny*W*(dgdy*dWdy + (-1 + x)*dgdx*(2*W + (-1 + x)*dWdx)) + 2*g*(4*(2 + x)*siny*W^2 - (-1 + x)^2*siny*W*dWdx + (1 + x)*siny*(dWdy^2 + (-1 + x)^2*dWdx^2) + (1 + x)*W*(3*cosy*dWdy + siny*(d2Wdy2 + (-1 + x)*(8*dWdx + (-1 + x)*d2Wdx2))))))/2))

        J[v, idx+2] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-1/16*(k*(-1 + x)^4*(4*rh*siny*(-(dBdy*dAdx) + dAdy*dBdx) + (-1 + x)*B^2*siny*2*siny*cosy*(2*W + (-1 + x)*dWdx) + B*(2*(-1 + x)*siny^3*W*dBdy - 8*rh*cosy*dAdx - (-1 + x)^2*siny^3*(dWdy*dBdx - dBdy*dWdx))))/(rh^4*(1 + x)*g^(3/2)*h) + (((-1 + x)^2*(1 + x)*(dgdy*dpdy + (-1 + x)^2*dgdx*dpdx))/2 + g*((-1 + x)^4*dpdx + (-1 + x)^2*(1 + x)*(cosy/siny*dpdy + d2pdy2 + (-1 + x)^2*d2pdx2)))/(4*rh^2*(1 + x)*g^2*h))

        J[v, idx+3] = -2*rh^4*(-1 + x)^2*(-1 + x^2)^2*f*g*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type] - 2*rh^4*(-1 + x^2)^2*f*g*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type] + Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(4*rh^4*g*(-1/2*((-1 + x^2)^2*cosy/siny*f) + (-1 + x^2)^2*dfdy) - rh^2*(1 + x)^2*(rh - rh*x)^2*f*dgdy) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(4*rh^4*g*(-1/2*((-1 + x)^4*(1 + x)*f) + (-1 + x)^2*(-1 + x^2)^2*dfdx) - rh^2*(-1 + x)^2*(1 + x)^2*(rh - rh*x)^2*f*dgdx) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-32*rh^6*(1 + x)^2*g^2*h*d0V + 2*(1 + x)*(rh - rh*x)^2*f*(4*(-1 + x)^2*(1 + x)*B^2*cosy^2*f + 2*(-1 + x)^2*(1 + x)*B*f*2*siny*cosy*dBdy + (-1 + x)^2*(1 + x)*f*siny^2*(dBdy^2 + (-1 + x)^2*dBdx^2) - rh^2*(-1 + x)^2*dgdx) + 2*(1 + x)*(rh - rh*x)^2*(2*(-1 + x)^2*(1 + x)*B^2*cosy^2*f^2 + (-1 + x)^2*(1 + x)*B*f^2*2*siny*cosy*dBdy + ((-1 + x)^2*(1 + x)*f^2*siny^2*(dBdy^2 + (-1 + x)^2*dBdx^2))/2 - rh^2*(-1 + x)^2*f*dgdx - (rh^2*(1 + x)*(dfdy*dgdy + (-1 + x)^2*dfdx*dgdx))/2) + rh^2*g*((-1 + x)^6*B^2*siny^4*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) - 8*rh*(-1 + x)^4*B*siny^2*(dAdy*dWdy + (-1 + x)*dAdx*(2*W + (-1 + x)*dWdx)) + 4*rh^2*(4*(-1 + x)^2*dAdy^2 + 4*(-1 + x)^4*dAdx^2 - ((-1 + x)^4*(1 + x)*dfdx)/2 - ((-1 + x^2)^2*(cosy/siny*dfdy + d2fdy2 + (-1 + x)^2*d2fdx2))/2)))

        J[v, idx+4] = (8*rh^3*(1 + x)*g*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*siny*dWdy)/(-1 + x)^2 + (8*rh^3*(1 + x)*g*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny*(2*W + (-1 + x)*dWdx))/(-1 + x) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(8*rh*(1 + x)*f*(-2*B^2*cosy*siny^2*dWdy + (4*rh*siny*(dAdy*dBdy + (-1 + x)^2*dAdx*dBdx))/(-1 + x)^2 + B*((8*rh*cosy*dAdy)/(-1 + x)^2 - siny^3*(dBdy*dWdy + (-1 + x)*dBdx*(2*W + (-1 + x)*dWdx)))) - (2*rh^3*(3*(1 + x)*siny*(dgdy*dWdy + (-1 + x)*dgdx*(2*W + (-1 + x)*dWdx)) + 2*g*(-4*x*siny*W + 3*(1 + x)*cosy*dWdy + siny*((1 + x)*d2Wdy2 + (-1 + x)*((3 + x)*dWdx + (-1 + x^2)*d2Wdx2)))))/(-1 + x)^2)

        J[v, idx+5] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-2*(-1 + x)*g^2 - ((1 + x)*(dgdy^2 + (-1 + x)^2*dgdx^2))/2 + g*(2*(1 + x)*cosy/siny*dgdy + (1 + x)*d2gdy2 + (-1 + x)*((-3 + x)*dgdx + (-1 + x^2)*d2gdx2)))

        if j==2
            J[v, idx+6] = 0
        else
            J[v, idx+6] = 2*rh^2*(-1 + x)^2*(1 + x)^2*g^2*h^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*dfdy + rh^2*(-1 + x)^2*g^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(4*(-1 + x)^2*(1 + x)*f*h^2 + 2*(-1 + x)^2*(1 + x)^2*h^2*dfdx) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-32*rh^4*(1 + x)^2*g^3*h^3*d0V - 2*(-1 + x)^2*(1 + x)*f^2*g*h^2*(2*(-1 + x)^2*(1 + x)*B^2*cosy^2 + (-1 + x)^2*(1 + x)*B*2*siny*cosy*dBdy + ((-1 + x)^2*(1 + x)*siny^2*(dBdy^2 + (-1 + x)^2*dBdx^2))/2) - 2*rh^2*(-1 + x)^2*(1 + x)^2*f*h^2*(dgdy^2 + (-1 + x)^2*dgdx^2) - 4*(-1 + x)^2*(1 + x)*f*g*h^2*(2*(-1 + x)^2*(1 + x)*B^2*cosy^2*f + (-1 + x)^2*(1 + x)*B*f*2*siny*cosy*dBdy + ((-1 + x)^2*(1 + x)*f*siny^2*(dBdy^2 + (-1 + x)^2*dBdx^2))/2 + 2*rh^2*(-1 + x)^2*dgdx - 2*rh^2*(1 + x)*(-(cosy/siny*dgdy) + (-1 + x)*dgdx)) - g^2*((-1 + x)^6*B^2*h^2*siny^4*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) - 8*rh*(-1 + x)^4*B*h^2*siny^2*(dAdy*dWdy + (-1 + x)*dAdx*(2*W + (-1 + x)*dWdx)) - rh^2*(-1 + x)^2*(-16*h^2*(dAdy^2 - ((-1 + x)^2*(-4*dAdx^2 + (1 + x)*dfdx))/4) - 4*(1 + x)*f*((1 + x)*(dhdy^2 + (-1 + x)^2*dhdx^2) + 4*h^2*(1 - x - ((1 + x)*(dpdy^2 + (-1 + x)^2*dpdx^2))/2) - (1 + x)*h*(d2hdy2 + (-1 + x)*(dhdx + (-1 + x)*d2hdx2))))))

        end
    end
    if funcidx==2
        J[v, idx+0] = Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*((2*(1 + x)*(rh - rh*x)^2*f*dAdy)/(-1 + x)^4 - (rh*(1 + x)*B*f*siny^2*dWdy)/2) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*((2*(1 + x)*(rh - rh*x)^2*f*dAdx)/(-1 + x)^2 - (rh*(-1 + x)*(1 + x)*B*f*siny^2*(2*W + (-1 + x)*dWdx))/2) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*((-4*rh^2*(1 + x)*(dAdy*dfdy + (-1 + x)^2*dAdx*dfdx))/(-1 + x)^2 - (k*rh*(1 + x)^2*f^2*siny*(dpdy*dBdx - dBdy*dpdx))/(2*g^(1/2)) + (rh*f*(-4*rh*(-1 + x)^4*dAdx + (1 + x)*(4*rh*(-1 + x)^2*cosy/siny*dAdy - (-1 + x)^4*siny^2*dBdy*dWdy + 4*rh*(-1 + x)^2*d2Ady2 - 8*rh*(-1 + x)^3*dAdx - 2*(-1 + x)^5*siny^2*W*dBdx - (-1 + x)^6*siny^2*dBdx*dWdx - 4*rh*(1 - x)^3*(2*dAdx + (-1 + x)*d2Adx2))))/(-1 + x)^4 + (rh*B*((2*k*(1 + x)^2*cosy*f^2*dpdx)/g^(1/2) + 2*(1 + x)*siny^2*(dfdy*dWdy + (-1 + x)*dfdx*(2*W + (-1 + x)*dWdx)) - 2*f*siny*(4*siny*W + 3*(1 + x)*cosy*dWdy + siny*((1 + x)*d2Wdy2 + (-1 + x)*((5 + 3*x)*dWdx + (-1 + x^2)*d2Wdx2)))))/2)

        J[v, idx+1] = Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(2*(1 + x)*(rh - rh*x)^2*f*g*W*dAdy + (rh*(-1 + x)^2*(1 + x)^3*f^3*dBdy)/2 + (-1 + x)^2*B*(rh*(1 + x)^3*cosy/siny*f^3 - (rh*(-1 + x)^2*(1 + x)*f*g*siny^2*W*dWdy)/2)) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(2*(-1 + x)^2*(1 + x)*(rh - rh*x)^2*f*g*W*dAdx + (rh*(-1 + x)^4*(1 + x)^3*f^3*dBdx)/2 - (rh*(-1 + x)^5*(1 + x)*B*f*g*siny^2*W*(2*W + (-1 + x)*dWdx))/2) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-8*rh^2*(-1 + x)^2*(1 + x)*g*W*(dAdy*dfdy + (-1 + x)^2*dAdx*dfdx) + (k*(-1 + x)^2*(1 + x)^2*f^2*g^(1/2)*siny*(rh*dpdy*(4*rh*1/siny^2*dAdx - (-1 + x)^2*W*dBdx) + rh*(-4*rh*1/siny^2*dAdy + (-1 + x)^2*W*dBdy)*dpdx))/2 - (-1 + x)^2*(1 + x)^2*f^2*(rh*(1 + x)*(dBdy*dfdy + (-1 + x)^2*dBdx*dfdx) - k*g^(1/2)*siny*(rh*dpdy*(4*rh*1/siny^2*dAdx - (-1 + x)^2*W*dBdx) + rh*(-4*rh*1/siny^2*dAdy + (-1 + x)^2*W*dBdy)*dpdx)) + 2*f*g*(-(rh*(-1 + x)^5*(1 + x)*siny^2*W^2*dBdx) + 2*(1 + x)*(rh - rh*x)^2*(dAdy*dWdy + (-1 + x)^2*dAdx*dWdx) + (rh*W*(-4*rh*(-1 + x)^4*dAdx + (1 + x)*(4*rh*(-1 + x)^2*cosy/siny*dAdy - (-1 + x)^4*siny^2*dBdy*dWdy + 4*rh*(-1 + x)^2*d2Ady2 - (-1 + x)^6*siny^2*dBdx*dWdx - 4*rh*(1 - x)^3*(2*dAdx + (-1 + x)*d2Adx2))))/2) + 2*f*((1 + x)*(rh - rh*x)^2*W*(dAdy*dgdy + (-1 + x)^2*dAdx*dgdx) + g*(-(rh*(-1 + x)^5*(1 + x)*siny^2*W^2*dBdx) + 2*(1 + x)*(rh - rh*x)^2*(dAdy*dWdy + (-1 + x)^2*dAdx*dWdx) + (rh*W*(-4*rh*(-1 + x)^4*dAdx + (1 + x)*(4*rh*(-1 + x)^2*cosy/siny*dAdy - (-1 + x)^4*siny^2*dBdy*dWdy + 4*rh*(-1 + x)^2*d2Ady2 - (-1 + x)^6*siny^2*dBdx*dWdx - 4*rh*(1 - x)^3*(2*dAdx + (-1 + x)*d2Adx2))))/2)) - rh*(1 + x)^2*f^3*((-1 + x)^4*dBdx + (-1 + x)^2*(1 + x)*(3*cosy/siny*dBdy + d2Bdy2 + (-1 + x)*(2*dBdx + (-1 + x)*d2Bdx2))) + (-1 + x)^2*B*(2*rh*(1 + x)^3*f^3 + 2*rh*(-1 + x)^2*(1 + x)*g*siny^2*W*(dfdy*dWdy + (-1 + x)*dfdx*(2*W + (-1 + x)*dWdx)) + (k*rh*(1 + x)^2*f^2*g^(1/2)*(2*(1 - x)*W*(siny*dpdy - (-1 + x)*cosy*dpdx) + (-1 + x)^2*siny*(dWdy*dpdx - dpdy*dWdx)))/2 - (1 + x)^2*f^2*(2*rh*(1 + x)*cosy/siny*dfdy - k*rh*g^(1/2)*(2*(1 - x)*W*(siny*dpdy - (-1 + x)*cosy*dpdx) + (-1 + x)^2*siny*(dWdy*dpdx - dpdy*dWdx))) - rh*(-1 + x)^2*f*g*siny*(4*(2 + x)*siny*W^2 - (-1 + x)^2*siny*W*dWdx + (1 + x)*siny*(dWdy^2 + (-1 + x)^2*dWdx^2) + (1 + x)*W*(3*cosy*dWdy + siny*(d2Wdy2 + (-1 + x)*(8*dWdx + (-1 + x)*d2Wdx2)))) - (rh*(-1 + x)^2*f*siny*((1 + x)*siny*W*(dgdy*dWdy + (-1 + x)*dgdx*(2*W + (-1 + x)*dWdx)) + 2*g*(4*(2 + x)*siny*W^2 - (-1 + x)^2*siny*W*dWdx + (1 + x)*siny*(dWdy^2 + (-1 + x)^2*dWdx^2) + (1 + x)*W*(3*cosy*dWdy + siny*(d2Wdy2 + (-1 + x)*(8*dWdx + (-1 + x)*d2Wdx2))))))/2))

        J[v, idx+2] = ((-1 + x)^2*f*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*dpdy)/(8*rh^2*g^2*h) + ((-1 + x)^4*f*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*dpdx)/(8*rh^2*g^2*h) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*((3*k*(-1 + x)^4*f*(4*rh*siny*(-(dBdy*dAdx) + dAdy*dBdx) + (-1 + x)*B^2*siny*2*siny*cosy*(2*W + (-1 + x)*dWdx) + B*(2*(-1 + x)*siny^3*W*dBdy - 8*rh*cosy*dAdx - (-1 + x)^2*siny^3*(dWdy*dBdx - dBdy*dWdx))))/(32*rh^4*(1 + x)*g^(5/2)*h) + (f*((-1 + x)^4*dpdx + (-1 + x)^2*(1 + x)*(cosy/siny*dpdy + d2pdy2 + (-1 + x)^2*d2pdx2)))/(4*rh^2*(1 + x)*g^2*h) - (f*(((-1 + x)^2*(1 + x)*(dgdy*dpdy + (-1 + x)^2*dgdx*dpdx))/2 + g*((-1 + x)^4*dpdx + (-1 + x)^2*(1 + x)*(cosy/siny*dpdy + d2pdy2 + (-1 + x)^2*d2pdx2))))/(2*rh^2*(1 + x)*g^3*h))

        J[v, idx+3] = -(rh^2*(1 + x)^2*(rh - rh*x)^2*f*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*dfdy) + 2*(1 + x)*(rh - rh*x)^2*f*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-(rh^2*(-1 + x)^2*f) - (rh^2*(-1 + x)^2*(1 + x)*dfdx)/2) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-4*rh^4*g*(16*rh^2*(1 + x)^2*f*h*d0V - 4*(-1 + x)^4*siny^2*W^2 - 4*(-1 + x)^5*siny^2*W*dWdx - (-1 + x)^4*siny^2*(dWdy^2 + (-1 + x)^2*dWdx^2)) + rh^2*((-1 + x)^6*B^2*f*siny^4*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) - 8*rh*(-1 + x)^4*B*f*siny^2*(dAdy*dWdy + (-1 + x)*dAdx*(2*W + (-1 + x)*dWdx)) + 4*rh^2*(((-1 + x^2)^2*(dfdy^2 + (-1 + x)^2*dfdx^2))/2 + f*(4*(-1 + x)^2*dAdy^2 + 4*(-1 + x)^4*dAdx^2 - ((-1 + x)^4*(1 + x)*dfdx)/2 - ((-1 + x^2)^2*(cosy/siny*dfdy + d2fdy2 + (-1 + x)^2*d2fdx2))/2))))

        J[v, idx+4] = (-6*rh^3*(1 + x)*f*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*siny*dWdy)/(-1 + x)^2 - (6*rh^3*(1 + x)*f*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny*(2*W + (-1 + x)*dWdx))/(-1 + x) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*((8*rh^3*(1 + x)*siny*(dfdy*dWdy + (-1 + x)*dfdx*(2*W + (-1 + x)*dWdx)))/(-1 + x)^2 - (4*rh^3*f*(-4*x*siny*W + 3*(1 + x)*cosy*dWdy + siny*((1 + x)*d2Wdy2 + (-1 + x)*((3 + x)*dWdx + (-1 + x^2)*d2Wdx2))))/(-1 + x)^2)

        J[v, idx+5] = (-1 + x)*(-1 + x^2)*f*g*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type] + (1 + x)*f*g*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type] + Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(2*(1 + x)*cosy/siny*f*g - (1 + x)*f*dgdy) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*((-3 + x)*(-1 + x)*f*g - (-1 + x)^2*(1 + x)*f*dgdx) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-4*(-1 + x)*f*g + (96*rh^2*(1 + x)*g^2*h*d0V)/(-1 + x)^2 + f*(2*(1 + x)*cosy/siny*dgdy + (1 + x)*d2gdy2 + (-1 + x)*((-3 + x)*dgdx + (-1 + x^2)*d2gdx2)))

        if j==2
            J[v, idx+6] = 0
        else
            J[v, idx+6] = Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(-4*rh^2*(-1 + x)^2*(1 + x)^2*cosy/siny*f^2*g*h^2 - 2*rh^2*(-1 + x)^2*(1 + x)^2*f^2*h^2*dgdy) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-2*(-1 + x)^2*(1 + x)*(2*rh^2*(-1 + x)^2 - 2*rh^2*(-1 + x)*(1 + x))*f^2*g*h^2 - 2*rh^2*(-1 + x)^4*(1 + x)^2*f^2*h^2*dgdx) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-2*(-1 + x)^2*(1 + x)*f^2*h^2*(2*(-1 + x)^2*(1 + x)*B^2*cosy^2*f + (-1 + x)^2*(1 + x)*B*f*2*siny*cosy*dBdy + ((-1 + x)^2*(1 + x)*f*siny^2*(dBdy^2 + (-1 + x)^2*dBdx^2))/2 + 2*rh^2*(-1 + x)^2*dgdx - 2*rh^2*(1 + x)*(-(cosy/siny*dgdy) + (-1 + x)*dgdx)) - 3*rh^2*g^2*h^2*(32*rh^2*(1 + x)^2*f*h*d0V + 12*(-1 + x)^4*siny^2*W^2 + 12*(-1 + x)^5*siny^2*W*dWdx + 3*(-1 + x)^4*siny^2*(dWdy^2 + (-1 + x)^2*dWdx^2)) - 2*g*((-1 + x)^6*B^2*f*h^2*siny^4*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) - 8*rh*(-1 + x)^4*B*f*h^2*siny^2*(dAdy*dWdy + (-1 + x)*dAdx*(2*W + (-1 + x)*dWdx)) - rh^2*(-1 + x)^2*((1 + x)^2*h^2*(dfdy^2 + (-1 + x)^2*dfdx^2) - 16*f*h^2*(dAdy^2 - ((-1 + x)^2*(-4*dAdx^2 + (1 + x)*dfdx))/4) - 2*(1 + x)*f^2*((1 + x)*(dhdy^2 + (-1 + x)^2*dhdx^2) + 4*h^2*(1 - x - ((1 + x)*(dpdy^2 + (-1 + x)^2*dpdx^2))/2) - (1 + x)*h*(d2hdy2 + (-1 + x)*(dhdx + (-1 + x)*d2hdx2))))))

        end
    end
    if funcidx==3
        J[v, idx+0] = 0

        J[v, idx+1] = 0

        J[v, idx+2] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*((k*(-1 + x)^4*f*(4*rh*siny*(-(dBdy*dAdx) + dAdy*dBdx) + (-1 + x)*B^2*siny*2*siny*cosy*(2*W + (-1 + x)*dWdx) + B*(2*(-1 + x)*siny^3*W*dBdy - 8*rh*cosy*dAdx - (-1 + x)^2*siny^3*(dWdy*dBdx - dBdy*dWdx))))/(16*rh^4*(1 + x)*g^(3/2)*h^2) - (f*(((-1 + x)^2*(1 + x)*(dgdy*dpdy + (-1 + x)^2*dgdx*dpdx))/2 + g*((-1 + x)^4*dpdx + (-1 + x)^2*(1 + x)*(cosy/siny*dpdy + d2pdy2 + (-1 + x)^2*d2pdx2))))/(4*rh^2*(1 + x)*g^2*h^2))

        J[v, idx+3] = -32*rh^6*(1 + x)^2*f*g^2*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*d0V

        J[v, idx+4] = 0

        J[v, idx+5] = (32*rh^2*(1 + x)*g^3*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*d0V)/(-1 + x)^2

        if j==2
            J[v, idx+6] = Mx[1, 1 + ordx, i]*My[1, 1 + ordy, 1, 1]
        else
            J[v, idx+6] = 2*rh^2*(-1 + x)^4*(1 + x)^2*f^2*g^2*h*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type] + 2*rh^2*(-1 + x)^2*(1 + x)^2*f^2*g^2*h*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type] - 4*rh^2*(-1 + x)^2*(1 + x)^2*f^2*g^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*dhdy - 2*rh^2*(-1 + x)^2*(1 + x)*f^2*g^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-((-1 + x)*(1 + x)*h) + 2*(-1 + x)^2*(1 + x)*dhdx) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-32*rh^4*(1 + x)^2*f*g^3*h^2*d0V - 2*rh^2*(-1 + x)^2*(1 + x)^2*f^2*h*(dgdy^2 + (-1 + x)^2*dgdx^2) - 4*(-1 + x)^2*(1 + x)*f^2*g*h*(2*(-1 + x)^2*(1 + x)*B^2*cosy^2*f + (-1 + x)^2*(1 + x)*B*f*2*siny*cosy*dBdy + ((-1 + x)^2*(1 + x)*f*siny^2*(dBdy^2 + (-1 + x)^2*dBdx^2))/2 + 2*rh^2*(-1 + x)^2*dgdx - 2*rh^2*(1 + x)*(-(cosy/siny*dgdy) + (-1 + x)*dgdx)) - 2*rh^2*g^3*h*(32*rh^2*(1 + x)^2*f*h*d0V + 12*(-1 + x)^4*siny^2*W^2 + 12*(-1 + x)^5*siny^2*W*dWdx + 3*(-1 + x)^4*siny^2*(dWdy^2 + (-1 + x)^2*dWdx^2)) - g^2*(2*(-1 + x)^6*B^2*f*h*siny^4*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) - 16*rh*(-1 + x)^4*B*f*h*siny^2*(dAdy*dWdy + (-1 + x)*dAdx*(2*W + (-1 + x)*dWdx)) - rh^2*(-1 + x)^2*(2*(1 + x)^2*h*(dfdy^2 + (-1 + x)^2*dfdx^2) - 32*f*h*(dAdy^2 - ((-1 + x)^2*(-4*dAdx^2 + (1 + x)*dfdx))/4) - 2*(1 + x)*f^2*(8*h*(1 - x - ((1 + x)*(dpdy^2 + (-1 + x)^2*dpdx^2))/2) - (1 + x)*(d2hdy2 + (-1 + x)*(dhdx + (-1 + x)*d2hdx2))))))

        end
    end
    if funcidx==4
        J[v, idx+0] = -(rh*(-1 + x)*(-1 + x^2)*B*f*g*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny^2) - rh*(1 + x)*B*f*g*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type]*siny^2 + Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(-(rh*(1 + x)*f*g*siny^2*dBdy) + (rh*B*(2*(1 + x)*g*siny^2*dfdy - f*siny*(6*(1 + x)*cosy*g + (1 + x)*siny*dgdy)))/2) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-2*rh*(-1 + x)*(1 + x)*f*g*siny^2*dBdx + (rh*B*(4*(-1 + x)*(1 + x)*g*siny^2*dfdx - f*siny*(8*g*siny + 2*(-1 + x)*(1 + x)*siny*dgdx)))/2) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-(rh*(-1 + x)^2*(1 + x)*f*g*siny^2*dBdx) + (rh*B*(2*(-1 + x)^2*(1 + x)*g*siny^2*dfdx - f*siny*(2*(-1 + x)*(5 + 3*x)*g*siny + (-1 + x)^2*(1 + x)*siny*dgdx)))/2)

        J[v, idx+1] = -(rh*(-1 + x)^6*(1 + x)*B*f*g^2*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny^2*W) - rh*(-1 + x)^4*(1 + x)*B*f*g^2*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type]*siny^2*W + Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(2*f*g^2*(2*(1 + x)*(rh - rh*x)^2*dAdy - (rh*(-1 + x)^4*(1 + x)*siny^2*W*dBdy)/2) + (-1 + x)^2*B*(rh*(-1 + x)^2*(1 + x)*g^2*siny^2*W*dfdy - (rh*(-1 + x)^2*f*g*siny*((1 + x)*siny*W*dgdy + 2*g*(3*(1 + x)*cosy*W + 2*(1 + x)*siny*dWdy)))/2 + k*rh*(-1 + x)^2*(1 + x)^2*f^2*g^(3/2)*siny*dpdx)) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(2*f*g^2*(2*(-1 + x)^2*(1 + x)*(rh - rh*x)^2*dAdx - (rh*(-1 + x)^6*(1 + x)*siny^2*W*dBdx)/2) + (-1 + x)^2*B*(-(k*rh*(-1 + x)^2*(1 + x)^2*f^2*g^(3/2)*siny*dpdy) + rh*(-1 + x)^4*(1 + x)*g^2*siny^2*W*dfdx - (rh*(-1 + x)^2*f*g*siny*((-1 + x)^2*(1 + x)*siny*W*dgdx + 2*g*(-((-1 + x)^2*siny*W) + 8*(-1 + x)*(1 + x)*siny*W + 2*(-1 + x)^2*(1 + x)*siny*dWdx)))/2)) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-4*rh^2*(-1 + x)^2*(1 + x)*g^2*(dAdy*dfdy + (-1 + x)^2*dAdx*dfdx) + k*(-1 + x)^2*(1 + x)^2*f^2*g^(3/2)*siny*(-(rh*(-1 + x)^2*dpdy*dBdx) + rh*(-1 + x)^2*dBdy*dpdx) + 2*f*g*((1 + x)*(rh - rh*x)^2*(dAdy*dgdy + (-1 + x)^2*dAdx*dgdx) + g*(-2*rh*(-1 + x)^5*(1 + x)*siny^2*W*dBdx + (rh*(-4*rh*(-1 + x)^4*dAdx + (1 + x)*(4*rh*(-1 + x)^2*cosy/siny*dAdy - (-1 + x)^4*siny^2*dBdy*dWdy + 4*rh*(-1 + x)^2*d2Ady2 - (-1 + x)^6*siny^2*dBdx*dWdx - 4*rh*(1 - x)^3*(2*dAdx + (-1 + x)*d2Adx2))))/2)) + (-1 + x)^2*B*(2*rh*(-1 + x)^3*(1 + x)*g^2*siny^2*W*dfdx + 2*k*rh*(1 - x)*(1 + x)^2*f^2*g^(3/2)*(siny*dpdy - (-1 + x)*cosy*dpdx) + rh*(-1 + x)^2*(1 + x)*g^2*siny^2*(dfdy*dWdy + (-1 + x)*dfdx*(2*W + (-1 + x)*dWdx)) - (rh*(-1 + x)^2*f*g*siny*(2*(-1 + x)*(1 + x)*siny*W*dgdx + (1 + x)*siny*(dgdy*dWdy + (-1 + x)*dgdx*(2*W + (-1 + x)*dWdx)) + 2*g*(8*(2 + x)*siny*W - (-1 + x)^2*siny*dWdx + (1 + x)*(3*cosy*dWdy + siny*(d2Wdy2 + (-1 + x)*(8*dWdx + (-1 + x)*d2Wdx2))))))/2))

        J[v, idx+2] = -1/16*(k*(-1 + x)^4*f*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(2*(-1 + x)*B^2*siny*2*siny*cosy + 2*(-1 + x)*B*siny^3*dBdy))/(rh^4*(1 + x)*g^(3/2)*h) - (k*(-1 + x)^4*f*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*((-1 + x)^2*B^2*siny*2*siny*cosy + (-1 + x)^2*B*siny^3*dBdy))/(16*rh^4*(1 + x)*g^(3/2)*h) + (k*(-1 + x)^6*B*f*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*siny^3*dBdx)/(16*rh^4*(1 + x)*g^(3/2)*h)

        J[v, idx+3] = Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(4*rh^4*(-1 + x)^4*g^2*siny^2*dWdy + rh^2*g*(-8*rh*(-1 + x)^4*B*f*siny^2*dAdy + 2*(-1 + x)^6*B^2*f*siny^4*dWdy)) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-2*rh^4*g^2*(-8*(-1 + x)^4*siny^2*W - 4*(-1 + x)^5*siny^2*dWdx) + rh^2*g*(-16*rh*(-1 + x)^5*B*f*siny^2*dAdx + (-1 + x)^6*B^2*f*siny^4*(8*W + 4*(-1 + x)*dWdx))) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-2*rh^4*g^2*(-4*(-1 + x)^5*siny^2*W - 2*(-1 + x)^6*siny^2*dWdx) + rh^2*g*(-8*rh*(-1 + x)^6*B*f*siny^2*dAdx + (-1 + x)^6*B^2*f*siny^4*(4*(-1 + x)*W + 2*(-1 + x)^2*dWdx)))

        J[v, idx+4] = (-4*rh^3*(-1 + x^2)*f*g*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny)/(-1 + x) - (4*rh^3*(1 + x)*f*g*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type]*siny)/(-1 + x)^2 + Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(4*rh*(1 + x)*f^2*(-2*B^2*cosy*siny^2 - B*siny^3*dBdy) + (8*rh^3*(1 + x)*g*siny*dfdy)/(-1 + x)^2 - (2*rh^3*f*(6*(1 + x)*cosy*g + 3*(1 + x)*siny*dgdy))/(-1 + x)^2) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-8*rh*(-1 + x)*(1 + x)*B*f^2*siny^3*dBdx + (16*rh^3*(1 + x)*g*siny*dfdx)/(-1 + x) - (2*rh^3*f*(-8*x*g*siny + 6*(-1 + x)*(1 + x)*siny*dgdx))/(-1 + x)^2) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-4*rh*(-1 + x)^2*(1 + x)*B*f^2*siny^3*dBdx + 8*rh^3*(1 + x)*g*siny*dfdx - (2*rh^3*f*(2*(-1 + x)*(3 + x)*g*siny + 3*(-1 + x)^2*(1 + x)*siny*dgdx))/(-1 + x)^2)

        J[v, idx+5] = 0

        if j==2
            J[v, idx+6] = 0
        else
            J[v, idx+6] = Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(-6*rh^2*(-1 + x)^4*g^3*h^2*siny^2*dWdy - g^2*(-8*rh*(-1 + x)^4*B*f*h^2*siny^2*dAdy + 2*(-1 + x)^6*B^2*f*h^2*siny^4*dWdy)) + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-(rh^2*g^3*h^2*(24*(-1 + x)^4*siny^2*W + 12*(-1 + x)^5*siny^2*dWdx)) - g^2*(-16*rh*(-1 + x)^5*B*f*h^2*siny^2*dAdx + (-1 + x)^6*B^2*f*h^2*siny^4*(8*W + 4*(-1 + x)*dWdx))) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-(rh^2*g^3*h^2*(12*(-1 + x)^5*siny^2*W + 6*(-1 + x)^6*siny^2*dWdx)) - g^2*(-8*rh*(-1 + x)^6*B*f*h^2*siny^2*dAdx + (-1 + x)^6*B^2*f*h^2*siny^4*(4*(-1 + x)*W + 2*(-1 + x)^2*dWdx)))

        end
    end
    if funcidx==5
        J[v, idx+0] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(2*k*rh*(1 + x)^2*B*cosy*f^2*g^(1/2) + k*rh*(1 + x)^2*f^2*g^(1/2)*siny*dBdy) - k*rh*(1 + x)^2*f^2*g^(1/2)*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*siny*dBdx

        J[v, idx+1] = Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(k*rh*(-1 + x)^2*(1 + x)^2*f^2*g^(3/2)*siny*(-4*rh*1/siny^2*dAdy + (-1 + x)^2*W*dBdy) + k*rh*(-1 + x)^2*(1 + x)^2*B*f^2*g^(3/2)*(-2*(1 - x)*(-1 + x)*cosy*W + (-1 + x)^2*siny*dWdy)) + Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(k*rh*(-1 + x)^2*(1 + x)^2*f^2*g^(3/2)*siny*(4*rh*1/siny^2*dAdx - (-1 + x)^2*W*dBdx) + k*rh*(-1 + x)^2*(1 + x)^2*B*f^2*g^(3/2)*(2*(1 - x)*siny*W - (-1 + x)^2*siny*dWdx))

        J[v, idx+2] = ((-1 + x)^4*f*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type])/(4*rh^2*g*h) + ((-1 + x)^2*f*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type])/(4*rh^2*g*h) - Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*d2V + (f*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*((-1 + x)^2*(1 + x)*cosy/siny*g + ((-1 + x)^2*(1 + x)*dgdy)/2))/(4*rh^2*(1 + x)*g^2*h) + (f*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*((-1 + x)^4*g + ((-1 + x)^4*(1 + x)*dgdx)/2))/(4*rh^2*(1 + x)*g^2*h)

        J[v, idx+3] = -32*rh^6*(1 + x)^2*f*g^2*h*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*d1V

        J[v, idx+4] = 0

        J[v, idx+5] = (32*rh^2*(1 + x)*g^3*h*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*d1V)/(-1 + x)^2

        if j==2
            J[v, idx+6] = 0
        else
            J[v, idx+6] = -32*rh^4*(1 + x)^2*f*g^3*h^3*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*d1V + 8*rh^2*(-1 + x)^2*(1 + x)^2*f^2*g^2*h^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*dpdy + 8*rh^2*(-1 + x)^4*(1 + x)^2*f^2*g^2*h^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*dpdx

        end
    end
    if funcidx==6
        J[v, idx+0] = (-4*rh^2*(1 - x)^3*(1 + x)*f*g*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type])/(-1 + x)^3 + (4*rh^2*(1 + x)*f*g*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type])/(-1 + x)^2 + Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*((-4*rh^2*(1 + x)*g*dfdy)/(-1 + x)^2 + (2*f*(2*rh^2*(-1 + x)^2*(1 + x)*cosy/siny*g + (1 + x)*(rh - rh*x)^2*dgdy))/(-1 + x)^4) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-4*rh^2*(1 + x)*g*dfdx + (2*f*((rh*(-4*rh*(-1 + x)^4 + (-8*rh*(1 - x)^3 - 8*rh*(-1 + x)^3)*(1 + x))*g)/2 + (-1 + x)^2*(1 + x)*(rh - rh*x)^2*dgdx))/(-1 + x)^4)

        J[v, idx+1] = -4*rh^2*(1 - x)^3*(-1 + x)*(1 + x)*f*g^2*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type]*W + 4*rh^2*(-1 + x)^2*(1 + x)*f*g^2*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type]*W + Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(-4*rh^2*(-1 + x)^2*(1 + x)*g^2*W*dfdy + 2*f*g*((1 + x)*(rh - rh*x)^2*W*dgdy + g*(2*rh^2*(-1 + x)^2*(1 + x)*cosy/siny*W + 2*(1 + x)*(rh - rh*x)^2*dWdy)) - 4*k*rh^2*(-1 + x)^2*(1 + x)^2*1/siny*f^2*g^(3/2)*dpdx) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(4*k*rh^2*(-1 + x)^2*(1 + x)^2*1/siny*f^2*g^(3/2)*dpdy - 4*rh^2*(-1 + x)^4*(1 + x)*g^2*W*dfdx + 2*f*g*((-1 + x)^2*(1 + x)*(rh - rh*x)^2*W*dgdx + g*((rh*(-4*rh*(-1 + x)^4 - 8*rh*(1 - x)^3*(1 + x))*W)/2 + 2*(-1 + x)^2*(1 + x)*(rh - rh*x)^2*dWdx)))

        J[v, idx+2] = -1/16*(k*(-1 + x)^4*f*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-8*rh*B*cosy - 4*rh*siny*dBdy))/(rh^4*(1 + x)*g^(3/2)*h) - (k*(-1 + x)^4*f*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*siny*dBdx)/(4*rh^3*(1 + x)*g^(3/2)*h)

        J[v, idx+3] = rh^2*g*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(32*rh^2*(-1 + x)^2*f*dAdy - 8*rh*(-1 + x)^4*B*f*siny^2*dWdy) + rh^2*g*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(32*rh^2*(-1 + x)^4*f*dAdx - 8*rh*(-1 + x)^5*B*f*siny^2*(2*W + (-1 + x)*dWdx))

        J[v, idx+4] = 4*rh*(1 + x)*f^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*((8*rh*B*cosy)/(-1 + x)^2 + (4*rh*siny*dBdy)/(-1 + x)^2) + 16*rh^2*(1 + x)*f^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny*dBdx

        J[v, idx+5] = 0

        if j==2
            J[v, idx+6] = 0
        else
            J[v, idx+6] = -(g^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(32*rh^2*(-1 + x)^2*f*h^2*dAdy - 8*rh*(-1 + x)^4*B*f*h^2*siny^2*dWdy)) - g^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(32*rh^2*(-1 + x)^4*f*h^2*dAdx - 8*rh*(-1 + x)^5*B*f*h^2*siny^2*(2*W + (-1 + x)*dWdx))

        end
    end
    if funcidx==7
        J[v, idx+0] = Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(-(rh*(1 + x)*f*g*siny^2*dWdy) + k*rh*(1 + x)^2*f^2*g^(1/2)*siny*dpdx) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-(k*rh*(1 + x)^2*f^2*g^(1/2)*siny*dpdy) + (rh*(1 + x)*f*g*(-2*(-1 + x)^5*siny^2*W - (-1 + x)^6*siny^2*dWdx))/(-1 + x)^4) + (rh*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(4*k*(1 + x)^2*cosy*f^2*g^(1/2)*dpdx + 2*(1 + x)*g*siny^2*(dfdy*dWdy + (-1 + x)*dfdx*(2*W + (-1 + x)*dWdx)) - f*siny*((1 + x)*siny*(dgdy*dWdy + (-1 + x)*dgdx*(2*W + (-1 + x)*dWdx)) + 2*g*(4*siny*W + 3*(1 + x)*cosy*dWdy + siny*((1 + x)*d2Wdy2 + (-1 + x)*((5 + 3*x)*dWdx + (-1 + x^2)*d2Wdx2))))))/2

        J[v, idx+1] = -(rh*(-1 + x)^4*(1 + x)^3*f^3*g*Mx[3, 1 + ordx, i]*My[1, 1 + ordy, j, type]) - rh*(-1 + x)^2*(1 + x)^3*f^3*g*Mx[1, 1 + ordx, i]*My[3, 1 + ordy, j, type] + Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(rh*(1 + x)^2*f^3*(-3*(-1 + x)^2*(1 + x)*cosy/siny*g + ((-1 + x)^2*(1 + x)*dgdy)/2) - rh*(-1 + x)^4*(1 + x)*f*g^2*siny^2*W*dWdy - (-1 + x)^2*(1 + x)^2*f^2*g*(rh*(1 + x)*dfdy - k*rh*(-1 + x)^2*g^(1/2)*siny*W*dpdx)) + Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-((-1 + x)^2*(1 + x)^2*f^2*g*(k*rh*(-1 + x)^2*g^(1/2)*siny*W*dpdy + rh*(-1 + x)^2*(1 + x)*dfdx)) + rh*(1 + x)^2*f^3*(-(((-1 + x)^4 + 2*(-1 + x)^3*(1 + x))*g) + ((-1 + x)^4*(1 + x)*dgdx)/2) + 2*f*g^2*(-(rh*(-1 + x)^5*(1 + x)*siny^2*W^2) - (rh*(-1 + x)^6*(1 + x)*siny^2*W*dWdx)/2)) + (-1 + x)^2*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(rh*(1 + x)^3*f^3*(2*g + cosy/siny*dgdy) + rh*(-1 + x)^2*(1 + x)*g^2*siny^2*W*(dfdy*dWdy + (-1 + x)*dfdx*(2*W + (-1 + x)*dWdx)) - (1 + x)^2*f^2*g*(2*rh*(1 + x)*cosy/siny*dfdy - k*rh*g^(1/2)*(2*(1 - x)*W*(siny*dpdy - (-1 + x)*cosy*dpdx) + (-1 + x)^2*siny*(dWdy*dpdx - dpdy*dWdx))) - (rh*(-1 + x)^2*f*g*siny*((1 + x)*siny*W*(dgdy*dWdy + (-1 + x)*dgdx*(2*W + (-1 + x)*dWdx)) + 2*g*(4*(2 + x)*siny*W^2 - (-1 + x)^2*siny*W*dWdx + (1 + x)*siny*(dWdy^2 + (-1 + x)^2*dWdx^2) + (1 + x)*W*(3*cosy*dWdy + siny*(d2Wdy2 + (-1 + x)*(8*dWdx + (-1 + x)*d2Wdx2))))))/2)

        J[v, idx+2] = -1/16*(k*(-1 + x)^4*f*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(4*rh*siny*dAdy - (-1 + x)^2*B*siny^3*dWdy))/(rh^4*(1 + x)*g^(3/2)*h) - (k*(-1 + x)^4*f*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*(-4*rh*siny*dAdx + B*(2*(-1 + x)*siny^3*W + (-1 + x)^2*siny^3*dWdx)))/(16*rh^4*(1 + x)*g^(3/2)*h) - (k*(-1 + x)^4*f*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(2*(-1 + x)*siny^3*W*dBdy - 8*rh*cosy*dAdx + 2*(-1 + x)*B*siny*2*siny*cosy*(2*W + (-1 + x)*dWdx) - (-1 + x)^2*siny^3*(dWdy*dBdx - dBdy*dWdx)))/(16*rh^4*(1 + x)*g^(3/2)*h)

        J[v, idx+3] = 2*(1 + x)*(rh - rh*x)^2*f*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*((-1 + x)^2*(1 + x)*B*f^2*2*siny*cosy + (-1 + x)^2*(1 + x)*f^2*siny^2*dBdy) + 2*(-1 + x)^4*(1 + x)^2*(rh - rh*x)^2*f^3*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny^2*dBdx + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(2*(1 + x)*(rh - rh*x)^2*f*(4*(-1 + x)^2*(1 + x)*B*cosy^2*f^2 + (-1 + x)^2*(1 + x)*f^2*2*siny*cosy*dBdy) + rh^2*g*(2*(-1 + x)^6*B*f*siny^4*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) - 8*rh*(-1 + x)^4*f*siny^2*(dAdy*dWdy + (-1 + x)*dAdx*(2*W + (-1 + x)*dWdx))))

        J[v, idx+4] = 4*rh*(1 + x)*f^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*((4*rh*siny*dAdy)/(-1 + x)^2 - B*siny^3*dWdy) + 4*rh*(1 + x)*f^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(4*rh*siny*dAdx - (-1 + x)*B*siny^3*(2*W + (-1 + x)*dWdx)) + 4*rh*(1 + x)*f^2*Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*((8*rh*cosy*dAdy)/(-1 + x)^2 - 4*B*cosy*siny^2*dWdy - siny^3*(dBdy*dWdy + (-1 + x)*dBdx*(2*W + (-1 + x)*dWdx)))

        J[v, idx+5] = 0

        if j==2
            J[v, idx+6] = 0
        else
            J[v, idx+6] = -2*(-1 + x)^2*(1 + x)*f^2*g*h^2*Mx[1, 1 + ordx, i]*My[2, 1 + ordy, j, type]*((-1 + x)^2*(1 + x)*B*f*2*siny*cosy + (-1 + x)^2*(1 + x)*f*siny^2*dBdy) - 2*(-1 + x)^6*(1 + x)^2*f^3*g*h^2*Mx[2, 1 + ordx, i]*My[1, 1 + ordy, j, type]*siny^2*dBdx + Mx[1, 1 + ordx, i]*My[1, 1 + ordy, j, type]*(-2*(-1 + x)^2*(1 + x)*f^2*g*h^2*(4*(-1 + x)^2*(1 + x)*B*cosy^2*f + (-1 + x)^2*(1 + x)*f*2*siny*cosy*dBdy) - g^2*(2*(-1 + x)^6*B*f*h^2*siny^4*(4*W^2 + dWdy^2 + 4*(-1 + x)*W*dWdx + (-1 + x)^2*dWdx^2) - 8*rh*(-1 + x)^4*f*h^2*siny^2*(dAdy*dWdy + (-1 + x)*dAdx*(2*W + (-1 + x)*dWdx))))

        end
    end
end
#DEFINE THE WHOLE SYSTEM (BCs+FIELD EQS) ON OUR GRID
function residual_jac!(R::Matrix{Float64},J::Matrix{Float64},a::Matrix{Float64},f::Field,g::Field,h::Field,W::Field,p::Field,A::Field,B::Field,WBC::Float64,rh::Float64,q::Float64,k::Float64,spin::Bool,x::Vector{Float64}=X,y::Vector{Float64}=Y)

    f.a=a[0*Nx+1:1*Nx,:]

    g.a=a[1*Nx+1:2*Nx,:]

    h.a=a[2*Nx+1:3*Nx,:]

    W.a=a[3*Nx+1:4*Nx,:]

    p.a=a[4*Nx+1:5*Nx,:]

    A.a=a[5*Nx+1:6*Nx,:]

    B.a=a[6*Nx+1:7*Nx,:]

    type=[f.type,g.type,h.type,W.type,p.type,A.type,B.type]
    ff=Matrix{Float64}(undef,Nx,Ny); dfdx=Matrix{Float64}(undef,Nx,Ny); d2fdx2=Matrix{Float64}(undef,Nx,Ny); dfdy=Matrix{Float64}(undef,Nx,Ny); d2fdy2=Matrix{Float64}(undef,Nx,Ny); d2fdxdy=Matrix{Float64}(undef,Nx,Ny);
    gg=Matrix{Float64}(undef,Nx,Ny); dgdx=Matrix{Float64}(undef,Nx,Ny); d2gdx2=Matrix{Float64}(undef,Nx,Ny); dgdy=Matrix{Float64}(undef,Nx,Ny); d2gdy2=Matrix{Float64}(undef,Nx,Ny); d2gdxdy=Matrix{Float64}(undef,Nx,Ny);
    hh=Matrix{Float64}(undef,Nx,Ny); dhdx=Matrix{Float64}(undef,Nx,Ny); d2hdx2=Matrix{Float64}(undef,Nx,Ny); dhdy=Matrix{Float64}(undef,Nx,Ny); d2hdy2=Matrix{Float64}(undef,Nx,Ny); d2hdxdy=Matrix{Float64}(undef,Nx,Ny);
    WW=Matrix{Float64}(undef,Nx,Ny); dWdx=Matrix{Float64}(undef,Nx,Ny); d2Wdx2=Matrix{Float64}(undef,Nx,Ny); dWdy=Matrix{Float64}(undef,Nx,Ny); d2Wdy2=Matrix{Float64}(undef,Nx,Ny); d2Wdxdy=Matrix{Float64}(undef,Nx,Ny);
    pp=Matrix{Float64}(undef,Nx,Ny); dpdx=Matrix{Float64}(undef,Nx,Ny); d2pdx2=Matrix{Float64}(undef,Nx,Ny); dpdy=Matrix{Float64}(undef,Nx,Ny); d2pdy2=Matrix{Float64}(undef,Nx,Ny); d2pdxdy=Matrix{Float64}(undef,Nx,Ny);
    AA=Matrix{Float64}(undef,Nx,Ny); dAdx=Matrix{Float64}(undef,Nx,Ny); d2Adx2=Matrix{Float64}(undef,Nx,Ny); dAdy=Matrix{Float64}(undef,Nx,Ny); d2Ady2=Matrix{Float64}(undef,Nx,Ny); d2Adxdy=Matrix{Float64}(undef,Nx,Ny);
    BB=Matrix{Float64}(undef,Nx,Ny); dBdx=Matrix{Float64}(undef,Nx,Ny); d2Bdx2=Matrix{Float64}(undef,Nx,Ny); dBdy=Matrix{Float64}(undef,Nx,Ny); d2Bdy2=Matrix{Float64}(undef,Nx,Ny); d2Bdxdy=Matrix{Float64}(undef,Nx,Ny);
    haxis=Array{Float64}(undef,Nx);
    d0V=Matrix{Float64}(undef,Nx,Ny); d1V=Matrix{Float64}(undef,Nx,Ny); d2V=Matrix{Float64}(undef,Nx,Ny); d3V=Matrix{Float64}(undef,Nx,Ny);

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
            pp[i,j]=p(i,j+1)
            dpdx[i,j]=p(i,j+1,dx=1)
            d2pdx2[i,j]=p(i,j+1,dx=2)
            dpdy[i,j]=p(i,j+1,dy=1)
            d2pdy2[i,j]=p(i,j+1,dy=2)
            d2pdxdy[i,j]=p(i,j+1,dx=1,dy=1)
            AA[i,j]=A(i,j+1)
            dAdx[i,j]=A(i,j+1,dx=1)
            d2Adx2[i,j]=A(i,j+1,dx=2)
            dAdy[i,j]=A(i,j+1,dy=1)
            d2Ady2[i,j]=A(i,j+1,dy=2)
            d2Adxdy[i,j]=A(i,j+1,dx=1,dy=1)
            BB[i,j]=B(i,j+1)
            dBdx[i,j]=B(i,j+1,dx=1)
            d2Bdx2[i,j]=B(i,j+1,dx=2)
            dBdy[i,j]=B(i,j+1,dy=1)
            d2Bdy2[i,j]=B(i,j+1,dy=2)
            d2Bdxdy[i,j]=B(i,j+1,dx=1,dy=1)
            d0V[i,j]=V(pp[i,j])
            d1V[i,j]=V(pp[i,j],1)
            d2V[i,j]=V(pp[i,j],2)
            d3V[i,j]=V(pp[i,j],3)
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
                        BC_Horizon!(i,j,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],pp[i,j],dpdx[i,j],d2pdx2[i,j],dpdy[i,j],d2pdy2[i,j],d2pdxdy[i,j],AA[i,j],dAdx[i,j],d2Adx2[i,j],dAdy[i,j],d2Ady2[i,j],d2Adxdy[i,j],BB[i,j],dBdx[i,j],d2Bdx2[i,j],dBdy[i,j],d2Bdy2[i,j],d2Bdxdy[i,j],R,idx,WBC,rh,q,k)
                    else
                        BC_Horizon_Spin!(i,j,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],pp[i,j],dpdx[i,j],d2pdx2[i,j],dpdy[i,j],d2pdy2[i,j],d2pdxdy[i,j],AA[i,j],dAdx[i,j],d2Adx2[i,j],dAdy[i,j],d2Ady2[i,j],d2Adxdy[i,j],BB[i,j],dBdx[i,j],d2Bdx2[i,j],dBdy[i,j],d2Bdy2[i,j],d2Bdxdy[i,j],R,idx,WBC,rh,q,k)
                    end
                elseif i==Nx
                    #DEFINE BCS AT INFINITY
                    if !spin
                        BC_Infinity!(i,j,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],pp[i,j],dpdx[i,j],d2pdx2[i,j],dpdy[i,j],d2pdy2[i,j],d2pdxdy[i,j],AA[i,j],dAdx[i,j],d2Adx2[i,j],dAdy[i,j],d2Ady2[i,j],d2Adxdy[i,j],BB[i,j],dBdx[i,j],d2Bdx2[i,j],dBdy[i,j],d2Bdy2[i,j],d2Bdxdy[i,j],R,idx,WBC,rh,q,k)
                    else
                        BC_Infinity_Spin!(i,j,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],pp[i,j],dpdx[i,j],d2pdx2[i,j],dpdy[i,j],d2pdy2[i,j],d2pdxdy[i,j],AA[i,j],dAdx[i,j],d2Adx2[i,j],dAdy[i,j],d2Ady2[i,j],d2Adxdy[i,j],BB[i,j],dBdx[i,j],d2Bdx2[i,j],dBdy[i,j],d2Bdy2[i,j],d2Bdxdy[i,j],R,idx,WBC,rh,q,k)
                    end
                else
                    #DEFINE FIELD EQUATIONS EVERYWHERE ELSE ON THE GRID
                    Field_Eqs!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],pp[i,j],dpdx[i,j],d2pdx2[i,j],dpdy[i,j],d2pdy2[i,j],d2pdxdy[i,j],AA[i,j],dAdx[i,j],d2Adx2[i,j],dAdy[i,j],d2Ady2[i,j],d2Adxdy[i,j],BB[i,j],dBdx[i,j],d2Bdx2[i,j],dBdy[i,j],d2Bdy2[i,j],d2Bdxdy[i,j],haxis[i],siny[j],cosy[j],d0V[i,j],d1V[i,j],d2V[i,j],d3V[i,j],R,idx,WBC,rh,q,k,X[i],Y[j])
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
                            BC_Horizon_Jac!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],pp[i,j],dpdx[i,j],d2pdx2[i,j],dpdy[i,j],d2pdy2[i,j],d2pdxdy[i,j],AA[i,j],dAdx[i,j],d2Adx2[i,j],dAdy[i,j],d2Ady2[i,j],d2Adxdy[i,j],BB[i,j],dBdx[i,j],d2Bdx2[i,j],dBdy[i,j],d2Bdy2[i,j],d2Bdxdy[i,j],J,idx,WBC,rh,q,k,v,funcidx,type[funcidx],ordx,ordy)
                        else
                            BC_Horizon_Jac_Spin!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],pp[i,j],dpdx[i,j],d2pdx2[i,j],dpdy[i,j],d2pdy2[i,j],d2pdxdy[i,j],AA[i,j],dAdx[i,j],d2Adx2[i,j],dAdy[i,j],d2Ady2[i,j],d2Adxdy[i,j],BB[i,j],dBdx[i,j],d2Bdx2[i,j],dBdy[i,j],d2Bdy2[i,j],d2Bdxdy[i,j],J,idx,WBC,rh,q,k,v,funcidx,type[funcidx],ordx,ordy)
                        end
                    elseif i==Nx
                        #DEFINE BCS AT INFINITY
                        if !spin
                            BC_Infinity_Jac!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],pp[i,j],dpdx[i,j],d2pdx2[i,j],dpdy[i,j],d2pdy2[i,j],d2pdxdy[i,j],AA[i,j],dAdx[i,j],d2Adx2[i,j],dAdy[i,j],d2Ady2[i,j],d2Adxdy[i,j],BB[i,j],dBdx[i,j],d2Bdx2[i,j],dBdy[i,j],d2Bdy2[i,j],d2Bdxdy[i,j],J,idx,WBC,rh,q,k,v,funcidx,type[funcidx],ordx,ordy)
                        else
                            BC_Infinity_Jac_Spin!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],pp[i,j],dpdx[i,j],d2pdx2[i,j],dpdy[i,j],d2pdy2[i,j],d2pdxdy[i,j],AA[i,j],dAdx[i,j],d2Adx2[i,j],dAdy[i,j],d2Ady2[i,j],d2Adxdy[i,j],BB[i,j],dBdx[i,j],d2Bdx2[i,j],dBdy[i,j],d2Bdy2[i,j],d2Bdxdy[i,j],J,idx,WBC,rh,q,k,v,funcidx,type[funcidx],ordx,ordy)
                        end
                    else
                        #DEFINE FIELD EQUATIONS EVERYWHERE ELSE ON THE GRID
                        Field_Eqs_Jac!(i,j+1,ff[i,j],dfdx[i,j],d2fdx2[i,j],dfdy[i,j],d2fdy2[i,j],d2fdxdy[i,j],gg[i,j],dgdx[i,j],d2gdx2[i,j],dgdy[i,j],d2gdy2[i,j],d2gdxdy[i,j],hh[i,j],dhdx[i,j],d2hdx2[i,j],dhdy[i,j],d2hdy2[i,j],d2hdxdy[i,j],WW[i,j],dWdx[i,j],d2Wdx2[i,j],dWdy[i,j],d2Wdy2[i,j],d2Wdxdy[i,j],pp[i,j],dpdx[i,j],d2pdx2[i,j],dpdy[i,j],d2pdy2[i,j],d2pdxdy[i,j],AA[i,j],dAdx[i,j],d2Adx2[i,j],dAdy[i,j],d2Ady2[i,j],d2Adxdy[i,j],BB[i,j],dBdx[i,j],d2Bdx2[i,j],dBdy[i,j],d2Bdy2[i,j],d2Bdxdy[i,j],haxis[i],siny[j],cosy[j],d0V[i,j],d1V[i,j],d2V[i,j],d3V[i,j],J,idx,WBC,rh,q,k,v,funcidx,type[funcidx],ordx,ordy,X[i],Y[j])
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
function solve_system(f::Field,g::Field,h::Field,W::Field,p::Field,A::Field,B::Field,tol::Float64,WBC::Float64,rh::Float64,q::Float64,k::Float64,spin::Bool=false,show_trace::Bool=true,iterations::Int64=30,method = :newton)
    return nlsolve(only_fj!((R,J,a)->residual_jac!(R,J,a,f,g,h,W,p,A,B,WBC,rh,q,k,spin)), vcat(f.a,g.a,h.a,W.a,p.a,A.a,B.a), ftol=0.0,xtol=tol, show_trace=show_trace, iterations = iterations, method = method)
end
nothing