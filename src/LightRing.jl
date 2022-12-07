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