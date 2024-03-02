function LightRing(f::Field, g::Field, h::Field, W::Field,rh::Float64,q::Float64=0.0; guess = nothing)
    y=pi/2
    M=GetMass(f,g,h,W,rh)
    j=GetJ(f,g,h,W,rh)
    chi = j/M^2

    compare_kerr = true

    if chi >= 1.0
        compare_kerr=false
    end

    #guess is of the form [x_co-rotating, x_counter-rotating] is most cases
    use_guess=false
    if !isnothing(guess)
        use_guess = true
    end


    if compare_kerr
        rhK=M/2*sqrt(1-chi^2-q^2)

        rbl_minus=M*find_zero(r->2*q^2+(-3+r)*r-2*chi*sqrt(r-q^2),2*(1+cos(2/3*acos(chi))))
        rbl_plus=M*find_zero(r->2*q^2+(-3+r)*r+2*chi*sqrt(r-q^2),2*(1+cos(2/3*acos(-chi))))

        xkerr_minus=rBLKerr2x(rbl_minus,rhK,chi,q)
        xkerr_plus=rBLKerr2x(rbl_plus,rhK,chi,q)

        wKerr_co= (rbl_plus^2*sqrt(M*(rbl_plus-M*q^2)) + M^3*q^2*chi - M^2*rbl_plus*chi) / (rbl_plus^4+M^3*(M*q^2-rbl_plus)*chi^2)
        wKerr_counter= (-rbl_minus^2*sqrt(M*(rbl_minus-M*q^2)) + M^3*q^2*chi - M^2*rbl_minus*chi) / (rbl_minus^4+M^3*(M*q^2-rbl_minus)*chi^2)

        RKerr_counter = 2*rhK*sqrt(gKerrN(xkerr_minus,y,rhK,0.0,chi,q)/((-1+xkerr_minus)^2*fKerrN(xkerr_minus,y,rhK,0.0,chi,q)))
        RKerr_co = 2*rhK*sqrt(gKerrN(xkerr_plus,y,rhK,0.0,chi,q)/((-1+xkerr_plus)^2*fKerrN(xkerr_plus,y,rhK,0.0,chi,q)))
    end

    #See arXiv:1609.01340
    hplus = x->( -gtphi(x,y,f,g,h,W,rh) + sqrt(gtphi(x,y,f,g,h,W,rh)^2 - gtt(x,y,f,g,h,W,rh)*gphiphi(x,y,f,g,h,W,rh)) )/gtt(x,y,f,g,h,W,rh)
    hminus = x->( -gtphi(x,y,f,g,h,W,rh) - sqrt(gtphi(x,y,f,g,h,W,rh)^2 - gtt(x,y,f,g,h,W,rh)*gphiphi(x,y,f,g,h,W,rh)) )/gtt(x,y,f,g,h,W,rh)

    eqminus = x->gphiphi(x,y,f,g,h,W,rh,dx=1) + 2*gtphi(x,y,f,g,h,W,rh,dx=1)*hminus(x) + gtt(x,y,f,g,h,W,rh,dx=1)*hminus(x)^2
    eqplus = x->gphiphi(x,y,f,g,h,W,rh,dx=1) + 2*gtphi(x,y,f,g,h,W,rh,dx=1)*hplus(x) + gtt(x,y,f,g,h,W,rh,dx=1)*hplus(x)^2

    hder=1e-6

    println()
    println("LIGHT RING DETAILS")

    to_return=[]
    to_returnx=[]

    #MINUS

    if !use_guess
        sols_minus=find_zeros(x->eqminus(x),-0.999,1.0)
        for sol in sols_minus
            w = 1/hminus(sol)
            d2h=( hminus(sol+hder) - 2*hminus(sol) +hminus(sol-hder) ) / (hder^2) #2nd derivative finite differences
            R = CircumferencialRadius(sol,f,g,h,W,rh)

            append!(to_return,[R/RKerr_co-1, w/wKerr_co-1])
            append!(to_returnx,[sol])

            if compare_kerr 
                if d2h > 0.0 && w > 0.0
                    println("Prograde (co-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_co-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_co-1)
                elseif d2h > 0.0 && w < 0.0
                    println("Retrograde (counter-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_counter-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_counter-1)
                elseif d2h < 0.0 && w > 0.0
                    println("Prograde (co-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_co-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_co-1)
                else
                    println("Retrograde (counter-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_counter-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_counter-1)
                end
            else
                if d2h > 0.0 && w > 0.0
                    println("Prograde (co-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
                elseif d2h > 0.0 && w < 0.0
                    println("Retrograde (counter-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
                elseif d2h < 0.0 && w > 0.0
                    println("Prograde (co-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
                else
                    println("Retrograde (counter-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
                end
            end
        end
    else
        sol = find_zero(x->eqminus(x),guess[1])
        w = 1/hminus(sol)
        d2h=( hminus(sol+hder) - 2*hminus(sol) +hminus(sol-hder) ) / (hder^2)
        R = CircumferencialRadius(sol,f,g,h,W,rh)

        append!(to_return,[R/RKerr_co-1, w/wKerr_co-1])
        append!(to_returnx,[sol])

        if compare_kerr 
            if d2h > 0.0 && w > 0.0
                println("Prograde (co-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_co-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_co-1)
                #append!(to_return,[R/RKerr_co-1, w/wKerr_co-1])
            elseif d2h > 0.0 && w < 0.0
                println("Retrograde (counter-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_counter-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_counter-1)
            elseif d2h < 0.0 && w > 0.0
                println("Prograde (co-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_co-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_co-1)
            else
                println("Retrograde (counter-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_counter-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_counter-1)
            end
        else
            if d2h > 0.0 && w > 0.0
                println("Prograde (co-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
            elseif d2h > 0.0 && w < 0.0
                println("Retrograde (counter-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
            elseif d2h < 0.0 && w > 0.0
                println("Prograde (co-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
            else
                println("Retrograde (counter-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
            end
        end
    end

    #PLUS

    if !use_guess
        sols_plus=find_zeros(x->eqplus(x),-0.999,1.0)
        for sol in sols_plus
            w = 1/hplus(sol)
            d2h=( hplus(sol+hder) - 2*hplus(sol) +hplus(sol-hder) ) / (hder^2)
            R = CircumferencialRadius(sol,f,g,h,W,rh)

            append!(to_return,[R/RKerr_counter-1, w/wKerr_counter-1])
            append!(to_returnx,[sol])

            if compare_kerr
                if d2h > 0.0 && w > 0.0
                    println("Prograde (co-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_co-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_co-1)
                elseif d2h > 0.0 && w < 0.0
                    println("Retrograde (counter-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_counter-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_counter-1)
                elseif d2h < 0.0 && w > 0.0
                    println("Prograde (co-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_co-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_co-1)
                else
                    println("Retrograde (counter-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_counter-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_counter-1)
                end
            else
                if d2h > 0.0 && w > 0.0
                    println("Prograde (co-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
                elseif d2h > 0.0 && w < 0.0
                    println("Retrograde (counter-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
                elseif d2h < 0.0 && w > 0.0
                    println("Prograde (co-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
                else
                    println("Retrograde (counter-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
                end
            end
        end
    else
        sol = find_zero(x->eqplus(x),guess[2])
        w = 1/hplus(sol)
        d2h=( hplus(sol+hder) - 2*hplus(sol) +hplus(sol-hder) ) / (hder^2)
        R = CircumferencialRadius(sol,f,g,h,W,rh)

        append!(to_return,[R/RKerr_counter-1, w/wKerr_counter-1])
        append!(to_returnx,[sol])

        if compare_kerr
            if d2h > 0.0 && w > 0.0
                println("Prograde (co-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_co-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_co-1)
            elseif d2h > 0.0 && w < 0.0
                println("Retrograde (counter-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_counter-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_counter-1)
            elseif d2h < 0.0 && w > 0.0
                println("Prograde (co-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_co-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_co-1)
            else
                println("Retrograde (counter-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". Difference to Comparable KerrN (R/RKerrN-1) = ", R/RKerr_counter-1, ". ω*M = ", w*M, ". ω/ωKerrN - 1 = ", w/wKerr_counter-1)
            end
        else
            if d2h > 0.0 && w > 0.0
                println("Prograde (co-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
            elseif d2h > 0.0 && w < 0.0
                println("Retrograde (counter-rotating) stable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
            elseif d2h < 0.0 && w > 0.0
                println("Prograde (co-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
            else
                println("Retrograde (counter-rotating) unstable Circular Photon Orbit Located at x = ", sol, ". r/rh = ", 2/(1-sol), ". Circumferencial Radius/M = ", R/M, ". ω*M = ", w*M)
            end
        end
    end

    println()

    append!(to_return,to_returnx)
    return to_return

    #return [R_plus/RKerr_plus-1, w_plus/wKerr_plus-1, R_minus/RKerr_minus-1, w_minus/wKerr_minus-1, sol_plus, sol_minus]
end