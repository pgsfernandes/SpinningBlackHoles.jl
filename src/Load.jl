function LoadSystem()
	
	if isdefined(Main,:Nx)
		@eval const Nx = Main.Nx
	else
		print("Enter x resolution: ")
		@eval const Nx = parse(Int, chomp(readline()))
		println()
	end

	if isdefined(Main,:Ny)
		@eval const Ny = Main.Ny
	else
		print("Enter y resolution: ")
		@eval const Ny = parse(Int, chomp(readline()))
		println()
	end

    @eval begin
        const X = xpoints(Nx)
        const Mx = GetChebyshev(Nx,X)
        const Y = ypoints(Ny)
        const My = GetTrig(Ny,Y)
    end
	println()
	println("Spinning Black Holes package loaded!")
	print("Resolution: Nx=",Nx,", Ny=",Ny)
	println()
	return X, Mx, Y, My
end