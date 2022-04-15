module Riemman   
    using Plots
    plotly()

    using Plots
    using WriteVTK
    plotly()


    abstract type AbstractWaveType end
    struct shock        <: AbstractWaveType end
    struct rarefaction  <: AbstractWaveType end
    struct dry          <: AbstractWaveType end

 
    tol = 1e-6
    gi = 0.1
    limit_Newton = 10000

    xmin = -1
    xmax = 1
    xd = (xmin + xmax)/2
    g  = 9.81 
    θ  = deg2rad(90.0)
    γ = 0.067
    S = γ * (cos(θ) - 1)
    #S = 0
    ncell = 50
    nite = 15
    Δt = 0.01 
    t = 0
    myzero = 1e-12


    hth = zeros(ncell)
    uth = zeros(ncell)
    αth = zeros(ncell)
    h = zeros(ncell)
    u = zeros(ncell)
    α = zeros(ncell)

    xcell = LinRange(xmin,xmax,ncell)
    VTK_folder = "dry_right_shock/"



    # Parameter
    hL = 3
    uL = 0
    αL = 1

    hR = 1
    uR = 0.00
    αR = 0.5


    function init!(h, u, α, hth, uth, αth)
        for icell in 1:ncell
            if xcell[icell] < xd
                h[icell] = hL
                α[icell] = αL
                hth[icell] = hL
                αth[icell] = αL
            else
                h[icell] = hR
                α[icell] = αR
                hth[icell] = hR
                αth[icell] = αR
            end
        end
    end

    function compute_state!(h, u, α, hL, hR, uL, uR, αL, αR, t)
        ΠL = g * hL^2 / 2 + αL * S
        ΠR = g * hR^2 / 2 + αR * S
        SL   =  uL - √(g*hL)
        SR   =  uR + √(g*hR)
        uStar = (ΠR - ΠL + αL * hL * (SL - uL) * uL - αR * hR * (SR - uR) * uR) / (αL * hL * (SL - uL) - αR * hR * (SR - uR))
        @show uStar SL SR
        hLStar = hL * (SL -uL) / (SL - uStar)
        hRStar = hR * (SR -uR) / (SR - uStar)
        for icell in 1:ncell
            if xcell[icell] / t ≤ SL  
                h[icell], u[icell], α[icell] = left_state()
            elseif SL ≤ xcell[icell] / t ≤ uStar
                h[icell], u[icell], α[icell] = star_left_state(uStar, hLStar, αL) 
            elseif uStar ≤ xcell[icell] / t ≤ SR
                h[icell], u[icell], α[icell] = star_right_state(uStar, hRStar, αR)
            elseif xcell[icell] / t ≥ SR
                h[icell], u[icell], α[icell] = right_state()
            end
        end
        return h, u, α
    end
    

    function left_state()
        h = hL
        u = uL
        α = αL
        return h, u, α
    end

    function right_state()
        h = hR
        u = uR
        α = αR
        return h, u, α
    end

    function star_left_state(uStar, hLStar, αL) 
        h = hLStar
        u = uStar
        α = αL
        return h, u, α
    end

    function star_right_state(uStar, hRStar, αR) 
        h = hRStar
        u = uStar
        α = αR
        return h, u, α
    end

    function run!(h, u, α, hth, uth, αth, t, xcell)
        for n in 1:nite 
            vtk_output(h, u, α, hth, uth, αth, xcell, n, t)
            t += Δt
            analytical_solution_dambreak!(hth, uth, αth, xcell, t)
            compute_state!(h, u, α, hL, hR, uL, uR, αL, αR, t)
            println("iter=",n,"  Δt=",Δt,"  t=",t)
        end

        return h, u, α
    end

    function vtk_output(h, u, α, hth, uth, αth, x, n, t)
        if t == 0
            isdir("post_process/"*VTK_folder) && rm("post_process/"*VTK_folder, force=true, recursive= true)
            !isdir("post_process/"*VTK_folder) && mkdir("post_process/"*VTK_folder)
        end

        # -- Output collection
        dir = string(@__DIR__ , "/")
        file = (dir * "post_process/" * VTK_folder * "riemann_problem")
        pvd = paraview_collection(file, append = n > 1)

        # -- Construct the unstructured grid
        xnodes = LinRange(xmin,xmax,ncell + 1)
        cells = [MeshCell(VTKCellTypes.VTK_LINE, [i, i+1]) for i in 1:ncell]

        # -- Output file
        t = round(t, digits = 5)
        file1 = (dir * "post_process/" * VTK_folder  * "$t.riemann_problem.vtu")
        vtkfile = vtk_grid(file1, xnodes, cells)
        
        # -- Values at the cells
        #@show h
        vtkfile["h", VTKCellData()] = h
        vtkfile["u", VTKCellData()] = u
        vtkfile["alpha", VTKCellData()] = α
        vtkfile["hth", VTKCellData()] = hth
        vtkfile["uth", VTKCellData()] = uth
        vtkfile["alphath", VTKCellData()] = αth

        # -- Save each file and collection
        vtk_save(vtkfile)
        pvd[float(t)] = vtkfile
        vtk_save(pvd)
    end


    function display_plot(h, u, α, x)

        plt1 = plot(x, α, title="α")
        scatter!(plt1, xcell, zeros(ncell), ms = 2.5, label=false) 

        plt2 = plot(x, h, title="Height")
        scatter!(plt2, xcell, zeros(ncell), ms = 2.5, label=false) 

        plt3 = plot(x, u, title="Velocity")
        scatter!(plt3, xcell, zeros(ncell), ms = 2.5, label=false) 

        subplot = plot(plt1, plt2, plt3, layout = (3, 1), label="")
        display(subplot)
    end

    function analytical_solution_dambreak!(hth, uth, αth, x, t) #demi solution du probleme car on prend abs(xi - xcenter)
        xa = xd - t*√(g*hL)
        xb = xd + 2*t * √(g*hL)
        for (i,xᵢ) in enumerate(x)
            if xᵢ ≤ xa
                αth[i] = one(hL)
                hth[i] = hL
                uth[i] = zero(hL)
            elseif xa ≤ xᵢ ≤ xb
                αth[i] = one(hL)
                hth[i] = (4/(9*g))*(√(g*hL)-(xᵢ-xd)/(2*t))^2
                uth[i] = (2/3) * ((xᵢ-xd)/t + √(g*hL))
            else
                αth[i] = zero(hL)
                hth[i] = zero(hL)
                uth[i] = zero(hL)
            end
        end
    end


    init!(h, u, α, hth, uth, αth)
    run!(h, u, α, hth, uth, αth, t, xcell)
    display_plot(h, u, α, xcell)

end
