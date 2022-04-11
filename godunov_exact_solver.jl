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
    hL = 0.5
    uL = 0.1
    αL = 1

    hR = 0.0
    uR = 0.00
    αR = 0


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

    function compute_pressure(hL, hR, uL, uR, αL, αR)
        ΠL = g * hL^2 / 2 + αL * S
        ΠR = g * hR^2 / 2 + αR * S
        # on traite le cas où il y a zone sèche car résultat immédiat on connait directement h* donc Π* pas besoin d'itérer
        if abs(ΠL) < myzero && abs(ΠR) > myzero
            ΠStar = zero(hL) 
        elseif abs(ΠR) < myzero && abs(ΠL) > myzero
            ΠStar = zero(hL) 
        elseif abs(ΠR) > myzero && abs(ΠL) > myzero
            ΠStar = Newton_Raphson(hL, hR, uL, uR, αL, αR, ΠL, ΠR, tol)
        else
            error("No liquid")
        end
        return ΠL, ΠR, ΠStar
    end

    function Newton_Raphson(hL, hR, uL, uR, αL, αR, ΠL, ΠR, tol)
        ig = 0.5 * (ΠL + ΠR)
        count_newton = 0
        left_wave, right_wave = compute_wave_state(ΠL, ΠR, ig)
        f, df, fL, fR = define_pressure_function(hL, hR, uL, uR, αL, αR, left_wave, right_wave)
        ΠStar = ig - f(ig)/(df(ig) + myzero)
        while abs(ig - ΠStar)/(0.5*(ig + ΠStar)) > tol 
            ig = ΠStar
            left_wave, right_wave = compute_wave_state(ΠL, ΠR, ig)
            f, df, fL, fR = define_pressure_function(hL, hR, uL, uR, αL, αR, left_wave, right_wave)
            ΠStar = ig - f(ig)/(df(ig) + myzero)
            count_newton += 1
            if count_newton > limit_Newton
                error("Newton do not converge")
            end
        end
        @show count_newton
        display(plot(f, -0.1 , 0.1))
        #error("")
        return ΠStar
    end




    function compute_wave_state(ΠL, ΠR, ΠStar)
        #if dry region ΠStar = myzero ce qui permet de simplifier les choses
        if abs(ΠR) < myzero
            right_wave = dry()
            # ΠStar = myzero
            if ΠL ≥ myzero
                left_wave = rarefaction()
            else
                left_wave = shock()
            end
        elseif abs(ΠL) < myzero
            left_wave = dry()
            # ΠStar = myzero
            if ΠR ≥ myzero
                right_wave = rarefaction()
            else
                right_wave = shock()
            end
        else 
            if ΠL ≥ ΠStar
                left_wave = rarefaction()
            else
                left_wave = shock()
            end
            if ΠR ≥ ΠStar
                right_wave = rarefaction()
            else
                right_wave = shock()
            end
        end
        return left_wave, right_wave
    end

    function define_pressure_function(hL, hR, uL, uR, αL, αR, ::shock, ::rarefaction)
        fL(p) = (√(2*(p-αL * S)/ (αL*g)) - hL) / √(hL * √(2*(p-αL * S)/ (αL*g))) * √(g/2 * (hL + √(2*(p-αL * S)/ (αL*g))))
        fR(p) = -2* (√(g * hR) - (2 * g * (p - αR * S) / αR))^1/4
        hLStar(p) = √(2*(p - αL * S) / (αL * g))
        dhLStar(p) = 1/ √(2 * αL * g * (p - αL * S))
        dfL(p) = dhLStar(p) * ((2 * hL * hLStar(p) - (hLStar(p) - hL)) / ( 2 * (hL * hLStar(p))^3/2) * √(g / 2 * (hLStar(p) + hL)) + g / 4 * (hLStar(p) - hL) / √(g / 2 * ( hL + hLStar(p)) * hLStar(p) * hL))
        dfR(p) = 0.5 * g / αR * (p - αR * S)^(-3/4)
        f(p) = fL(p) + fR(p) + uR - uL
        df(p) = dfL(p) + dfR(p)
        return f, df, fL, fR
    end

    function define_pressure_function(hL, hR, uL, uR, αL, αR, ::shock, ::shock)
        fL(p) = (√(2*(p-αL * S)/ (αL*g)) - hL) / √(hL * √(2*(p-αL * S)/ (αL*g))) * √(g/2 * (hL + √(2*(p-αL * S)/ (αL*g))))
        fR(p) = (√(2*(p-αR * S)/ (αR*g)) - hR) / √(hR * √(2*(p-αR * S)/ (αR*g))) * √(g/2 * (hR + √(2*(p-αR * S)/ (αR*g))))
        hLStar(p) = √(2*(p - αL * S) / (αL * g))
        dhLStar(p) = 1/ √(2 * αL * g * (p - αL * S))
        hRStar(p) = √(2*(p - αR * S) / (αR * g))
        dhRStar(p) = 1/ √(2 * αR * g * (p - αR * S))
        dfL(p) = dhLStar(p) * ((2 * hL * hLStar(p) - (hLStar(p) - hL)) / ( 2 * (hL * hLStar(p))^3/2) * √(g / 2 * (hLStar(p) + hL)) + g / 4 * (hLStar(p) - hL) / √(g / 2 * ( hL + hLStar(p)) * hLStar(p) * hL))
        dfR(p) = dhRStar(p) * ((2 * hR * hRStar(p) - (hRStar(p) - hR)) / ( 2 * (hR * hRStar(p))^3/2) * √(g / 2 * (hRStar(p) + hR)) + g / 4 * (hRStar(p) - hR) / √(g / 2 * ( hR + hRStar(p)) * hRStar(p) * hR))
        f(p) = fL(p) + fR(p) + uR - uL
        df(p) = dfL(p) + dfR(p)
        return f, df, fL, fR
    end

    function define_pressure_function(hL, hR, uL, uR, αL, αR, ::rarefaction, ::rarefaction)
        fL(p) = -2* (√(g * hL) - (2 * g * (p - αL * S) / αL))^1/4
        fR(p) = -2* (√(g * hR) - (2 * g * (p - αR * S) / αR))^1/4
        dfL(p) = 0.5 * g / αL * (p - αL * S)^(-3/4)
        dfR(p) = 0.5 * g / αR * (p - αR * S)^(-3/4)
        f(p) = fL(p) + fR(p) + uR - uL
        df(p) = dfL(p) + dfR(p)
        return f, df, fL, fR
    end

    function define_pressure_function(hL, hR, uL, uR, αL, αR, ::rarefaction, ::shock)
        fL(p) = -2* (√(g * hL) - (2 * g * (p - αL * S) / αL))^1/4
        fR(p) = (√(2*(p-αR * S)/ (αR*g)) - hR) / √(hR * √(2*(p-αR * S)/ (αR*g))) * √(g/2 * (hR + √(2*(p-αR * S)/ (αR*g))))
        hRStar(p) = √(2*(p - αR * S) / (αR * g))
        dhRStar(p) = 1/ √(2 * αR * g * (p - αR * S))
        dfL(p) = 0.5 * g / αL * (p - αL * S)^(-3/4)
        dfR(p) = dhRStar(p) * ((2 * hR * hRStar(p) - (hRStar(p) - hR)) / ( 2 * (hR * hRStar(p))^3/2) * √(g / 2 * (hRStar(p) + hR)) + g / 4 * (hRStar(p) - hR) / √(g / 2 * ( hR + hRStar(p)) * hRStar(p) * hR))
        f(p) = fL(p) + fR(p) + uR - uL
        df(p) = dfL(p) + dfR(p)
        return f, df, fL, fR
    end



    function compute_state!(h, u, α, hL, hR, uL, uR, αL, αR, ΠStar, t, ::rarefaction, ::dry)
        uStar = uL + 2*(√(g*hL) - (-2*αL * g * S)^1/4)
        hStar = √(-2*αL*S/g)
        STL   =  uStar - √(g * hStar)
        SHL   = uL - √(g * hL)  
        for icell in 1:ncell
            if xcell[icell] / t ≤ SHL  
                h[icell], u[icell], α[icell] = left_state()
            elseif SHL ≤ xcell[icell] / t ≤ STL
                h[icell], u[icell], α[icell] = fan_left_state(t, xcell[icell])
            elseif STL ≤ xcell[icell] / t ≤ uStar
                h[icell], u[icell], α[icell] = star_right_dry_state(left_wave)
            elseif xcell[icell] / t ≥ uStar
                h[icell], u[icell], α[icell] = right_state()
            end
        end
    end

    function compute_state!(h, u, α, hL, hR, uL, uR, αL, αR, ΠStar, t, ::shock, ::dry)
        hStar = √(-2*αL*S/g)
        uStar = uL - (hStar - hL) / √(hStar * hL + myzero) * √(g/2*(hStar + hL))
        S1   =  uL - √(g*hStar*(hStar+hL)/(2*hL))
        for icell in 1:ncell
            if xcell[icell] / t ≤ S1  
                h[icell], u[icell], α[icell] = left_state()
            elseif S1 ≤ xcell[icell] / t ≤ uStar
                h[icell], u[icell], α[icell] = star_right_dry_state(left_wave)
            elseif xcell[icell] / t ≥ uStar
                h[icell], u[icell], α[icell] = right_state()
            end
        end
    end


    function compute_state!(h, u, α, hL, hR, uL, uR, αL, αR, ΠStar, t, ::dry, ::shock)
        hStar = √(-2*αR*S/g)
        uStar = uR + (hStar - hR) / √(hStar * hR + myzero) * √(g/2*(hStar + hR))
        S3   =  uR + √(g*hStar*(hStar+hR)/(2*hR))
        for icell in 1:ncell
            if xcell[icell] / t ≤ uStar  
                h[icell], u[icell], α[icell] = left_state()
            elseif uStar ≤ xcell[icell] / t ≤ S3
                h[icell], u[icell], α[icell] = star_left_dry_state(right_wave)
            elseif xcell[icell] / t ≥ S3
                h[icell], u[icell], α[icell] = right_state()
            end
        end
    end

    function compute_state!(h, u, α, hL, hR, uL, uR, αL, αR, ΠStar, t, ::dry, ::rarefaction)
        uStar =  uR - 2*(√(g*hR) - (-2*αR*S*g)^1/4)
        hStar = √(-2*αR*S/g)
        STR   =  uStar + √(g * hStar)
        SHR   = uR + √(g * hR)  
        for icell in 1:ncell
            if xcell[icell] / t ≤ uStar  
                h[icell], u[icell], α[icell] = left_state()
            elseif uStar ≤ xcell[icell] / t ≤ STR
                h[icell], u[icell], α[icell] = star_left_dry_state(right_wave)
            elseif STR ≤ xcell[icell] / t ≤ SHR
                h[icell], u[icell], α[icell] = fan_right_state(t, xcell[icell])
            elseif xcell[icell] / t ≥ SHR
                h[icell], u[icell], α[icell] = right_state()
            end
        end
    end

    function compute_state!(h, u, α, hL, hR, uL, uR, αL, αR, ΠStar, t, ::rarefaction, ::shock)
        fL = define_pressure_function(hL, hR, uL, uR, αL, αR, rarefaction(), shock())[3]
        fR = define_pressure_function(hL, hR, uL, uR, αL, αR, rarefaction(), shock())[4]
        uStar = 0.5 * (uR + uL) + 0.5* (fR(ΠStar) + fL(ΠStar)) 
        hLStar = √(2*(ΠStar - αL * S) / (αL * g))
        hRStar = √(2*(ΠStar - αR * S) / (αR * g))
        STL   =  uStar + √(g * hLStar)
        SHL   = uL + √(g * hL)  
        S3   =  uR + √(g*hRStar*(hRStar+hR)/(2*hR))
        for icell in 1:ncell
            if xcell[icell] / t ≤ SHL  
                h[icell], u[icell], α[icell] = left_state()
            elseif SHL ≤ xcell[icell] / t ≤ STL
                h[icell], u[icell], α[icell] = fan_left_state(t, xcell[icell])
            elseif STL ≤ xcell[icell] / t ≤ uStar
                h[icell], u[icell], α[icell] = star_left_state(uStar, hLStar, αL) 
            elseif uStar ≤ xcell[icell] / t ≤ S3
                h[icell], u[icell], α[icell] = star_right_state(uStar, hRStar, αR)
            elseif xcell[icell] / t ≥ S3
                h[icell], u[icell], α[icell] = right_state()
            end
        end
        return h, u, α
    end
    
    function compute_state!(h, u, α, hL, hR, uL, uR, αL, αR, ΠStar, t, ::shock, ::shock)
        fL = define_pressure_function(hL, hR, uL, uR, αL, αR, rarefaction(), shock())[3]
        fR = define_pressure_function(hL, hR, uL, uR, αL, αR, rarefaction(), shock())[4]
        uStar = 0.5 * (uR + uL) + 0.5* (fR(ΠStar) + fL(ΠStar)) 
        hLStar = √(2*(ΠStar - αL * S) / (αL * g))
        hRStar = √(2*(ΠStar - αR * S) / (αR * g))
        S1   =  uL - √(g*hLStar*(hLStar+hL)/(2*hL))
        S3   =  uR + √(g*hRStar*(hRStar+hR)/(2*hR))
        for icell in 1:ncell
            if xcell[icell] / t ≤ S1  
                h[icell], u[icell], α[icell] = left_state()
            elseif S1 ≤ xcell[icell] / t ≤ uStar
                h[icell], u[icell], α[icell] = star_left_state(uStar, hLStar, αL) 
            elseif uStar ≤ xcell[icell] / t ≤ S3
                h[icell], u[icell], α[icell] = star_right_state(uStar, hRStar, αR)
            elseif xcell[icell] / t ≥ S3
                h[icell], u[icell], α[icell] = right_state()
            end
        end
        return h, u, α
    end
    
    function compute_state!(h, u, α, hL, hR, uL, uR, αL, αR, ΠStar, t, ::rarefaction, ::rarefaction)
        fL = define_pressure_function(hL, hR, uL, uR, αL, αR, rarefaction(), shock())[3]
        fR = define_pressure_function(hL, hR, uL, uR, αL, αR, rarefaction(), shock())[4]
        uStar = 0.5 * (uR + uL) + 0.5* (fR(ΠStar) + fL(ΠStar)) 
        hLStar = √(2*(ΠStar - αL * S) / (αL * g))
        hRStar = √(2*(ΠStar - αR * S) / (αR * g))
        STL   =  uStar + √(g * hLStar)
        SHL   = uL + √(g * hL)  
        STR   =  uStar + √(g * hRStar)
        SHR   = uR + √(g * hR)  
        for icell in 1:ncell
            if xcell[icell] / t ≤ SHL  
                h[icell], u[icell], α[icell] = left_state()
            elseif SHL ≤ xcell[icell] / t ≤ STL
                h[icell], u[icell], α[icell] = fan_left_state(t, xcell[icell])
            elseif STL ≤ xcell[icell] / t ≤ uStar
                h[icell], u[icell], α[icell] = star_left_state(uStar, hLStar, αL) 
            elseif uStar ≤ xcell[icell] / t ≤ STR
                h[icell], u[icell], α[icell] = star_right_state(uStar, hRStar, αR)
            elseif STR ≤ xcell[icell] / t ≤ SHR
                h[icell], u[icell], α[icell] = fan_right_state(t, xcell[icell])
            elseif xcell[icell] / t ≥ SHR
                h[icell], u[icell], α[icell] = right_state()
            end
        end
        return h, u, α
    end
    
    function compute_state!(h, u, α, hL, hR, uL, uR, αL, αR, ΠStar, t, ::shock, ::rarefaction)
        fL = define_pressure_function(hL, hR, uL, uR, αL, αR, rarefaction(), shock())[3]
        fR = define_pressure_function(hL, hR, uL, uR, αL, αR, rarefaction(), shock())[4]
        uStar = 0.5 * (uR + uL) + 0.5* (fR(ΠStar) + fL(ΠStar)) 
        hLStar = √(2*(ΠStar - αL * S) / (αL * g))
        hRStar = √(2*(ΠStar - αR * S) / (αR * g))
        S1   =  uL - √(g*hLStar*(hLStar+hL)/(2*hL))
        STR   =  uStar + √(g * hRStar)
        SHR   = uR + √(g * hR) 
        for icell in 1:ncell
            if xcell[icell] / t ≤ S1  
                h[icell], u[icell], α[icell] = left_state()
            elseif S1 ≤ xcell[icell] / t ≤ uStar
                h[icell], u[icell], α[icell] = star_left_state(uStar, hLStar, αL) 
            elseif uStar ≤ xcell[icell] / t ≤ STR
                h[icell], u[icell], α[icell] = star_right_state(uStar, hRStar, αR)
            elseif STR ≤ xcell[icell] / t ≤ SHR
                h[icell], u[icell], α[icell] = fan_right_state(t, xcell[icell])
            elseif xcell[icell] / t ≥ SHR
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

    function fan_left_state(t, x)
        h = 1/(9*g) * (uL + 2*√(g * hL) - x/t)^2
        u = 1/3 * (uL + 2 * √(g * hL) + 2*x/t)
        α = αL
        return h, u, α
    end

    function fan_right_state(t, x)
        h = 1/(9*g) * (-uR + 2*√(g * hR) + x/t)^2
        u = 1/3 * (uR - 2 * √(g * hR) + 2*x/t)
        α = αR
        return h, u, α
    end


    function star_right_dry_state(::rarefaction)
        h = √(-2*αL*S/g)
        u = uL + 2*(√(g*hL) - (-2*αL*S*g)^1/4)
        α = αL
        return h, u, α
    end

    function star_right_dry_state(::shock)
        h = √(-2*αL*S/g)
        u = uL - (h - hL) / √(h * hL + myzero) * √(g/2*(h + hL))
        α = αL
        return h, u, α
    end

    function star_left_dry_state(::rarefaction)
        h = √(-2*αR*S/g)
        u = uR - 2*(√(g*hR) - (-2*αR*S*g)^1/4)
        α = αR
        return h, u, α
    end

    function star_left_dry_state(::shock)
        h = √(-2*αR*S/g)
        u = uR + (h - hR) / √(h * hR) * √(g/2*(h + hR))
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

    function run!(h, u, α, hth, uth, αth, ΠStar, t, xcell, left_wave, right_wave)
        for n in 1:nite 
            vtk_output(h, u, α, hth, uth, αth, xcell, n, t)
            t += Δt
            if left_wave isa rarefaction && right_wave isa dry
                analytical_solution_dambreak!(hth, uth, αth, xcell, t, left_wave, right_wave)
            end
            compute_state!(h, u, α, hL, hR, uL, uR, αL, αR, ΠStar, t, left_wave, right_wave)
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

    function analytical_solution_dambreak!(hth, uth, αth, x, t, ::rarefaction, ::dry) #demi solution du probleme car on prend abs(xi - xcenter)
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
    ΠL, ΠR, ΠStar = compute_pressure(hL, hR, uL, uR, αL, αR)
    left_wave, right_wave = compute_wave_state(ΠL, ΠR, ΠStar)
    @show ΠStar, ΠL, ΠR
    @show left_wave right_wave
    run!(h, u, α, hth, uth, αth, ΠStar, t, xcell, left_wave, right_wave)
    display_plot(h, u, α, xcell)

end
