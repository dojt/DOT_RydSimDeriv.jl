# Activate .

using Revise

# Activate DOT_RydSimDeriv

using DOT_RydSim
using DOT_RydSim: μs_t, Rad_per_μs_t, Radperμs_per_μs_t
using DOT_RydSimDeriv

# Now activate ./test/

using Test

using DOT_NiceMath
using DOT_StatsHelp
using DOT_StatsHelp.Numbers64

using LinearAlgebra: eigvals, eigvecs, Hermitian, Diagonal, normalize
using Random: shuffle!

using Unitful
using Unitful: 	Time, Frequency, μs

using Zygote, UnitfulChainRules


using GenericFFT

function make_Pauli!( H ::Hermitian{ℂ,Matrix{ℂ}} ) :: Nothing
    let d  = size(H) |> first,
	U  = eigvecs(H),
        𝜆⃗  = [  ℝ( (-1)^j  )    for j=1:d  ]
        shuffle!( 𝜆⃗ )
        H .= Hermitian( U⋅Diagonal(𝜆⃗)⋅U' )
    end
    nothing;
end

hw = load_hw(;
             Ω_downslew_factor = 1//1  #=1//3=#,
	     Δ_downslew_factor = 1//1  #=1//2=#)

N_ATOMS  = 2
R_STDDEV = 0 # 64
ε        = 1e-3

ψ₀ = randn(ℂ,2^N_ATOMS) |> normalize
R  = let A = randn(ℂ,2^N_ATOMS,2^N_ATOMS) ; Hermitian( (A+A')⋅R_STDDEV/2 ) end
𝚷  = let A = randn(ℂ,2^N_ATOMS,2^N_ATOMS) ; Hermitian( (A+A')         /2 ) end
make_Pauli!(𝚷)

println("𝜆⃗(R) = ", eigvals(R) )
println("𝜆⃗(𝚷) = ", eigvals(𝚷) )

#---------------------------------------------------------------------------------------
#                                                                                      -
#---------------------------------------------------------------------------------------

using Plots
plotly();

#---------------------------------------------------------------------------------------
# Pulse shape plots                                                                    -
#---------------------------------------------------------------------------------------
plotΩ =
let plt = plot()
    (;𝛺ₘₐₓ, 𝛺ᵣₑₛ, 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤, φᵣₑₛ, 𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ) = hw

    𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺 = (1/𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤 + 1/𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤)⋅𝛺ₘₐₓ
    𝑡ᵒⁿ  = (0//1)μs
    𝑡ᵒᶠᶠ = 𝑡ᵒⁿ + (25//8)𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺
    𝑇    = 𝑡ᵒᶠᶠ + 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺

    p = Pulse__Ω_BangBang{ℚ,ℝ}( 𝑡ᵒⁿ, 𝑡ᵒᶠᶠ, 𝑇, 9//10⋅𝛺ₘₐₓ
				;   𝛺ₘₐₓ, 𝛺ᵣₑₛ,
				𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
				φᵣₑₛ,
				𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ)
    DOT_RydSim._check(p)

    plot!(plt,
	  𝑡 -> p(𝑡) , 0.0μs: 0.0001μs :𝑇,
	  ; label="𝛺",
	  color=:blue)

    (;x⃗,y⃗) = plotpulse(p)
    scatter!(plt,
	     x⃗, y⃗
	     ; label="",
	     color=:blue,
	     markerstrokewidth=0)
    println("Phase = $(phase(p))")

    plt
end

#---------------------------------------------------------------------------------------
# Um-Parameterization                                                                  -
#---------------------------------------------------------------------------------------

# from Pluto notebook "ramped-functions.jl" (in ShiftRules)
function um_parameterize(;
                        𝛿𝑥   ::Rad_per_μs_t{ℚ},
                        𝑠ꜛ   ::Radperμs_per_μs_t{ℚ},
                        𝑠ꜜ   ::Radperμs_per_μs_t{ℚ},
                        𝛥𝑡   ::μs_t{ℚ}                                # point where pulse touches 0 -- not compatible
                                                    ) ::NamedTuple

    𝑟 = 1/𝑠ꜛ + 1/𝑠ꜜ

    (
        u(𝑥 ::Rad_per_μs_t{𝕂}) ::𝕂
    ) where{𝕂}                                    = ustrip(NoUnits,   𝛥𝑡⋅𝑥 − sign(𝑥)/2 ⋅ 𝑟⋅𝑥^2   )

    (
        𝑥(u ::𝕂) ::Rad_per_μs_t{ℝ}
    ) where{𝕂}                                    = sign(u)⋅(   𝛥𝑡  −  √( 𝛥𝑡^2 - 2𝑟⋅abs(u) )   ) / 𝑟

    (
        ∂u(𝑥 ::Rad_per_μs_t{𝕂}) ::μs_t{𝕂}
    ) where{𝕂}                                    = 𝛥𝑡 − 𝑟⋅abs(𝑥)

    (
    ∂𝑥(u ::𝕂) ::Rad_per_μs_t
    ) where{𝕂}                                    = 1/∂u(𝑥(u))

    round(𝑥) = δround(𝑥;𝛿=𝛿𝑥)
    return (   𝑥=round∘𝑥,   u,  ∂u, ∂𝑥, ␣𝑥=𝑥   )
end

# test it
@testset "Testing um-parameterization" begin
    let randℚ()::ℚ = (  num=rand(1:1_000); den=rand(num:1_000); num//den  )
        for _o_iter = 1:1
            𝛥𝑡      = randℚ()μs
            𝑠ꜛ      = randℚ()/1000μs^2
            𝑠ꜜ      = randℚ()/1000μs^2
            𝑟       = 1/𝑠ꜛ + 1/𝑠ꜜ
            twor    = 2𝑟

            for _i_iter = 1:2
                let 𝑥 = rand()*2𝛥𝑡/𝑟 - 𝛥𝑡/𝑟
	            @assert -𝛥𝑡/𝑟 ≤ 𝑥 ≤ 𝛥𝑡/𝑟
	            u = um.u(𝑥)
	            @assert -𝛥𝑡^2/2𝑟 ≤ u ≤ 𝛥𝑡^2/2𝑟
                    @test  1.0/μs+ um.𝑥(u)  ≈ 1.0/μs+ 𝑥             atol=1e-3/μs
		    @test  1.0μs+  um.∂u(𝑥) ≈ 1.0μs+  (um.u)'(𝑥)
                end
                let
		    u = ustrip(NoUnits,  rand()*2𝛥𝑡^2/2𝑟 − 𝛥𝑡^2/2𝑟  )
		    @assert -𝛥𝑡^2/2𝑟 ≤ u ≤ 𝛥𝑡^2/2𝑟
                    𝑥  = um.𝑥(u)
                    @test  1.0+    um.u(𝑥)  ≈ 1.0+    u             atol=1e-3
		    @test  1.0/μs+ um.∂𝑥(u) ≈ 1.0/μs+ (um.␣𝑥)'(u)
                end
            end
        end
    end
end

let
    global hw
	(; 𝑡ᵣₑₛ,𝛥𝑡ₘᵢₙ,𝑡ᵒᶠᶠₘₐₓ,𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺,  𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)

    𝑡ᵒⁿ  = (0//1)μs
    𝑡ᵒᶠᶠ = 𝑡ᵒⁿ + (25//8)𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺

    um_parameterize(; 𝛿𝑥=𝛺ᵣₑₛ, 𝑠ꜛ=hw.𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝑠ꜜ=hw.𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
                    𝛥𝑡=𝑡ᵒᶠᶠ + 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺)
end


#---------------------------------------------------------------------------------------
# Some tests                                                                           -
#---------------------------------------------------------------------------------------
let
    global hw
    (  ; 𝑡ᵣₑₛ,𝛥𝑡ₘᵢₙ,𝑡ᵒᶠᶠₘₐₓ,
       𝛿𝑡ᵉᶠᶠₘᵢₙ𝛺,𝛿𝑡ᵉᶠᶠₘₐₓ𝛺, 𝛿𝑡ᵉᶠᶠₘᵢₙ𝛥,𝛿𝑡ᵉᶠᶠₘₐₓ𝛥,
       𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ                  ) = get_hw_data(hw)

    @test is_δrounded(𝛿𝑡ᵉᶠᶠₘᵢₙ𝛺 ; 𝛿=𝑡ᵣₑₛ)
    @test is_δrounded(𝛿𝑡ᵉᶠᶠₘₐₓ𝛺 ; 𝛿=𝑡ᵣₑₛ)
    @test is_δrounded(𝛿𝑡ᵉᶠᶠₘᵢₙ𝛥 ; 𝛿=𝑡ᵣₑₛ)
    @test is_δrounded(𝛿𝑡ᵉᶠᶠₘₐₓ𝛥 ; 𝛿=𝑡ᵣₑₛ)

    @test 𝛿𝑡ᵉᶠᶠₘᵢₙ𝛺 < 𝛿𝑡ᵉᶠᶠₘₐₓ𝛺
    @test 𝛿𝑡ᵉᶠᶠₘᵢₙ𝛥 < 𝛿𝑡ᵉᶠᶠₘₐₓ𝛥
end

let
    function rand_range( rg ::StepRange)
        case = rand(1:3)
        if      case == 1
            return rg.start
        elseif  case == 2
            rand(rg)
        elseif  case == 3
            return rg.stop
        end
    end


    global hw
    (  ; 𝑡ᵣₑₛ,𝛥𝑡ₘᵢₙ,𝑡ᵒᶠᶠₘₐₓ,
       𝛿𝑡ᵉᶠᶠₘᵢₙ𝛺,𝛿𝑡ᵉᶠᶠₘₐₓ𝛺, 𝛿𝑡ᵉᶠᶠₘᵢₙ𝛥,𝛿𝑡ᵉᶠᶠₘₐₓ𝛥,
       𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ                    ) = get_hw_data(hw)
    (  ; 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
       𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤            ) = hw

    for _1iter = 1:100
        𝑡ᵒⁿ   = rand_range( (0//1)μs : 𝑡ᵣₑₛ : 3𝛥𝑡ₘᵢₙ )     ;  𝑡ᵒⁿ  > 𝛥𝑡ₘᵢₙ  ||  ( 𝑡ᵒⁿ  = 𝛥𝑡ₘᵢₙ )
        𝑡ᵒᶠᶠ  = rand_range( 𝑡ᵒⁿ+𝛥𝑡ₘᵢₙ : 𝑡ᵣₑₛ : 𝑡ᵒᶠᶠₘₐₓ )   ;  𝑡ᵒᶠᶠ > 𝛥𝑡ₘᵢₙ  ||  ( 𝑡ᵒᶠᶠ = 𝛥𝑡ₘᵢₙ )
        𝛿𝑡ᵉᶠᶠ = rand_range( 𝛿𝑡ᵉᶠᶠₘᵢₙ𝛥 : 𝑡ᵣₑₛ : 𝛿𝑡ᵉᶠᶠₘₐₓ𝛥/10 )
        𝛺     = rand_range( (0//1)/μs : 𝛺ᵣₑₛ : 𝛺ₘₐₓ )
        𝑇     = max( 𝑡ᵒⁿ + 𝛿𝑡ᵉᶠᶠ + 𝛥ₘₐₓ⋅(1/𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤 + 1/𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤)/2 ,
                     𝑡ᵒᶠᶠ + 𝛺/𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤                                 )

        evΔ = Evolution_Δ( 𝑡ᵒⁿ , 𝛿𝑡ᵉᶠᶠ
                           ; 𝛺, Ω_𝑡ᵒᶠᶠ=𝑡ᵒᶠᶠ,
                           ε, hw, 𝑇)
        @test true

        for _2iter = 1:100
            𝛥 = rand_range( (0//1)/μs : 𝛥ᵣₑₛ : 𝛥ₘₐₓ )
            ψ = copy(ψ₀)
            evf(𝛥, evΔ ; 𝚷,R,ψ)
            @test true
        end
    end
end


#---------------------------------------------------------------------------------------
# Simple Plots                                                                         -
#---------------------------------------------------------------------------------------
plotΩ =
let
    global hw
    (; 𝑡ᵣₑₛ,𝛥𝑡ₘᵢₙ,𝑡ᵒᶠᶠₘₐₓ, 𝛿𝑡ᵉᶠᶠₘᵢₙ𝛺,𝛿𝑡ᵉᶠᶠₘₐₓ𝛺,  𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)

    # 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺 = (1/𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤 + 1/𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤)⋅𝛺ₘₐₓ
    𝑡ᵒⁿ  = (0//1)μs
    𝛿𝑡ᵉᶠᶠ = ..................................................................................
    𝑇    = 𝑡ᵒᶠᶠ + 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺
    evΩ = Evolution_Ω( 𝑡ᵒⁿ , 𝑡ᵒᶠᶠ
                      ; 𝛥 = (0//1)/μs#=𝛥ₘₐₓ=#, ε, hw ,
                       𝑇)

    𝑓(𝛺) = let ψ = copy(ψ₀)
               evf(𝛺, evΩ ; 𝚷,R,ψ)
           end

    plt1 = scatter( 𝛺 -> 𝑓(𝛺),
                    -𝛺ₘₐₓ: 7𝛺ᵣₑₛ :+𝛺ₘₐₓ
	            ; label="",
                    color=:blue,
	            markersize=0.5, markerstrokewidth=0,
	            xaxis="𝛺")

    #
    #"""This doesn't work: It needs other `evf()`"""
    #
    # um = um_parameterize(; 𝛿𝑥=𝛺ᵣₑₛ, 𝑠ꜛ=hw.𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝑠ꜜ=hw.𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
    #                      𝛥𝑡=𝑇)
    # plt2 = scatter(u -> 𝑓(um.𝑥(u)),
    #                 um.u(-𝛺ₘₐₓ): 0.01 :um.u(+𝛺ₘₐₓ)
    #                 ; label="",
    #                 color=:red,
    #                 markersize=0.5, markerstrokewidth=0,
    #                 xaxis="u(𝛺)")
    plt1
end

plotΔ =
let
    global hw
    (; 𝑡ᵣₑₛ,𝛥𝑡ₘᵢₙ,𝑡ᵒᶠᶠₘₐₓ,  𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)

    𝛺   = -𝛺ₘₐₓ/100
    𝑡ᵒⁿ  = max(𝛥𝑡ₘᵢₙ, 10⋅𝑡ᵣₑₛ)
    𝑡ᵒᶠᶠ = min( 10000⋅𝑡ᵣₑₛ,
                𝑡ᵒᶠᶠₘₐₓ - get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺(hw;𝛺,𝛥=𝛥ₘₐₓ))

    @show 𝑡ᵒⁿ,𝑡ᵒᶠᶠ,𝑡ᵒᶠᶠₘₐₓ


    evΔ = Evolution_Δ( 𝑡ᵒⁿ , 𝑡ᵒᶠᶠ
                      ; 𝛺, ε, hw )

    scatter( 𝛥 -> let ψ = copy(ψ₀)
                evf(𝛥, evΔ ; 𝚷,R,ψ)
             end,
             -𝛥ₘₐₓ: 100001𝛥ᵣₑₛ :+𝛥ₘₐₓ
	     ;
             label="",
	     markersize=0.5, markerstrokewidth=0,
	     xaxis="𝛥")
end

#---------------------------------------------------------------------------------------
# Fourier Transform
#---------------------------------------------------------------------------------------

plotΩ_f̂ =
let
    global hw
    (; 𝑡ᵣₑₛ,𝛥𝑡ₘᵢₙ,𝑡ᵒᶠᶠₘₐₓ, 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺,  𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)

    𝑡ᵒⁿ  = (0//1)μs
    𝑡ᵒᶠᶠ = 𝑡ᵒⁿ + (25//8)𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺
    evΩ  = Evolution_Ω( 𝑡ᵒⁿ , 𝑡ᵒᶠᶠ
                        ;
                        𝛥 = 𝛥ₘₐₓ, ε=0.001, hw,
                        𝑇 = 𝑡ᵒᶠᶠ + 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺 )

    Set_of_𝛺s = (  -𝛺ₘₐₓ: 100𝛺ᵣₑₛ :+𝛺ₘₐₓ  )
    N = length(Set_of_𝛺s)
    f⃗ = [ let ψ = copy(ψ₀)
             evf(𝛺, evΩ ; 𝚷,R,ψ)
          end
          for 𝛺 ∈ Set_of_𝛺s       ]

    f̂ = fft(f⃗) ./ N

    Ωₘₐₓ = ustrip(u"μs^(-1)", 𝛺ₘₐₓ)
    scatter(  [ ( k==0 ? 0.0 : 2π/(2Ωₘₐₓ/k) )   for k = 0:length(f̂)÷2 ],
              abs.(@view f̂[1:1+length(f̂)÷2])
              ; label="",
              color=:blue,
	      markersize=0.5, markerstrokewidth=0 )

    λ = DOT_RydSimDeriv.λ(evΩ)
    K = 2π/λ
    scatter!([K],[0.0]
             ; label="",
             color=:red,
	     markersize=1, markerstrokewidth=0 )
end

plotΔ_f̂ =
let
    global hw
    (; 𝑡ᵣₑₛ,𝛥𝑡ₘᵢₙ,𝑡ᵒᶠᶠₘₐₓ,  𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)

    𝛺   = 𝛺ₘₐₓ/10
    𝑡ᵒⁿ  = max(𝛥𝑡ₘᵢₙ, 10⋅𝑡ᵣₑₛ)
    𝑡ᵒᶠᶠ = min(       10000⋅𝑡ᵣₑₛ, 𝑡ᵒᶠᶠₘₐₓ- get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺(hw;𝛺,𝛥=𝛥ₘₐₓ))
    evΔ  = Evolution_Δ( 𝑡ᵒⁿ , 𝑡ᵒᶠᶠ
                        ; 𝛺, ε=1e-4, hw )

    Set_of_𝛥s = (   𝛥ᵣₑₛ⋅ 𝛥  for 𝛥 ∈ range(start=-𝛥ₘₐₓ/𝛥ᵣₑₛ, stop=+𝛥ₘₐₓ/𝛥ᵣₑₛ, length=5001)   )
    N = length(Set_of_𝛥s)
    f⃗ = [ let ψ = copy(ψ₀)
             evf(𝛥, evΔ ; 𝚷,R,ψ)
          end
          for 𝛥 ∈ Set_of_𝛥s       ]

    f̂ = fft(f⃗) ./ N

    𝛥ₘₐₓ = ustrip(u"μs^(-1)", 𝛥ₘₐₓ)
    scatter(  [ ( k==0 ? 0.0 : 2π/(2𝛥ₘₐₓ/k) )   for k = 0:length(f̂)÷2 ],
              abs.(@view f̂[1:1+length(f̂)÷2])
              ; label="",
              color=:blue,
	      markersize=0.5, markerstrokewidth=0 )

    λ = DOT_RydSimDeriv.λ(evΔ)
    K = 2π/λ
    scatter!([K],[0.0]
             ; label="",
             color=:red,
	     markersize=1, markerstrokewidth=0 )
end

#EOF
