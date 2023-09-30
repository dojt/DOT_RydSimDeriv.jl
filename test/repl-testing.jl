# Activate .

using Revise

# Activate DOT_RydSimDeriv

using DOT_RydSim
using DOT_RydSim: μs_t, Rad_per_μs_t, Radperμs_per_μs_t
using DOT_RydSimDeriv

# Now activate ./test/

using DOT_NiceMath
using DOT_StatsHelp
using DOT_StatsHelp.Numbers64

using LinearAlgebra: eigvals, eigvecs, Hermitian, Diagonal, normalize
using Random: shuffle!

using Unitful
using Unitful: 	Time, Frequency, μs

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
R_STDDEV = 64
ε        = 1e-4

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
# Pulse shape plots                                                                                     -
#---------------------------------------------------------------------------------------
plotΩ =
let plt = plot()
    (;𝛺ₘₐₓ, 𝛺ᵣₑₛ, 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤, φᵣₑₛ, 𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ) = hw
    (; 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺) = get_hw_data(hw)


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
# Um-Parameterization
#---------------------------------------------------------------------------------------

function um_parameterize(;
                        𝛿𝑥   ::Rad_per_μs_t{ℚ},
                        𝑠ꜛ   ::Radperμs_per_μs_t{ℚ},
                        𝑠ꜜ   ::Radperμs_per_μs_t{ℚ},
                        𝛥𝑡   ::μs_t{ℚ}              ) ::Function

    𝑟 = 1/𝑠ꜛ + 1/𝑠ꜜ
    𝑥(𝑦 ::ℚ) ::Rad_per_μs_t{ℝ} = sign(𝑦)⋅(   𝛥𝑡  -  √( 𝛥𝑡^2 - 2𝑟⋅abs(𝑦) )   ) / 𝑟

    return 𝑥
end

𝑥 =
let
    global hw
	(; 𝑡ᵣₑₛ,𝛥𝑡ₘᵢₙ,𝑡ᵒᶠᶠₘₐₓ,𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺,  𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)

    𝑡ᵒⁿ  = (0//1)μs
    𝑡ᵒᶠᶠ = 𝑡ᵒⁿ + (25//8)𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺

    𝑥 = um_parameterize(; 𝛿𝑥=𝛺ᵣₑₛ, 𝑠ꜛ=hw.𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝑠ꜜ=hw.𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
                        𝛥𝑡=𝑡ᵒᶠᶠ + 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺)
end

#---------------------------------------------------------------------------------------
# Simple Plots
#---------------------------------------------------------------------------------------
plotΩ =
let
    global hw
	(; 𝑡ᵣₑₛ,𝛥𝑡ₘᵢₙ,𝑡ᵒᶠᶠₘₐₓ,𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺,  𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)

    𝑡ᵒⁿ  = (0//1)μs
    𝑡ᵒᶠᶠ = 𝑡ᵒⁿ + (25//8)𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺
    evΩ = Evolution_Ω( 𝑡ᵒⁿ , 𝑡ᵒᶠᶠ
                      ; 𝛥 = 𝛥ₘₐₓ, ε, hw ,
                       𝑇 = 𝑡ᵒᶠᶠ + 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ𝛺)

    scatter( 𝛺 -> let ψ = copy(ψ₀)
                evf(𝛺, evΩ ; 𝚷,R,ψ)
             end,
             -𝛺ₘₐₓ: 7𝛺ᵣₑₛ :+𝛺ₘₐₓ
	     ;
             label="",
	     markersize=0.5, markerstrokewidth=0,
	     xaxis="𝛺")
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
