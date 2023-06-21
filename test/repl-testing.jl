using .DOT_RydSimDeriv

using DOT_NiceMath
using DOT_NiceMath.NumbersF64


using Unitful
using Unitful: 	Time, Frequency, μs

using Plots
plotly();


using LinearAlgebra: eigvals, Hermitian, normalize


hw = load_hw(;  Ω_downslew_factor = 1//3,
				Δ_downslew_factor = 1//2)


N_ATOMS  = 1
R_STDDEV = 0#64
ε        = 1e-1

ϕ = randn(ℂ,2^N_ATOMS) |> normalize
ψ = randn(ℂ,2^N_ATOMS) |> normalize

R = let A = randn(ℂ,2^N_ATOMS,2^N_ATOMS) ; Hermitian( (A+A')⋅R_STDDEV/2 ) end

println("λ⃗ = ", eigvals(R) )

#---------------------------------------------------------------------------------------
# Plots
#---------------------------------------------------------------------------------------
plotΩ =
let
    global hw
	(; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)
	𝛥 = (0//1)/μs # 𝛥ᵣₑₛ #-𝛥ₘₐₓ/2
	scatter( 𝛺 -> evf_Ω(𝛺;𝛥,ϕ,R,ψ,ε,hw) , -𝛺ₘₐₓ: 7𝛺ᵣₑₛ :+𝛺ₘₐₓ
			; label="",
			markersize=0.5, markerstrokewidth=0,
			xaxis="𝛺")
end

plotΔ =
let
	(; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)
	𝛺 = -𝛺ₘₐₓ/100
	scatter( 𝛥 -> evf_Δ(𝛥; 𝛺,ϕ,R,ψ,ε,hw) , -𝛥ₘₐₓ: 100001𝛥ᵣₑₛ :+𝛥ₘₐₓ
			; label="",
			markersize=0.5, markerstrokewidth=0,
			xaxis="𝛥")
end

#---------------------------------------------------------------------------------------
# Fourier Transform
#---------------------------------------------------------------------------------------

plotΩ_f̂ =
let
    global hw
	(; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)
	𝛥 = (0//1)/μs # 𝛥ᵣₑₛ #-𝛥ₘₐₓ/2

    Set_of_𝛺s = (  -𝛺ₘₐₓ: 𝛺ᵣₑₛ :+𝛺ₘₐₓ  )
    N = length(Set_of_𝛺s)
    f⃗ = [ evf_Ω(𝛺;𝛥,ϕ,R,ψ,ε,hw)   for 𝛺 ∈ Set_of_𝛺s ]

    f̂ = fft(f⃗) ./ N

    scatter(abs.(f̂)
        ;
        label="",
		markersize=0.5, markerstrokewidth=0 )
end

