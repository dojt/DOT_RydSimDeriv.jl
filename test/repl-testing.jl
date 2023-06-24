#using DOT_RydSimDeriv

using DOT_NiceMath
using DOT_NiceMath.NumbersF64


using Unitful
using Unitful: 	Time, Frequency, μs

using Plots
plotly();


using LinearAlgebra: eigvals, Hermitian, normalize
using GenericFFT


hw = load_hw(;  Ω_downslew_factor = 1//1  #=1//3=#,
				Δ_downslew_factor = 1//1  #=1//2=#)


N_ATOMS  = 1
R_STDDEV = 0#64
ε        = 1e-1

ϕ  = randn(ℂ,2^N_ATOMS) |> normalize
ψ₀ = randn(ℂ,2^N_ATOMS) |> normalize

R = let A = randn(ℂ,2^N_ATOMS,2^N_ATOMS) ; Hermitian( (A+A')⋅R_STDDEV/2 ) end

println("λ⃗ = ", eigvals(R) )

#---------------------------------------------------------------------------------------
# Plots
#---------------------------------------------------------------------------------------
plotΩ =
let
    global hw
	(; 𝑡ᵣₑₛ,𝛥𝑡ₘᵢₙ,𝑡ᵒᶠᶠₘₐₓ,  𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)

    𝑡ᵒⁿ  = max(𝛥𝑡ₘᵢₙ, 10⋅𝑡ᵣₑₛ)
    𝑡ᵒᶠᶠ = min(       10000⋅𝑡ᵣₑₛ, 𝑡ᵒᶠᶠₘₐₓ)
    evΩ = Evolution_Ω( 𝑡ᵒⁿ , 𝑡ᵒᶠᶠ
                      ; 𝛥 = (0//1)/μs, ε, hw )

	scatter( 𝛺 -> let ψ = copy(ψ₀)
                      evf(𝛺, evΩ ; ϕ,R,ψ)
                  end,
             -𝛺ₘₐₓ: 7𝛺ᵣₑₛ :+𝛺ₘₐₓ
			 ;
             label="",
			 markersize=0.5, markerstrokewidth=0,
			 xaxis="𝛺")
end

plotΔ =
let
	(; 𝑡ᵣₑₛ,𝛥𝑡ₘᵢₙ,𝑡ᵒᶠᶠₘₐₓ,  𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)

    𝑡ᵒⁿ  = max(𝛥𝑡ₘᵢₙ, 10⋅𝑡ᵣₑₛ)
    𝑡ᵒᶠᶠ = min(       10000⋅𝑡ᵣₑₛ, 𝑡ᵒᶠᶠₘₐₓ)
    evΔ = Evolution_Δ( 𝑡ᵒⁿ , 𝑡ᵒᶠᶠ
                      ; 𝛺 = -𝛺ₘₐₓ/100, ε, hw )
	scatter( 𝛥 -> let ψ = copy(ψ₀)
                      evf(𝛥, evΔ ; ϕ,R,ψ)
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
	(; 𝑡ᵣₑₛ,𝛥𝑡ₘᵢₙ,𝑡ᵒᶠᶠₘₐₓ,  𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(hw)

    𝑡ᵒⁿ  = max(𝛥𝑡ₘᵢₙ, 10⋅𝑡ᵣₑₛ)
    𝑡ᵒᶠᶠ = min(       10000⋅𝑡ᵣₑₛ, 𝑡ᵒᶠᶠₘₐₓ)
    evΩ  = Evolution_Ω( 𝑡ᵒⁿ , 𝑡ᵒᶠᶠ
                        ; 𝛥 = (0//1)/μs, ε=1e-4, hw )

    Set_of_𝛺s = (  -𝛺ₘₐₓ: 𝛺ᵣₑₛ :+𝛺ₘₐₓ  )
    N = length(Set_of_𝛺s)
    f⃗ = [ let ψ = copy(ψ₀)
              evf(𝛺, evΩ ; ϕ,R,ψ)
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


#EOF
