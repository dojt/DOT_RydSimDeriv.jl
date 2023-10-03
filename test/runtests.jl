########################################################################
#                                                                      #
# DOT_RydSim/test/runtests.jl                                          #
#                                                                      #
# (c) Dirk Oliver Theis 2023                                           #
#                                                                      #
# License:                                                             #
#                                                                      #
#             Apache 2.0                                               #
#                                                                      #
########################################################################

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0. Packages & Helpers

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0.1. Packages

using DOT_RydSimDeriv

using Test
using JET

# using Logging

using Unitful
using Unitful: μs

using LinearAlgebra: Hermitian, Diagonal, eigvecs, normalize
using Random: shuffle!

using DOT_NiceMath
using DOT_NiceMath.NumbersF64

using DOT_RydSim: is_δrounded


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.1. Helpers

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

function make_Pauli!( H ::Hermitian{ℂ,Matrix{ℂ}} ) :: Nothing
    let d  = size(H) |> first,
	U  = eigvecs(H),
        𝜆⃗  = [  ℝ( (-1)^j  )    for j=1:d  ]
        shuffle!( 𝜆⃗ )
        H .= Hermitian( U⋅Diagonal(𝜆⃗)⋅U' )
    end
    nothing;
end

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1. JET.jl

using JSON # Only for ignoring by JET

@testset verbose=true "Testing DOT_RydSimDeriv.jl" begin


    @testset verbose=true "JET.jl package test" begin
        #
        # Basic JET-based package test only:

        test_package(DOT_RydSimDeriv,
                     ignored_modules=(
                         AnyFrameModule(JSON.Parser) ,
                         AnyFrameModule(Base)  # This is the most idiotic line in the history of computer programs...
                                               # Why not declare vars?!  I hate millenials!
                     ) )
    end

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2. Misc Tests


# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.2. 𝛿𝑡ᵉᶠᶠ

    @testset verbose=true "Effectiv time" begin
        hw = load_hw()

        (  ; 𝑡ᵣₑₛ,𝛥𝑡ₘᵢₙ,𝑡ᵒᶠᶠₘₐₓ,
           𝛿𝑡ᵉᶠᶠₘᵢₙ𝛺,𝛿𝑡ᵉᶠᶠₘₐₓ𝛺, 𝛿𝑡ᵉᶠᶠₘᵢₙ𝛥,𝛿𝑡ᵉᶠᶠₘₐₓ𝛥,
           𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ                   ) = get_hw_data(hw)

        @test is_δrounded(𝛿𝑡ᵉᶠᶠₘᵢₙ𝛺 ; 𝛿=𝑡ᵣₑₛ)
        @test is_δrounded(𝛿𝑡ᵉᶠᶠₘₐₓ𝛺 ; 𝛿=𝑡ᵣₑₛ)
        @test is_δrounded(𝛿𝑡ᵉᶠᶠₘᵢₙ𝛥 ; 𝛿=𝑡ᵣₑₛ)
        @test is_δrounded(𝛿𝑡ᵉᶠᶠₘₐₓ𝛥 ; 𝛿=𝑡ᵣₑₛ)

        @test 𝛿𝑡ᵉᶠᶠₘᵢₙ𝛺 < 𝛿𝑡ᵉᶠᶠₘₐₓ𝛺
        @test 𝛿𝑡ᵉᶠᶠₘᵢₙ𝛥 < 𝛿𝑡ᵉᶠᶠₘₐₓ𝛥
    end

    @testset verbose=true "Make some Ω-evolutions" begin
        hw = load_hw()

        (  ; 𝑡ᵣₑₛ,𝛥𝑡ₘᵢₙ,𝑡ᵒᶠᶠₘₐₓ,
           𝛿𝑡ᵉᶠᶠₘᵢₙ𝛺,𝛿𝑡ᵉᶠᶠₘₐₓ𝛺, 𝛿𝑡ᵉᶠᶠₘᵢₙ𝛥,𝛿𝑡ᵉᶠᶠₘₐₓ𝛥,
           𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ                    ) = get_hw_data(hw)
        (  ; 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
           𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤            ) = hw

        N_ATOMS  = 3
        R_STDDEV = 64
        ε        = 1e-3

        ψ₀ = randn(ℂ,2^N_ATOMS) |> normalize
        R  = let A = randn(ℂ,2^N_ATOMS,2^N_ATOMS) ; Hermitian( (A+A')⋅R_STDDEV/2 ) end
        𝚷  = let A = randn(ℂ,2^N_ATOMS,2^N_ATOMS) ; Hermitian( (A+A')         /2 ) end
        make_Pauli!(𝚷)

        for _1iter = 1:1_000
            𝑡ᵒⁿ   = rand_range( (0//1)μs : 𝑡ᵣₑₛ : 3𝛥𝑡ₘᵢₙ )     ;  𝑡ᵒⁿ  > 𝛥𝑡ₘᵢₙ  ||  ( 𝑡ᵒⁿ  = 𝛥𝑡ₘᵢₙ )
            𝑡ᵒᶠᶠ  = rand_range( 𝑡ᵒⁿ+𝛥𝑡ₘᵢₙ : 𝑡ᵣₑₛ : 𝑡ᵒᶠᶠₘₐₓ )   ;  𝑡ᵒᶠᶠ > 𝛥𝑡ₘᵢₙ  ||  ( 𝑡ᵒᶠᶠ = 𝛥𝑡ₘᵢₙ )
            𝛿𝑡ᵉᶠᶠ = rand_range( 𝛿𝑡ᵉᶠᶠₘᵢₙ𝛺 : 𝑡ᵣₑₛ : 𝛿𝑡ᵉᶠᶠₘₐₓ𝛺/10 )
            𝛥     = rand_range( (0//1)/μs : 𝛥ᵣₑₛ : 𝛥ₘₐₓ )
            𝑇     = max( 𝑡ᵒⁿ + 𝛿𝑡ᵉᶠᶠ + 𝛺ₘₐₓ⋅(1/𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤 + 1/𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤)/2 ,
                         𝑡ᵒᶠᶠ + 𝛥/𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤                                  )

            evΩ = Evolution_Ω( 𝑡ᵒⁿ , 𝛿𝑡ᵉᶠᶠ
                               ; 𝛥, Δ_𝑡ᵒᶠᶠ=𝑡ᵒᶠᶠ,
                               ε, hw, 𝑇)
            @test true

            for _2iter = 1:100
                𝛺 = rand_range( (0//1)/μs : 𝛺ᵣₑₛ : 𝛺ₘₐₓ )
                ψ = copy(ψ₀)
                evf(𝛺, evΩ ; 𝚷,R,ψ)
                @test true
            end
        end
    end

    @testset verbose=true "Make some Δ-evolutions" begin
        hw = load_hw()

        (  ; 𝑡ᵣₑₛ,𝛥𝑡ₘᵢₙ,𝑡ᵒᶠᶠₘₐₓ,
           𝛿𝑡ᵉᶠᶠₘᵢₙ𝛺,𝛿𝑡ᵉᶠᶠₘₐₓ𝛺, 𝛿𝑡ᵉᶠᶠₘᵢₙ𝛥,𝛿𝑡ᵉᶠᶠₘₐₓ𝛥,
           𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ                    ) = get_hw_data(hw)
        (  ; 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
           𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤            ) = hw

        N_ATOMS  = 3
        R_STDDEV = 64
        ε        = 1e-3

        ψ₀ = randn(ℂ,2^N_ATOMS) |> normalize
        R  = let A = randn(ℂ,2^N_ATOMS,2^N_ATOMS) ; Hermitian( (A+A')⋅R_STDDEV/2 ) end
        𝚷  = let A = randn(ℂ,2^N_ATOMS,2^N_ATOMS) ; Hermitian( (A+A')         /2 ) end
        make_Pauli!(𝚷)

        for _1iter = 1:1_000
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

# ———————————————————————————————————————————————————————————————————————————————————————————————————
end #^ testset all of it

#  @testset "A broken test:" begin
#      @test DOODELDIDOO skip=true
#  end

#runtests.jl
#EOF
