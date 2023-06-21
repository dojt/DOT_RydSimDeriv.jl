########################################################################
#                                                                      #
# DOT_RydSimDeriv.jl                                                   #
#                                                                      #
# (c) Dirk Oliver Theis 2023                                           #
#                                                                      #
# License:                                                             #
#                                                                      #
#             Apache 2.0                                               #
#                                                                      #
########################################################################

# ***************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 0. ToC  +
#                                                                                                             |
#  Table of Contents                                                                                          |
#  -----------------                                                                                          |
#                                                                                                             |
#    1.  Module header & imports                                                                              |
#                                                                                                             |
#        1.1.  Exports                                                                                        |
#        1.2.  Imports                                                                                        |
#                                                                                                             |
#                                                                                                             |
#    2.  Interface to lower level (RydSim)                                                                    |
#                                                                                                             |
#        2.1.  `load_hw()`                                                                                    |
#        2.2.  `get_hw_data()`                                                                                |
#        2.3.  `get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺()`                                                                         |
#                                                                                                             |
#                                                                                                             |
#    3.  Evolutions                                                                                           |
#                                                                                                             |
#        3.1.  Ω Evolution                                                                                    |
#              3.1.a. Ω Constructor                                                                           |
#              3.1.b. Ω Callable                                                                              |
#                                                                                                             |
#        3.2.  Δ Evolution                                                                                    |
#              3.2.a. Δ Constructor                                                                           |
#              3.2.b. Δ Callable                                                                              |
#                                                                                                             |
#                                                                                                             |
#    4.  EVF — Expectation Value Function                                                                     |
#                                                                                                             |
#    5.  ..............                                                                                       |
#                                                                                                             |
#—————————————————————————————————————————————————————————————————————————————————————————————————————————————+

# ******************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1. Module header & imports

"""
Module `DOT_RydSimDeriv`

Build on top of package `DOT_RydSim` to provide functionality to run the type of quantum evolutions
that are used in derivatives based on shift rules such as symmetric difference quotients,
Banchi-Crooks's "stochastic" shift rules, and the "Nyquist" shift rules.

# Exports

### Interface with lower stack

* Function [`load_hw`](@ref)`()` — load hardware configuration from file
* Functions [`get_hw_data`](@ref)`()`, [`get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺`](@ref)`()` — extract hw data

### Defining and evaluation EVFs

EVF = Expectation Value Function

* Types with constructors [`Evolution_Ω`](@ref), [`Evolution_Δ`](@ref).



"""
module DOT_RydSimDeriv

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.1. Exports
export load_hw
export get_hw_data, get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺
export evf_Ω, evf_Ω

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.1. Imports
using DOT_NiceMath
using DOT_NiceMath.NumbersF64

using DOT_RydSim
using DOT_RydSim:
    μs_t,
    Rad_per_μs_t


using DOT_RydSim.HW_Descriptions:
    HW_Descr,
	default_HW_Descr,
	fileread_HW_Descr,
	HW_AWS_QuEra


using Unitful: μs
using LinearAlgebra: Hermitian


# ******************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2. Interface to lower level

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.1. load_hw()
"""
Function
```julia
load_hw(filename = :default
		; Ω_downslew_factor = 1//1,
		  Δ_downslew_factor = 1//1  ) ::HW_Descr
```

The `filename` can be either a string identifying a file name, or the symbol `:default`.

The only file type currently supported is AWS-QuEra's (`HW_AWS_QuEra`).
"""
function load_hw(filename ::String
				;
				Ω_downslew_factor = 1//1,
				Δ_downslew_factor = 1//1              ) ::HW_Descr{ℚ}

	return fileread_HW_Descr(HW_AWS_QuEra
							;   filename,
								ℤ,
								Ω_downslew_factor,
								Δ_downslew_factor )
end

function load_hw(select ::Symbol =:default
				;
				Ω_downslew_factor = 1//1,
				Δ_downslew_factor = 1//1              ) ::HW_Descr{ℚ}

	@assert select == :default  "What?!??"
	return default_HW_Descr(;
							ℤ,
							Ω_downslew_factor,
							Δ_downslew_factor)
end

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.2. get_hw_data()

_NT = @NamedTuple{
				_blah::Nothing,
				𝛺ₘₐₓ        ::Rad_per_μs_t{ℚ},
				𝛺ᵣₑₛ        ::Rad_per_μs_t{ℚ},
				𝛥ₘₐₓ        ::Rad_per_μs_t{ℚ}, 𝛥ᵣₑₛ::Rad_per_μs_t{ℚ},
				𝑡ᵈᵒʷⁿ       ::μs_t{ℚ},
				𝑡ᵒᶠᶠₘₐₓ     ::μs_t{ℚ},
				𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ ::μs_t{ℚ},
				𝑡ᵣₑₛ        ::μs_t{ℚ},
				𝑡ₘₐₓ        ::μs_t{ℚ}
		}

@doc raw"""
Function `get_hw_data(::HW_Descr) ::NamedTuple`

Returns a named tuple with the following fields, all of
unitful rational number types:
* `𝛺ₘₐₓ`, `𝛺ᵣₑₛ`;
* `𝛥ₘₐₓ` `𝛥ᵣₑₛ`;
* `𝑡ᵣₑₛ`;
* `𝑡ₘₐₓ`           — max total evolution time
* `𝑡ᵈᵒʷⁿ`          — time needed between 𝑡ᵒᶠᶠ and EOEv to allow for 
  full range of 𝛺 and 𝛥.
* `𝑡ᵒᶠᶠₘₐₓ`        — largest switch-off time which allows full range of 𝛺 and 𝛥
* `𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ`    — smallest duration ``t^{\text{off}}-t^{\text{on}}``
  which allows full range of 𝛺 and 𝛥
"""
function get_hw_data(hw ::HW_Descr{ℚ}) ::_NT

	( ; 𝛺ₘₐₓ, 𝛺ᵣₑₛ, 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤, φᵣₑₛ,
		𝛥ₘₐₓ, 𝛥ᵣₑₛ, 𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
		𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ                               ) = hw

	𝛺_𝑢𝑝𝑡𝑖𝑚𝑒 = 𝛺ₘₐₓ / 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤 ; 𝛺_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒 = 𝛺ₘₐₓ / 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤
	𝛥_𝑢𝑝𝑡𝑖𝑚𝑒 = 𝛥ₘₐₓ / 𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤 ; 𝛥_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒 = 𝛥ₘₐₓ / 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤

    𝑡ᵈᵒʷⁿ = δround_up(   max(𝛺_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒,
			                 𝛥_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒);
			             𝛿=𝑡ᵣₑₛ )

	return (_blah=nothing,
			𝛺ₘₐₓ, 𝛺ᵣₑₛ,
			𝛥ₘₐₓ, 𝛥ᵣₑₛ,
            𝑡ᵈᵒʷⁿ,
			𝑡ᵒᶠᶠₘₐₓ     = 𝑡ₘₐₓ - 𝑡ᵈᵒʷⁿ,
			𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ = δround_up(
							max( 𝛺_𝑢𝑝𝑡𝑖𝑚𝑒,
								 𝛥_𝑢𝑝𝑡𝑖𝑚𝑒,
								 𝛥𝑡ₘᵢₙ);
							𝛿=𝑡ᵣₑₛ),
			𝑡ᵣₑₛ, 𝑡ₘₐₓ)
end

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 2.3. get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺()

@doc raw"""
Function `get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺(hw::HW_Descr ;  𝛺 =hw.𝛺ₘₐₓ, 𝛥 =hw.𝛥ₘₐₓ) `

𝛺-pulse must end this quantity *later* than 𝛥-pulse in order not to break the RWA with max 𝛺,𝛥
"""
function
get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺( hw ::HW_Descr{ℚ}
					;
					𝛺 :: Rad_per_μs_t{ℚ} = hw.𝛺ₘₐₓ,
					𝛥 :: Rad_per_μs_t{ℚ} = hw.𝛥ₘₐₓ ) ::μs_t{ℚ}

	( ; 𝛺ₘₐₓ, 𝛺ᵣₑₛ, 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤, φᵣₑₛ,
		𝛥ₘₐₓ, 𝛥ᵣₑₛ, 𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
		𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ                               ) = hw

	𝛺 = abs(𝛺)
	𝛥 = abs(𝛥)
	@assert 𝛺 > 0/μs

	𝛺_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒 = 𝛺 / 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤
	𝛥_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒 = 𝛥 / 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤

	return δround_up( max((0//1)μs, 𝛥_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒 - 𝛺_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒 )
					  ; 𝛿=𝑡ᵣₑₛ)
end

# ******************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 3. Evolutions

abstract type Evolution_t end

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 3.1. Ω Evolution

@kwdef struct Evolution_Ω <: Evolution_t
    pΔ     :: Pulse__Δ_BangBang{ℚ}

    𝑡₀     ::μs_t{ℚ}
	Ω_𝑡ᵒⁿ  ::μs_t{ℚ}
    Ω_𝑡ᵒᶠᶠ ::μs_t{ℚ}
	𝑇      ::μs_t{ℚ}

    ε      ::ℝ

	hw     ::HW_Descr
end

# -      -      -      -      -      -      -      -      -      -      -      -      -      -      - 3.1.a. Ω Constructor
@doc raw"""
Constructor
```julia
Evolution_Ω( 𝑡ᵒⁿ  ::μs_t{ℚ},
             𝑡ᵒᶠᶠ ::μs_t{ℚ}
             ;
             𝛥    ::Rad_per_μs_t,
             ε    ::ℝ,
             hw   ::HW_Descr,
             𝑇    ::μs_t{ℚ}     = ...)   ::Evolution_Ω
```

Creates evolution data for variable 𝛺, with fixed 𝛥, starting at time 𝑡ᵒⁿ, and ending at time 𝑇.

A lower bound for the end-time 𝑇 of the evolution is 𝑡ᵒᶠᶠ + 𝑡ᵈᵒʷⁿ; an upper bound is 𝑡ₘₐₓ.

Lower bounds for 𝑡ᵒᶠᶠ-𝑡ᵒⁿ are 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ and 𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺; an upper bound for 𝑡ᵒᶠᶠ is 𝑡ᵒᶠᶠₘₐₓ.

The quantities mentioned above are defined in the named tuple returned by
[`get_hw_data`](@ref)`()`, except for 𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺, which is returned by
[`get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺`](@ref)`()`.
"""
function Evolution_Ω( 𝑡ᵒⁿ  ::μs_t{ℚ},
                      𝑡ᵒᶠᶠ ::μs_t{ℚ}
                      ;
                      𝛥  ::Rad_per_μs_t,
                      ε  ::ℝ,
                      hw ::HW_Descr,
                      𝑇  ::μs_t{ℚ}     =
                          𝑡ᵒᶠᶠ +
                          # ;
                          get_hw_data(hw).𝑡ᵈᵒʷⁿ             )   ::Evolution_Ω

	(; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ, 𝑡ᵈᵒʷⁿ,𝑡ᵒᶠᶠₘₐₓ,𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ, 𝑡ᵣₑₛ,𝑡ₘₐₓ) =
        get_hw_data(hw)
    𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺 =
        get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺(hw;𝛺=𝛺ₘₐₓ, 𝛥)

    #
    # Check args
    #
    𝑡ᵒⁿ  == δround_down( 𝑡ᵒⁿ ;𝛿=hw.𝑡ᵣₑₛ)  || throw(ArgumentError(
                                             "𝑡ᵒⁿ must be multiple of HW 𝑡ᵣₑₛ"))
    𝑡ᵒᶠᶠ == δround_down( 𝑡ᵒᶠᶠ ;𝛿=hw.𝑡ᵣₑₛ) || throw(ArgumentError(
                                             "𝑡ᵒᶠᶠ must be multiple of HW 𝑡ᵣₑₛ"))
    𝑇    == δround_down( 𝑇 ;𝛿=hw.𝑡ᵣₑₛ)    || throw(ArgumentError(
                                             "𝑇 must be multiple of HW 𝑡ᵣₑₛ"))
    𝑡ᵒⁿ + 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ ≤ 𝑡ᵒᶠᶠ              || throw(ArgumentError(
                                             "𝑡ᵒᶠᶠ-𝑡ᵒⁿ must be ≥ 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ"))
    𝑡ᵒⁿ + 𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺 ≤ 𝑡ᵒᶠᶠ              || throw(ArgumentError(
                                             "𝑡ᵒᶠᶠ-𝑡ᵒⁿ must be ≥ 𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺"))
    𝑡ᵒᶠᶠ ≤ 𝑡ᵒᶠᶠₘₐₓ                        || throw(ArgumentError(
                                             "Need 𝑡ᵒᶠᶠ ≤ 𝑡ᵒᶠᶠₘₐₓ"))
    𝑇 ≤ 𝑡ₘₐₓ                              || throw(ArgumentError(
                                             "Need 𝑇 ≤ 𝑡ₘₐₓ"))
    𝑡ᵒᶠᶠ + 𝑡ᵈᵒʷⁿ ≤ 𝑇                      || throw(ArgumentError(
                                             "Need 𝑇 ≥ 𝑡ᵒᶠᶠ + 𝑡ᵈᵒʷⁿ"))
    #
    # Make Δ pulse
    #
	Δ_𝑡ᵒⁿ  = 𝑡ᵒⁿ
    Δ_𝑡ᵒᶠᶠ = 𝑡ᵒᶠᶠ - 𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺

	pΔ = Pulse__Δ_BangBang{ℚ}(Δ_𝑡ᵒⁿ, Δ_𝑡ᵒᶠᶠ, 𝑇, 𝛥
							  ;  hw.𝛥ₘₐₓ, hw.𝛥ᵣₑₛ,
							     hw.𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, hw.𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
                                 # ;
								 hw.𝑡ₘₐₓ, hw.𝑡ᵣₑₛ, hw.𝛥𝑡ₘᵢₙ)
	DOT_RydSim._check(pΔ)

    #
    # Construct it:
    #
    Evolution_Ω(;pΔ,
                Ω_𝑡ᵒⁿ  = 𝑡ᵒⁿ,
                Ω_𝑡ᵒᶠᶠ = 𝑡ᵒᶠᶠ,
                𝑡₀     = min(𝑡ᵒⁿ, Δ_𝑡ᵒⁿ),
                𝑇,
                ε,
                hw)
end #^ Evolution_Ω()

# -      -      -      -      -      -      -      -      -      -      -      -      -      -      - 3.1.b. Ω Callable
function (ev::Evolution_Ω)(𝛺 ::Rad_per_μs_t{ℚ}
                           ;
                           ϕ ::Vector{ℂ},
                           R ::Hermitian{ℂ,Matrix{ℂ}},
                           ψ ::Vector{ℂ}              ) ::ℂ

    @assert length(ϕ) == length(ψ)
    @assert ( length(ϕ) , length(ψ) ) == size(R)


	(; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(ev.hw)
	# (; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ, 𝑡ᵒᶠᶠₘₐₓ,𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ) = get_hw_info(hw)

	pΩ = Pulse__Ω_BangBang{ℚ,ℝ}(ev.Ω_𝑡ᵒⁿ, ev.Ω_𝑡ᵒᶠᶠ, ev.𝑇,
                                𝛺
								;   hw.𝛺ₘₐₓ, hw.𝛺ᵣₑₛ,
									hw.𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, hw.𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
									hw.φᵣₑₛ,
									hw.𝑡ₘₐₓ, hw.𝑡ᵣₑₛ, hw.𝛥𝑡ₘᵢₙ)
	DOT_RydSim._check(pΩ)


	ψᵤₛₑ = copy(ψ)

	schröd!(  ψᵤₛₑ, ℝ(ev.𝑇)
			  ;
              𝑡₀ = ev.𝑡₀,
              Ω  = pΩ,
			  Δ  = ev.pΔ,
			  R,
			  ε )

	return ϕ' ⋅ ψᵤₛₑ
end

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 3.2. Δ Evolution

@kwdef struct Evolution_Δ <: Evolution_t
    pΩ     :: Pulse__Ω_BangBang{ℚ,ℝ}

    𝑡₀     ::μs_t{ℚ}
	Δ_𝑡ᵒⁿ  ::μs_t{ℚ}
    Δ_𝑡ᵒᶠᶠ ::μs_t{ℚ}
	𝑇      ::μs_t{ℚ}

    ε      ::ℝ

	hw     ::HW_Descr
end

# -      -      -      -      -      -      -      -      -      -      -      -      -      -      - 3.2.a. Δ Constructor
@doc raw"""
Constructor
```julia
Evolution_Δ( 𝑡ᵒⁿ  ::μs_t{ℚ},
             𝑡ᵒᶠᶠ ::μs_t{ℚ}
             ;
             𝛺    ::Rad_per_μs_t,
             ε    ::ℝ,
             hw   ::HW_Descr,
             𝑇    ::μs_t{ℚ}     = ...)   ::Evolution_Ω
```

Creates evolution data for variable 𝛥, with fixed 𝛺, starting at time 𝑡ᵒⁿ, and ending at time 𝑇.

A lower bound for the end-time 𝑇 of the evolution is 𝑡ᵒᶠᶠ + 𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺 + 𝑡ᵈᵒʷⁿ; an upper bound is
𝑡ₘₐₓ.

A lower bound for 𝑡ᵒᶠᶠ is 𝑡ᵒⁿ + 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ; an upper bound is 𝑡ᵒᶠᶠₘₐₓ.

The quantities mentioned above are defined in the named tuple returned by
[`get_hw_data`](@ref)`()`, except for 𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺, which is returned by
[`get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺`](@ref)`()`.
"""
function Evolution_Δ( 𝑡ᵒⁿ  ::μs_t{ℚ},
                      𝑡ᵒᶠᶠ ::μs_t{ℚ}
                      ;
                      𝛺   ::Rad_per_μs_t,
                      ε   ::ℝ,
                      hw  ::HW_Descr,
                      𝑇   ::μs_t{ℚ}     =
                          𝑡ᵒᶠᶠ +
                          get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺(hw;𝛺, 𝛥=𝛥ₘₐₓ) +
                          get_hw_data(hw).𝑡ᵈᵒʷⁿ             )   ::Evolution_Δ

	(; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ, 𝑡ᵈᵒʷⁿ,𝑡ᵒᶠᶠₘₐₓ,𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ, 𝑡ᵣₑₛ,𝑡ₘₐₓ) =
        get_hw_data(hw)
    𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺 =
        get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺(hw;𝛺, 𝛥=𝛥ₘₐₓ)

    #
    # Check args
    #
    𝑡ᵒⁿ  == δround_down( 𝑡ᵒⁿ ;𝛿=hw.𝑡ᵣₑₛ)  || throw(ArgumentError(
                                             "𝑡ᵒⁿ must be multiple of HW 𝑡ᵣₑₛ"))
    𝑡ᵒᶠᶠ == δround_down( 𝑡ᵒᶠᶠ ;𝛿=hw.𝑡ᵣₑₛ) || throw(ArgumentError(
                                             "𝑡ᵒᶠᶠ must be multiple of HW 𝑡ᵣₑₛ"))
    𝑇    == δround_down( 𝑇 ;𝛿=hw.𝑡ᵣₑₛ)    || throw(ArgumentError(
                                             "𝑇 must be multiple of HW 𝑡ᵣₑₛ"))
    𝑡ᵒⁿ + 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ ≤ 𝑡ᵒᶠᶠ              || throw(ArgumentError(
                                             "𝑡ᵒᶠᶠ-𝑡ᵒⁿ must be ≥ 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ"))
    # ;
    # ;
    𝑡ᵒᶠᶠ ≤ 𝑡ᵒᶠᶠₘₐₓ                        || throw(ArgumentError(
                                             "Need 𝑡ᵒᶠᶠ ≤ 𝑡ᵒᶠᶠₘₐₓ"))
    𝑇 ≤ 𝑡ₘₐₓ                              || throw(ArgumentError(
                                             "Need 𝑇 ≤ 𝑡ₘₐₓ"))
    𝑡ᵒᶠᶠ + 𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺  + 𝑡ᵈᵒʷⁿ ≤ 𝑇       || throw(ArgumentError(
                                             "Need 𝑇 ≥ 𝑡ᵒᶠᶠ + 𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺 + 𝑡ᵈᵒʷⁿ"))
    #
    # Make Ω pulse
    #
    Ω_𝑡ᵒⁿ  = 𝑡ᵒⁿ
    Ω_𝑡ᵒᶠᶠ = 𝑡ᵒᶠᶠ + get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺(hw;𝛺, 𝛥=𝛥ₘₐₓ)

	pΩ = Pulse__Ω_BangBang{ℚ,ℝ}(Ω_𝑡ᵒⁿ, Ω_𝑡ᵒᶠᶠ, 𝑇, 𝛺
								;   hw.𝛺ₘₐₓ, hw.𝛺ᵣₑₛ,
									hw.𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, hw.𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
									hw.φᵣₑₛ,
									hw.𝑡ₘₐₓ, hw.𝑡ᵣₑₛ, hw.𝛥𝑡ₘᵢₙ)
	DOT_RydSim._check(pΩ)

    #
    # Construct it:
    #
    Evolution_Δ(;pΩ,
                Δ_𝑡ᵒⁿ  = 𝑡ᵒⁿ,
                Δ_𝑡ᵒᶠᶠ = 𝑡ᵒᶠᶠ,
                𝑡₀     = min(𝑡ᵒⁿ, Δ_𝑡ᵒⁿ),
                𝑇,
                ε,
                hw)
end #^ Evolution_Δ()

# -      -      -      -      -      -      -      -      -      -      -      -      -      -      - 3.2.b. Δ Callable
function (ev::Evolution_Δ)(𝛥 ::Rad_per_μs_t{ℚ}
                           ;
                           ϕ ::Vector{ℂ},
                           R ::Hermitian{ℂ,Matrix{ℂ}},
                           ψ ::Vector{ℂ}              ) ::ℂ

    @assert length(ϕ) == length(ψ)
    @assert ( length(ϕ) , length(ψ) ) == size(R)


	(; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ) = get_hw_data(ev.hw)
	# (; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ, 𝑡ᵒᶠᶠₘₐₓ,𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ) = get_hw_info(hw)

	pΔ = Pulse__Δ_BangBang{ℚ}(ev.Δ_𝑡ᵒⁿ, ev.Δ_𝑡ᵒᶠᶠ, ev.𝑇,
                              𝛥
							  ;   hw.𝛥ₘₐₓ, hw.𝛥ᵣₑₛ,
							      hw.𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, hw.𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
							      hw.𝑡ₘₐₓ, hw.𝑡ᵣₑₛ, hw.𝛥𝑡ₘᵢₙ)
	DOT_RydSim._check(pΔ)


	ψᵤₛₑ = copy(ψ)

	schröd!(  ψᵤₛₑ, ℝ(ev.𝑇)
			  ;
              𝑡₀ = ev.𝑡₀,
			  Δ  = pΔ,
              Ω  = ev.pΩ,
			  R,
			  ε )

	return ϕ' ⋅ ψᵤₛₑ
end

# ******************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 4. EVF

@doc raw"""
Function
```julia
    evf(𝑥  ::Rad_per_μs_t{ℚ},
        ev ::EVO
        ;
        ϕ ::Vector{ℂ},
        R ::Hermitian{ℂ,Matrix{ℂ}},
        ψ ::Vector{ℂ}              ) ::ℝ   where{EVO<:Evolution_t}
```
Calls the callable of the given Evolution object, `ev`.
"""
function evf(𝑥  ::Rad_per_μs_t{ℚ},
             ev ::EVO
             ;
             ϕ ::Vector{ℂ},
             R ::Hermitian{ℂ,Matrix{ℂ}},
             ψ ::Vector{ℂ}              ) ::ℝ   where{EVO<:Evolution_t}

        1 - 2⋅abs²( ev(x; ϕ,R,ψ, kwargs...) )
end

#----------------------------------------------------------------------------------------------------#
#                                                                                                    #
#----------------------------------------------------------------------------------------------------#


function evf_Ω(𝛺  ::Rad_per_μs_t
               ;
               𝛥   ::Rad_per_μs_t,
               ϕ  ::Vector{ℂ},
               R  ::Hermitian{ℂ,Matrix{ℂ}},
               ψ  ::Vector{ℂ},
               ε  ::ℝ,
               hw ::HW_Descr              ) ::ℝ

	(; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ, 𝑡ᵒᶠᶠₘₐₓ,𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ, 𝑡ᵣₑₛ,𝑡ₘₐₓ) = get_hw_data(hw)
	(;𝛥𝑡ₘᵢₙ) = hw

	evf(ϕ,R,ψ ; 𝛺, 𝛥,
	 	Ω_𝑡ᵒⁿ=𝛥𝑡ₘᵢₙ, Ω_𝑡ᵒᶠᶠ=𝑡ᵒᶠᶠₘₐₓ,
	 	Δ_𝑡ᵒⁿ=𝛥𝑡ₘᵢₙ, Δ_𝑡ᵒᶠᶠ=𝑡ᵒᶠᶠₘₐₓ-get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺(hw;𝛺=𝛺ₘₐₓ, 𝛥),
	 	𝑇=𝑡ₘₐₓ,
        ε,
	 	hw)
end

function evf_Δ(𝛥  ::Rad_per_μs_t
               ;
               𝛺  ::Rad_per_μs_t,
               ϕ  ::Vector{ℂ},
               R  ::Hermitian{ℂ,Matrix{ℂ}},
               ψ  ::Vector{ℂ},
               ε  ::ℝ,
               hw ::HW_Descr              ) ::ℝ
    (; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ, 𝑡ᵒᶠᶠₘₐₓ,𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ, 𝑡ᵣₑₛ,𝑡ₘₐₓ) = get_hw_data(hw)
	(;𝛥𝑡ₘᵢₙ) = hw

    evf(ϕ,R,ψ ; 𝛺, 𝛥,
	 	Ω_𝑡ᵒⁿ=𝛥𝑡ₘᵢₙ, Ω_𝑡ᵒᶠᶠ=𝑡ᵒᶠᶠₘₐₓ+get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺(hw;𝛺, 𝛥=𝛥ₘₐₓ),
	 	Δ_𝑡ᵒⁿ=𝛥𝑡ₘᵢₙ, Δ_𝑡ᵒᶠᶠ=𝑡ᵒᶠᶠₘₐₓ,
	 	𝑇=𝑡ₘₐₓ,
        ε,
	 	hw)
end


end #^ module DOT_RydSimDeriv
# EOF
