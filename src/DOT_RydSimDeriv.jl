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
#        3.3.  EVF — Expectation Value Function                                                               |
#                                                                                                             |
#        3.4.  Fourier spectrum bound                                                                         |
#                                                                                                             |
#                                                                                                             |
#    4.  Shift ruling                                                                                         |
#                                                                                                             |
#        4.1.  EVF-eval based (non-physical)                                                                  |
#                                                                                                             |
#              4.1.a. Type `Shift_Rule{PType_*}`                                                              |
#              4.1.b. Callables                                                                               |
#              4.1.c. Instances                                                                               |
#                     • Symmetric Difference Quotient                                                         |
#                                                                                                             |
#        4.2.  Shot-based                                                                                     |
#                                                                                                             |
#                                                                                                             |
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

* Functions [`get_hw_data`](@ref)`()` (returns struct [`HW_Data`](@ref)) and
  [`get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺`](@ref)`()` — extract hw data

### Defining and evaluation EVFs

EVF = Expectation Value Function

* Types [`Evolution_Ω`](@ref), [`Evolution_Δ`](@ref) with constructors and callables;
  the callables (almost) compute the expectation value function
* Function [`evf`](@ref)`()` — based on the callable for the given evolution object.
* Function [`λ`](@ref)`()`   — approx. lower bound on wavelength in the Fourier spectrum.

### Shift rules — EVF-eval based

* Struct [`Shift_Rule`](@ref)`{PType_`\\*`}` where "\\*" is one of "Ω" or "Δ"
  * Objects are callable.
* Symmetric Difference Quotient:
  * Function [`make_SymDiffQuot`](@ref)`() ::Shift_Rule` — make SR with given ``\varepsilon``
  * Helper fn [`get_𝑥ₘₐₓ_SymDiffQuot`](@ref)`()`
"""
module DOT_RydSimDeriv

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 1.1. Exports
export load_hw,
       HW_Data, get_hw_data, get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺
export Evolution_Ω, Evolution_Δ,
       evf,
       λ
export Shift_Rule, PType_Ω, PType_Δ



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


using Unitful
using Unitful:       μs, ustrip
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

@doc raw"""
Struct `HW_Data`
Abstraction layer for info about HW needed for parameter arithmetic.  The fields are of
unitful rational number types:
* `𝛺ₘₐₓ`, `𝛺ᵣₑₛ`;
* `𝛥ₘₐₓ` `𝛥ᵣₑₛ`;
* `𝑡ᵣₑₛ`;
* `𝛥𝑡ₘᵢₙ`          — smallest positive time
* `𝑡ₘₐₓ`           — max total evolution time
* `𝑡ᵈᵒʷⁿ`          — time needed between 𝑡ᵒᶠᶠ and EOEv to allow for full range of 𝛺 and 𝛥.
* `𝑡ᵒᶠᶠₘₐₓ`        — largest switch-off time which allows full range of 𝛺 and 𝛥
* `𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ`    — smallest duration ``t^{\text{off}}-t^{\text{on}}`` which allows full range
  of 𝛺 and 𝛥
"""
@kwdef struct HW_Data
    𝛺ₘₐₓ        ::Rad_per_μs_t{ℚ}
    𝛺ᵣₑₛ        ::Rad_per_μs_t{ℚ}
    𝛥ₘₐₓ        ::Rad_per_μs_t{ℚ}
    𝛥ᵣₑₛ        ::Rad_per_μs_t{ℚ}
    𝑡ᵈᵒʷⁿ       ::μs_t{ℚ}
    𝑡ᵒᶠᶠₘₐₓ     ::μs_t{ℚ}
    𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ ::μs_t{ℚ}
    𝑡ᵣₑₛ        ::μs_t{ℚ}
    𝛥𝑡ₘᵢₙ       ::μs_t{ℚ}
    𝑡ₘₐₓ        ::μs_t{ℚ}
end

@doc raw"""
Function `get_hw_data(::HW_Descr) ::HW_Data`

"""
function get_hw_data(hw ::HW_Descr{ℚ}) ::HW_Data

    ( ; 𝛺ₘₐₓ, 𝛺ᵣₑₛ, 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤, φᵣₑₛ,
      𝛥ₘₐₓ, 𝛥ᵣₑₛ, 𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
      𝑡ₘₐₓ, 𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ                               ) = hw

    𝛺_𝑢𝑝𝑡𝑖𝑚𝑒 = 𝛺ₘₐₓ / 𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤 ; 𝛺_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒 = 𝛺ₘₐₓ / 𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤
    𝛥_𝑢𝑝𝑡𝑖𝑚𝑒 = 𝛥ₘₐₓ / 𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤 ; 𝛥_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒 = 𝛥ₘₐₓ / 𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤

    𝑡ᵈᵒʷⁿ = δround_up(   max(𝛺_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒,
			     𝛥_𝑑𝑜𝑤𝑛𝑡𝑖𝑚𝑒);
			 𝛿=𝑡ᵣₑₛ )

    return HW_Data(;
	    𝛺ₘₐₓ, 𝛺ᵣₑₛ,
	    𝛥ₘₐₓ, 𝛥ᵣₑₛ,
        𝑡ᵈᵒʷⁿ,
	    𝑡ᵒᶠᶠₘₐₓ     = 𝑡ₘₐₓ - 𝑡ᵈᵒʷⁿ,
	    𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ = δround_up(
		                    max( 𝛺_𝑢𝑝𝑡𝑖𝑚𝑒,
					 𝛥_𝑢𝑝𝑡𝑖𝑚𝑒,
					 𝛥𝑡ₘᵢₙ);
		          𝛿=𝑡ᵣₑₛ),
	    𝑡ᵣₑₛ, 𝛥𝑡ₘᵢₙ, 𝑡ₘₐₓ
    )
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
Type `Evolution_Ω`

Subtype of `Evolution_t`.  This doc doc docs the constructor and callable.

## Constructor
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


## Callable
```julia
    (ev::Evolution_Ω)(𝛺 ::Rad_per_μs_t{ℚ}
                      ;
                      ϕ ::Vector{ℂ},
                      R ::Hermitian{ℂ,Matrix{ℂ}},
                      ψ ::Vector{ℂ}              ) ::ℂ
```

!!! warning "Warning: ψ is updated!"

    The argument `ψ` gives the initial state of the evolution.
    After the function returns, the **vector** `ψ` **contains the final state** of the evolution!

Evaluates ``(\phi \mid U_R(\Omega) \psi)``, where ``U_R(\Omega)`` stands for the quantum evolution
with Rabi frequency ``\Omega`` and the detuning given in the evolution object, with the
"Rydberg"-term ``\hbar R`` in the Hamiltonian, i.e.,
```math
H/\hbar = \frac{\Omega}{2} X - \Delta |1\rangle\langle1| + R,
```
where |1⟩ is the Rydberg state vs |0⟩ the ground state.
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
    is_δrounded( 𝑡ᵒⁿ  ;𝛿=hw.𝑡ᵣₑₛ) || throw(ArgumentError("𝑡ᵒⁿ must be multiple of HW 𝑡ᵣₑₛ"))
    is_δrounded( 𝑡ᵒᶠᶠ ;𝛿=hw.𝑡ᵣₑₛ) || throw(ArgumentError("𝑡ᵒᶠᶠ must be multiple of HW 𝑡ᵣₑₛ"))
    is_δrounded( 𝑇    ;𝛿=hw.𝑡ᵣₑₛ) || throw(ArgumentError("𝑇 must be multiple of HW 𝑡ᵣₑₛ"))
    𝑡ᵒⁿ + 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ ≤ 𝑡ᵒᶠᶠ      || throw(ArgumentError("𝑡ᵒᶠᶠ-𝑡ᵒⁿ must be ≥ 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ"))
    𝑡ᵒⁿ + 𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺 ≤ 𝑡ᵒᶠᶠ      || throw(ArgumentError("𝑡ᵒᶠᶠ-𝑡ᵒⁿ must be ≥ 𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺"))
    𝑡ᵒᶠᶠ ≤ 𝑡ᵒᶠᶠₘₐₓ                || throw(ArgumentError("Need 𝑡ᵒᶠᶠ ≤ 𝑡ᵒᶠᶠₘₐₓ"))
    𝑇 ≤ 𝑡ₘₐₓ                      || throw(ArgumentError("Need 𝑇 ≤ 𝑡ₘₐₓ"))
    𝑡ᵒᶠᶠ + 𝑡ᵈᵒʷⁿ ≤ 𝑇              || throw(ArgumentError("Need 𝑇 ≥ 𝑡ᵒᶠᶠ + 𝑡ᵈᵒʷⁿ"))

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
								;   ev.hw.𝛺ₘₐₓ, ev.hw.𝛺ᵣₑₛ,
									ev.hw.𝛺_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, ev.hw.𝛺_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
									ev.hw.φᵣₑₛ,
									ev.hw.𝑡ₘₐₓ, ev.hw.𝑡ᵣₑₛ, ev.hw.𝛥𝑡ₘᵢₙ)
	DOT_RydSim._check(pΩ)


	schröd!(  ψ,   ℝ( ev.𝑇 )
			  ;
              𝑡₀ = ℝ( ev.𝑡₀ ),
              Ω  = pΩ,
			  Δ  = ev.pΔ,
			  ε  = ev.ε,
			  R             )

	return ϕ' ⋅ ψ
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
Type `Evolution_Ω`

Subtype of `Evolution_t`.  This doc doc docs the constructor and callable.

## Constructor
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

## Callable
```julia
    (ev::Evolution_Δ)(𝛥 ::Rad_per_μs_t{ℚ}
                      ;
                      ϕ ::Vector{ℂ},
                      R ::Hermitian{ℂ,Matrix{ℂ}},
                      ψ ::Vector{ℂ}              ) ::ℂ
```

!!! warning "Warning: ψ is updated!"

    The argument `ψ` gives the initial state of the evolution.
    After the function returns, the **vector** `ψ` **contains the final state** of the evolution!

Evaluates ``(\phi \mid U_R(\Delta) \psi)``, where ``U_R(\Delta)`` stands for the quantum evolution
with detuning ``\Delta`` and the Rabi frequency given in the evolution object, with the
"Rydberg"-term ``\hbar R`` in the Hamiltonian, i.e.,
```math
H/\hbar = \frac{\Omega}{2} X - \Delta |1\rangle\langle1| + R,
```
where |1⟩ is the Rydberg state vs |0⟩ the ground state.
"""
function Evolution_Δ( 𝑡ᵒⁿ  ::μs_t{ℚ},
                      𝑡ᵒᶠᶠ ::μs_t{ℚ}
                      ;
                      𝛺   ::Rad_per_μs_t,
                      ε   ::ℝ,
                      hw  ::HW_Descr,
                      𝑇   ::μs_t{ℚ}     =
                          𝑡ᵒᶠᶠ +
                          get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺(hw;𝛺,𝛥=get_hw_data(hw).𝛥ₘₐₓ) +
                          get_hw_data(hw).𝑡ᵈᵒʷⁿ                            )   ::Evolution_Δ

	(; 𝛺ₘₐₓ,𝛺ᵣₑₛ, 𝛥ₘₐₓ,𝛥ᵣₑₛ, 𝑡ᵈᵒʷⁿ,𝑡ᵒᶠᶠₘₐₓ,𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ, 𝑡ᵣₑₛ,𝑡ₘₐₓ) =
        get_hw_data(hw)
    𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺 =
        get_hw_𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺(hw;𝛺, 𝛥=𝛥ₘₐₓ)

    #
    # Check args
    #
    is_δrounded( 𝑡ᵒⁿ  ;𝛿=hw.𝑡ᵣₑₛ)   || throw(ArgumentError("𝑡ᵒⁿ must be multiple of HW 𝑡ᵣₑₛ"))
    is_δrounded( 𝑡ᵒᶠᶠ ;𝛿=hw.𝑡ᵣₑₛ)   || throw(ArgumentError("𝑡ᵒᶠᶠ must be multiple of HW 𝑡ᵣₑₛ"))
    is_δrounded( 𝑇    ;𝛿=hw.𝑡ᵣₑₛ)   || throw(ArgumentError("𝑇 must be multiple of HW 𝑡ᵣₑₛ"))
    𝑡ᵒⁿ + 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ ≤ 𝑡ᵒᶠᶠ        || throw(ArgumentError("𝑡ᵒᶠᶠ-𝑡ᵒⁿ must be ≥ 𝑡ᵒⁿ_𝑡ᵒᶠᶠₘᵢₙ"))
    # ;
    # ;
    𝑡ᵒᶠᶠ ≤ 𝑡ᵒᶠᶠₘₐₓ                  || throw(ArgumentError("Need 𝑡ᵒᶠᶠ ≤ 𝑡ᵒᶠᶠₘₐₓ"))
    𝑇 ≤ 𝑡ₘₐₓ                        || throw(ArgumentError("Need 𝑇 ≤ 𝑡ₘₐₓ"))
    𝑡ᵒᶠᶠ + 𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺  + 𝑡ᵈᵒʷⁿ ≤ 𝑇 || throw(ArgumentError("Need 𝑇 ≥ 𝑡ᵒᶠᶠ + 𝑡ᵒᶠᶠ⁻ᵈⁱᶠᶠ𝛥𝛺 + 𝑡ᵈᵒʷⁿ"))

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
                𝑡₀     = min(𝑡ᵒⁿ, Ω_𝑡ᵒⁿ),
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
							  ;   ev.hw.𝛥ₘₐₓ, ev.hw.𝛥ᵣₑₛ,
							      ev.hw.𝛥_𝑚𝑎𝑥_𝑢𝑝𝑠𝑙𝑒𝑤, ev.hw.𝛥_𝑚𝑎𝑥_𝑑𝑜𝑤𝑛𝑠𝑙𝑒𝑤,
							      ev.hw.𝑡ₘₐₓ, ev.hw.𝑡ᵣₑₛ, ev.hw.𝛥𝑡ₘᵢₙ)
	DOT_RydSim._check(pΔ)


	schröd!(  ψ, ℝ(ev.𝑇)
			  ;
              𝑡₀ = ev.𝑡₀,
			  Δ  = pΔ,
              Ω  = ev.pΩ,
			  ε  = ev.ε,
			  R             )

	return ϕ' ⋅ ψ
end

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 3.3. EVF

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
Calls the callable of the given Evolution object, `ev`, with initial state ψ and observable
```math
    2|\phi\rangle\langle\phi| - 1
```

!!! warning "Warning: ψ is updated!"

    The argument `ψ` gives the initial state of the evolution.
    After the function returns, the **vector** `ψ` **contains the final state** of the evolution!
"""
function evf(𝑥  ::Rad_per_μs_t{ℚ},
             ev ::EVO
             ;
             ϕ ::Vector{ℂ},
             R ::Hermitian{ℂ,Matrix{ℂ}},
             ψ ::Vector{ℂ},
             kwargs...                  ) ::ℝ   where{EVO<:Evolution_t}

    2⋅abs²( ev(𝑥; ϕ,R,ψ, kwargs...) )   - 1
end

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 3.4. Fourier band bound

@doc raw"""
Function `λ( ev ::Evolution )`

Returns a simple approximate lower bound to the wavelengths occuring in the Fourier spectrum.
"""
λ(ev ::Evolution_Ω) =   2π/ustrip(μs, ev.Ω_𝑡ᵒᶠᶠ - ev.Ω_𝑡ᵒⁿ )
λ(ev ::Evolution_Δ) =   2π/ustrip(μs, ev.Δ_𝑡ᵒᶠᶠ - ev.Δ_𝑡ᵒⁿ )


# ******************************************************************************************************************************
# ——————————————————————————————————————————————————————————————————————————————————————————————————— 4. Shift ruling

# ——————————————————————————————————————————————————————————————————————————————————————————————————— 4.1. EVF-eval based

# -      -      -      -      -      -      -      -      -      -      -      -      -      -      - 4.1.a. Type `Shift_Rule`

abstract type ParameterType_t end
struct PType_Ω <: ParameterType_t end
struct PType_Δ <: ParameterType_t end


@doc raw"""
Struct `Shift_Rule{PType}`

Data to define
```math
    \sum_{j=1}^m  a_j \cdot f(x - s_j)
```

The type parameter `PType` can be one of: `PType_Ω`, `PType_Δ`.

## Fields

* `𝑥 ::Rad_per_μs_t{ℚ}`           — parameter value where shift rule is anchored
* `𝑠 ::Vector{ Rad_per_μs_t{ℚ} }` — list of shifts from 𝑥
* `a ::Vector{ ℝ }`               — coefficients of the shifts

## Constructor and callable


  * Use the keyword-argument constructor:

    ```julia
        Shift_Rule{PType}(; 𝑥, 𝑠, a ) ::Shift_Rule{PType}  where{PType}
    ```

  * The callables

    ```julia
        function (sr::Shift_Rule{PType_Ω})(ev ::Evolution_Ω ; ϕ,R,ψ) ::ℝ
        function (sr::Shift_Rule{PType_Δ})(ev ::Evolution_Δ ; ϕ,R,ψ) ::ℝ
    ```

    compute the shift-rule by evaluating the expectation-value function [`evf`](@ref)`()`.

    Unlike `evf()`, the **callables do not modify the vector `ψ`.**
"""
struct Shift_Rule{PType}
    𝑥 ::Rad_per_μs_t{ℚ}
    𝑠 ::Vector{ Rad_per_μs_t{ℚ} }
    a ::Vector{ ℝ }

    #
    # Use kw-arg to try to ensure that D.A.U. doesn't call this by accident:
    #
    Shift_Rule{PType}(𝑥,𝑠,a;_checking::Bool) where{PType} = ( @assert _checking ; new(𝑥,𝑠,a) )
end

#                                                                                                   # `check_throw()`
@doc raw"""
Function
```julia
    check_throw(sr  ::Shift_Rule{PType},
                hwd ::HW_Data           ) ::Nothing  where{PType <: Union{PType_Ω,PType_Δ}}
```

Checks if the shift rule is conform with the hardware.  If a problem is found, an exception is
*thrown*; otherwise `nothing` is returned.
"""
function check_throw(sr  ::Shift_Rule{PType},
                     hwd ::HW_Data           ) ::Nothing   where{PType<:Union{PType_Ω,PType_Δ}}

    m =  length(sr.𝑠)
    m == length(sr.a) ||  throw(ArgumentError("Lengths of vector `𝑠` ($(m)) and \
                                               `a` ($(length(sr.a))) differ."))
#    m ≥ 1             ||  throw(ArgumentError("Empty shift rule = silly")

    #
    # Check rounding and bounds
    #

    if     PType===PType_Ω      𝛿 = hwd.𝛺ᵣₑₛ ; 𝑥ₘₐₓ = hwd.𝛺ₘₐₓ
    elseif PType===PType_Δ      𝛿 = hwd.𝛥ᵣₑₛ ; 𝑥ₘₐₓ = hwd.𝛥ₘₐₓ
    else                        error("How did you manage to get here?!??")
    end

    (;𝑥) = sr

    is_δrounded(𝑥;𝛿)                          || throw(ArgumentError(
                              "`𝑥` is not aligned to HW resolution."))

    all(  is_δrounded(𝑠;𝛿)    for 𝑠 ∈ sr.𝑠  ) || throw(ArgumentError(
                              "Not all shifts `𝑠` are aligned to HW resolution."))

    all(  -𝑥ₘₐₓ ≤ 𝑥-𝑠 ≤ 𝑥ₘₐₓ  for 𝑠 ∈ sr.𝑠  ) || throw(ArgumentError(
                              "Not all shifts land in the HW parameter range."))

    -𝑥ₘₐₓ ≤ 𝑥 ≤ 𝑥ₘₐₓ                          || throw(ArgumentError(
                              "User is a fucking idiot."))

    return nothing   # All okay!
end #^ check_throw()



function Shift_Rule{PType}(;
                           𝑥   ::Rad_per_μs_t{ℚ},                                                   # Constructor for Shift_Rule 
                           𝑠   ::Vector{ Rad_per_μs_t{ℚ} },
                           a   ::Vector{ ℝ },
                           hwd ::HW_Data
                           ) ::Shift_Rule{PType}  where{PType<:Union{PType_Ω,PType_Δ}}

    sr = Shift_Rule{PType}(𝑥,𝑠,a
                           ; _checking=true)
    check_throw(sr,hwd)

    return sr
end

# -      -      -      -      -      -      -      -      -      -      -      -      -      -      - 4.1.b. Callables

(sr::Shift_Rule{PType_Ω})(ev ::Evolution_Ω ; ϕ,R,ψ) ::ℝ =
    let 𝑥      = sr.𝑥,
        ψᶜᵒᵖʸ  = similar(ψ),
        f(𝑢)   = evf(𝑢, ev ; ϕ,R, ψ=(ψᶜᵒᵖʸ .= ψ))

        sum(   a⋅f(𝑥-𝑠)   for (a,𝑠) ∈ zip( sr.a, sr.𝑠 )   )
    end

(sr::Shift_Rule{PType_Δ})(ev ::Evolution_Δ ; ϕ,R,ψ) ::ℝ =
    let 𝑥      = sr.𝑥,
        ψᶜᵒᵖʸ  = similar(ψ),
        f(𝑢)   = evf(𝑢, ev ; ϕ,R, ψ=(ψᶜᵒᵖʸ .= ψ))

        sum(   a⋅f(𝑥-𝑠)   for (a,𝑠) ∈ zip( sr.a, sr.𝑠 )   )
    end

# -      -      -      -      -      -      -      -      -      -      -      -      -      -      - • Symmetric Difference Quotient

@doc raw"""

"""
function get_𝑥ₘₐₓ_SymDiffQuot(::  Type{PType_Ω},
                              ;
                              n   ::Int,
                              hwd ::HW_Data     ) ::Rad_per_μs_t{ℚ}
    𝜖 = n⋅hwd.𝛺ᵣₑₛ
    return hwd.𝛺ₘₐₓ - 𝜖
end


@doc raw"""
"""
function make_SymDiffQuot(::  Type{PType_Ω},
                          ;
                          𝛺   ::Rad_per_μs_t{ℚ},
                          n   ::Int,
                          hwd ::HW_Data        ) ::Shift_Rule{PType_Ω}
    n ≥ 1 || throw(ArgumentError("Need n ≥ 1"))

    𝜀         = n⋅hwd.𝛺ᵣₑₛ
    let 𝑥ₘₐₓ = get_𝑥ₘₐₓ_SymDiffQuot(PType_Ω;n,hwd),
        𝛿    = hwd.𝛺ᵣₑₛ

        -𝑥ₘₐₓ    ≤ abs(𝛺) ≤    +𝑥ₘₐₓ   || throw(ArgumentError("𝛺 out hw range."))
        -hwd.𝛺ₘₐₓ ≤ abs(𝛺) ≤ +hwd.𝛺ₘₐₓ || throw(ArgumentError("𝛺 out of shift range."))
        -hwd.𝛺ₘₐₓ ≤ 𝛺-𝜀 &&
            𝛺+𝜀 ≤           +hwd.𝛺ₘₐₓ  || throw(ArgumentError("Paranoid for a 𝒓𝒆𝒂𝒔𝒐𝒏!!"))

        is_δrounded(𝛺     ; 𝛿)   || throw(ArgumentError("𝛺 not aligned with HW resolution."))
        is_δrounded(𝛺 - 𝜀 ; 𝛿)   || throw(ArgumentError("𝛺-𝜀 not aligned with HW resolution."))
        is_δrounded(𝛺 + 𝜀 ; 𝛿)   || throw(ArgumentError("𝛺+𝜀 not aligned with HW resolution."))
    end

    α = ℝ(  1 / ustrip(u"μs^(-1)", 2𝜀)  )
    return Shift_Rule{PType_Ω}( ;  𝑥  = 𝛺,
                                𝑠     = Rad_per_μs_t{ℚ}[ -𝜀 , +𝜀 ],
                                a     =                [ +α , -α ]  )
end #^ make_SymDiffQuot




end #^ module DOT_RydSimDeriv
# *************************************************************************************************** EOF
#EOF
