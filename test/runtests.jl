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

using Test
using JET

using Logging

using Unitful
using Unitful: μs


using DOT_RydSimDeriv





@testset verbose=true "Testing DOT_RydSimDeriv.jl" begin
    @test  :donk  skip=true
end

#  @testset "A broken test:" begin
#      @test DOODELDIDOO skip=true
#  end

#runtests.jl
#EOF
