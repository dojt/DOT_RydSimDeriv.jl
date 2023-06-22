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

# using Logging
#
# using Unitful
# using Unitful: μs


using DOT_RydSimDeriv




using JSON # Only for ignoring by JET

@testset verbose=true "Testing DOT_RydSimDeriv.jl" begin

    #
    # Basic JET-based package test only:

    test_package(DOT_RydSimDeriv, ignored_modules=(AnyFrameModule(JSON.Parser),) )

end


#  @testset "A broken test:" begin
#      @test DOODELDIDOO skip=true
#  end

#runtests.jl
#EOF
