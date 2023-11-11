using SNPlots
using Test

@testset "SNPlots.jl" begin
    # Write your tests here.
    @test SNPlots.greet_SNPlots() == "Hello I am the SNPlots package!"
end
