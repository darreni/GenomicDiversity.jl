using SNPlots
using Test

@testset "SNPlots.jl greeting" begin
    # Write your tests here.
    @test SNPlots.greet_SNPlots() == "Hello I am the SNPlots package!"
end

@testset "SNPlots.jl basic data" begin
    # This is to test a simple dataset with 3 pops and 
    # all 3 genotypes and two kinds of missing data.
    basicData_genoData = [0 1 2 -1 missing;
                    0 1 2 2 -1;
                    1 0 2 0 0;
                    -1 0 2 1 -1;
                    missing 1 2 -1 missing;
                    1 0 2 0 0;
                    -1 0 2 1 -1]
    basicData_indGroup = ["pop1",
                "pop1",
                "pop2",
                "pop2",
                "pop3",
                "pop3",
                "pop3"]
    basicData_groupsToCalc = ["pop1", "pop2", "pop3"]

    # test getFreqsAndSampleSizes()

    freqs, sampleSizes = getFreqsAndSampleSizes(basicData_genoData, basicData_indGroup, basicData_groupsToCalc)

    @test isequal(freqs, Float64.([0.0  0.5  1.0  1.0   NaN;
                                   0.5  0.0  1.0  0.25  0.0;
                                   0.5  1/6  1.0  0.25  0.0]))
    # note that isequal() treats NaN elements as equal,
    # whereas `==` does not.

    @test isequal(sampleSizes, Int16.([2  2  2  1  0
                                       1  2  2  2  1
                                       1  3  3  2  1]) 
                           
    # test getSitePi()

    sitePi = getSitePi(freqs, sampleSizes)

    @test isequal(sitePi, Float64.([0.0  2/3  0.0  0.0  NaN;
                           1.0  0.0       0.0  0.5    0.0;
                           1.0  1/3  0.0  0.5    0.0]))
end
