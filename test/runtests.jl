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

    @test isequal(freqs, Float32.([0.0  0.5  1.0  1.0   NaN;
                                   0.5  0.0  1.0  0.25  0.0;
                                   0.5  0.16666667f0  1.0  0.25  0.0]))
    # note that isequal() treats NaN elements as equal,
    # whereas `==` does not.

    @test isequal(sampleSizes, Int16.([2  2  2  1  0
                                       1  2  2  2  1
                                       1  3  3  2  1])) 
                           
    # test getSitePi()

    sitePi = getSitePi(freqs, sampleSizes)

    @test isequal(sitePi, Float64.([0.0  2/3  0.0  0.0  NaN;
                           1.0  0.0       0.0  0.5    0.0;
                           1.0  0.33333334922790525  0.0  0.5    0.0]))
    # the strange number close to 1/3 in the matrix above is a result of freqs being Float32 and sitePi being Float64  
    
    # test getRegionPi()

    meanPi = getRegionPi(sitePi)

    @test isequal(meanPi, [0.16666666666666666;
                            0.3;
                            0.36666666984558105])

    # test getPairwiseNames() and getDxy()

    Dxy, pairwiseNames = getDxy(freqs, basicData_groupsToCalc)

    @test isequal(pairwiseNames, ["pop1_pop2";
                                  "pop1_pop3";
                                  "pop2_pop3"])

    @test isequal(Dxy, Float32.([0.5  0.5       0.0  0.75   NaN
                                0.5  0.5       0.0  0.75   NaN
                                0.5  0.16666667f0  0.0  0.375    0.0]))

    # test getRegionDxy()

    meanDxy = getRegionDxy(Dxy)

    @test isequal(meanDxy, [0.4375;
                            0.4375;
                            0.2083333432674408])

    # test getFst()

    Fst, FstNum, FstDen, pairwiseNames2 = getFst(freqs, sampleSizes, basicData_groupsToCalc)

    @test isequal(Fst, Float32.([0.38461542f0   0.33333334f0   NaN   0.5294118f0  NaN;
                                 0.38461542f0   0.032967024f0   NaN   0.5294118f0  NaN;
                                -1.0       -0.08108107f0  NaN  -0.33333334f0  NaN]))

    @test isequal(FstNum, Float32.([0.06944445f0   0.083333336f0  0.0   0.1875f0  NaN;
                                    0.06944445f0   0.0074999984f0     0.0   0.1875f0  NaN;
                                    -0.25       -0.007499999f0     0.0  -0.0625    0.0]))

    @test isequal(FstDen, Float32.([0.18055555f0  0.25    0.0  0.35416666f0  NaN;
                                    0.18055555f0  0.2275  0.0  0.35416666f0  NaN;
                                    0.25      0.0925  0.0  0.1875      0.0]))

    @test isequal(pairwiseNames2, ["pop1_pop2";
                                    "pop1_pop3";
                                    "pop2_pop3"])

end
