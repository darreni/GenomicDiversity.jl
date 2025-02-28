# This file serves as an example script that I am using during development of GenomicDiversity.jl V0.3.0

# load required packages
# using GenomicDiversity, MultivariateStats, CairoMakie, DataFrames, CSV, DelimitedFiles, Impute

using GenomicDiversity, MultivariateStats, CairoMakie, DataFrames, CSV, DelimitedFiles, Impute

# set data file names
pathName = "demoData/SparrowDemo_data_McCallumetal2024/"
metadataFile = pathName * "SparrowDemo.Fst_groups.txt"
genotypesFile = pathName * "SparrowDemo_genotypes.012"
individualsFile = genotypesFile * ".indv"
positionsFile = genotypesFile * ".pos"

# Now we load the data:

sparrowGenoData = GenomicDiversity.loadGenoData(metadataFile, individualsFile, genotypesFile, positionsFile)


groupsToPlot = ["GCSP","PSWS","GWCS"]

groupColorKey = Dict("GCSP" => "gold",
                    "PSWS" => "red",
                    "GWCS" => "blue")

groupColors = [groupColorKey[i] for i in groupsToPlot]

#= If we want to use Principal Components Analysis (PCA) to visualize a good 2-dimensional representation
of genomic relationships among individuals, we first need to impute missing genotypes.
We will modify the type of matrix so that missing values are represented by `missing` rather than `-1`, 
and then we'll impute: =#

sparrowGenoDataImputed = GenomicDiversity.imputeMissingGenoData(sparrowGenoData)

# Now we are ready to call the `GenomicDiversity.plotPCA()` function to do the PCA and make the plot:

PCAmodel = GenomicDiversity.plotPCA(sparrowGenoDataImputed,
        groupsToPlot, groupColors;
        sampleSet="Zonotrichia sparrows", regionText="whole_genome",
        flip1=false, flip2=false, showPlot=true)
PCAmodel.PCAfig 

# The above shows genomic variation illustrated in 2 dimensions, showing Golden-crowned sparrows in yellow, _pugetensis_ White-crowned Sparrows in red, and _gambelii_ White-crowned Sparrows in blue.

### Genotype-by-individual plot of highly differentiated SNPs

#= Our overall dataset has close to 46,000 SNPs, far too many to visualize meaningfully in a single plot.
We can come up with a meaningful subset by examining variation on a single scaffold (i.e., chromosome in this case)
and only showing the SNPs that strongly differ in allele frequency between chosen groups.

We'll calculate allele frequencies and sample sizes for each group, 
and then genetic differentiation (known as F<sub>ST</sub>) between the groups:
 =#

freqs, sampleSizes = GenomicDiversity.getFreqsAndSampleSizes(sparrowGenoDataImputed, groupsToPlot)


Fst, FstNumerator, FstDenominator, pairwiseNamesFst = GenomicDiversity.getFst(freqs,
                                            sampleSizes, groupsToPlot)


# Now we will choose a scaffold and specify some other parameters for the plotting algorithm:

chr = "CM018231.2" # the name of a scaffold in the reference genome
regionInfo = chooseChrRegion(sparrowGenoData, chr) # this gets the maximum position for the chromosome
group1 = "GCSP"   # the alleles most common in this  group will be assigned the same color in the graph
groupsToCompare = "GCSP_PSWS" # The groups to compare for the Fst filter below
Fst_cutoff =  0.8
missingFractionAllowed = 0.2

# Finally, we can actually make the plot:

plotInfo = plotGenotypeByIndividualWithFst(sparrowGenoDataImputed, groupsToCompare, 
    Fst_cutoff, missingFractionAllowed, 
    regionInfo, Fst, pairwiseNamesFst, 
    freqs, groupsToPlot, groupColors)

# The above way of calling the function has the `GenoData` object as the first argument. 
# Another way is to pass the the three components of the `GenoData` separately: 

plotInfo = plotGenotypeByIndividualWithFst(groupsToCompare, 
    Fst_cutoff, missingFractionAllowed, regionInfo, 
    sparrowGenoDataImputed.positions, Fst, pairwiseNamesFst, 
    sparrowGenoDataImputed.genotypes, sparrowGenoDataImputed.indInfo, freqs, 
    groupsToPlot, groupColors)
plotInfo[1] # this outputs the plot


