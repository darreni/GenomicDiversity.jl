# This file serves as an example script that I am using during development of GenomicDiversity.jl V0.3.0

# load required packages
using GenomicDiversity, MultivariateStats, CairoMakie, DataFrames, CSV, DelimitedFiles

# set data file names
genotype_file_name = "SparrowDemo_genotypes.012"
individuals_file_name = genotype_file_name * ".indv"
position_file_name = genotype_file_name * ".pos"
metadataFile = "SparrowDemo.Fst_groups.txt"


# Now we load the data:

# load the metadata
metadata = DataFrame(CSV.File(metadataFile))

# load the list of individuals
ind = DataFrame(CSV.File(individuals_file_name; header=["ind"], types=[String]))

# check that rows in each imported object are the same:
indNum = size(ind, 1) # number of individuals
if nrow(metadata) != indNum
    println("WARNING: number of rows in metadata file different than number of individuals in .indv file")
else println("Good news: number of rows in metadata file matches the number of individuals in .indv file")
end

# put individual names and metadata together, to enable easy confirmation that they match
ind_with_metadata = hcat(ind, metadata)

# load the genomic positions of loci
pos_whole_genome = DataFrame(CSV.File(position_file_name; header=["chrom", "position"], types=[String, Int]))

# load the genotype matrix
geno = readdlm(genotype_file_name, '\t', Int16, '\n')
loci_count = size(geno, 2) - 1   # because the first column is not a SNP (just a count from zero)
print(string("Read in genotypic data at ", loci_count," loci for ", indNum, " individuals. \n"))
genosOnly = geno[:, Not(1)] #remove first column, which was just a row index

groups_to_plot = ["GCSP","PSWS","GWCS"]
group_colors = ["gold","red","blue"]


### PCA

#= If we want to use Principal Components Analysis (PCA) to visualize a good 2-dimensional representation
of genomic relationships among individuals, we first need to impute missing genotypes.
We will modify the type of matrix so that missing values are represented by `missing` rather than `-1`, 
and then we'll impute: =#

# change missing data format
genosOnly_with_missing = Matrix{Union{Missing, Int16}}(genosOnly)
genosOnly_with_missing[genosOnly_with_missing .== -1] .= missing;
genosOnly_with_missing = Matrix{Union{Missing, Float32}}(genosOnly_with_missing) 

# impute the missing genotypes
imputed_genos = Impute.svd(genosOnly_with_missing)

# Now we are ready to call the `GenomicDiversity.plotPCA()` function to do the PCA and make the plot:

PCAmodel = GenomicDiversity.plotPCA(imputed_genos, ind_with_metadata,
        groups_to_plot, group_colors;
        sampleSet="Zonotrichia sparrows", regionText="whole_genome",
        flip1=false, flip2=false, showPlot=true)
PCAmodel.PCAfig  # shows the figure

# The above shows genomic variation illustrated in 2 dimensions, showing Golden-crowned sparrows in yellow, _pugetensis_ White-crowned Sparrows in red, and _gambelii_ White-crowned Sparrows in blue.

### Genotype-by-individual plot of highly differentiated SNPs

#= Our overall dataset has close to 46,000 SNPs, far too many to visualize meaningfully in a single plot.
We can come up with a meaningful subset by examining variation on a single scaffold (i.e., chromosome in this case)
and only showing the SNPs that strongly differ in allele frequency between chosen groups.

We'll calculate allele frequencies and sample sizes for each group, 
and then genetic differentiation (known as F<sub>ST</sub>) between the groups:
 =#

freqs, sampleSizes = GenomicDiversity.getFreqsAndSampleSizes(genosOnly_with_missing,
                            ind_with_metadata.Fst_group, groups_to_plot)

Fst, FstNumerator, FstDenominator, pairwiseNamesFst = GenomicDiversity.getFst(freqs,
                                            sampleSizes, groups_to_plot)


# Now we will choose a scaffold and specify some other parameters for the plotting algorithm:

chr = "CM018231.2" # the name of a scaffold in the reference genome
regionInfo = GenomicDiversity.chooseChrRegion(pos_whole_genome, chr) # this gets the maximum position for the chromosome
group1 = "GCSP"   # the alleles most common in this  group will be assigned the same color in the graph
groupsToCompare = "GCSP_PSWS" # The groups to compare for the Fst filter below
Fst_cutoff =  0.8
missingFractionAllowed = 0.2

# Finally, we can actually make the plot:

plotInfo = GenomicDiversity.plotGenotypeByIndividualWithFst(groupsToCompare, 
    Fst_cutoff, missingFractionAllowed, regionInfo, 
    pos_whole_genome, Fst, pairwiseNamesFst, 
    genosOnly_with_missing, ind_with_metadata, freqs, 
    groups_to_plot, group_colors)
plotInfo[1] # this outputs the plot

