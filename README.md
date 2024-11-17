# GenomicDiversity

[![Build Status](https://github.com/darreni/GenomicDiversity.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/darreni/GenomicDiversity.jl/actions/workflows/CI.yml?query=branch%3Amain)

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://darreni.github.io/GenomicDiversity.jl/dev)

## Description

GenomicDiversity is a Julia package for the analysis of genomic variation within and among multiple populations. 
Using input data that is a matrix of biallelic genotypes, with individuals in rows and loci in columns, 
GenomicDiversity can calculate per-locus and windowed within-group diversity, 
between group genetic distance, and relative differentiation between populations.
It can produce plots that illustrate various aspects of genomic diversity.

## Citation

If you use GenomicDiversity.jl, please cite this paper which will be the first to formally mention this package:

Irwin, D., S. Bensch, C. Charlebois, G. David, A. Geraldes, S.K. Gupta, B. Harr, P. Holt, J.H. Irwin, V.V. Ivanitskii, I.M. Marova, Y. Niu, S. Seneviratne, A. Singh, Y. Wu, S. Zhang, T.D. Price. 2024. The distribution and dispersal of large haploblocks in a superspecies. Preprint (not yet peer-reviewed) on _Authorea_: [https://doi.org/10.22541/au.173114974.42895933/v1](https://doi.org/10.22541/au.173114974.42895933/v1)

## Disclaimer

The functions in this package were developed by Darren Irwin over a number of years, for use in various analyses. 
They are now mature enough that I think it is an appropriate time to make them available as a Julia package.
However, I am certain they could be improved in various ways. My intent is to revise and update this package over time, 
and I welcome feedback and suggestions regarding how to best do so.

## Example analysis using GenomicDiversity.jl

The example here will be based on a dataset from this paper:

McCallum, Q., K. Askelson, F.F. Fogarty, L. Natola, E. Nikelski, A. Huang, and D. Irwin. 2024. Pronounced differentiation on the Z chromosome and parts of the autosomes in crowned sparrows contrasts with mitochondrial paraphyly: implications for speciation. _Journal of Evolutionary Biology_, voae004, [https://doi.org/10.1093/jeb/voae004](https://doi.org/10.1093/jeb/voae004)

We will load in a large dataset of genotypes throughout the genomes of Golden-crowned Sparrows and two forms of White-crowned Sparrows: _pugetensis_ and _gambelii_. 
The data files are within the folder `demoData/SparrowDemo_data_McCallumetal2024/` of the GenomicDiversity.jl repository.

These input files were produced by an analysis pipeline that took genotyping-by-sequencing (GBS) reads, 
mapped them to a reference genome, identified single nucleotide polymorphisms (SNPs), called genotypes of individuals,
and outputted the information in three files: 
* a genotype file containing a matrix of genotypes (rows representing individuals and columns representing loci), where 0 represents homozygous reference, 1 represents heterozygous, 2 represents homozygous alternate, and -1 represents missing data
* a `.indv` file containing individual identifiers (corresponding to rows of the genotype matrix)
* a `.pos` file providing the locations of loci in the genome (with rows representing loci, column 1 indicating the scaffold, and column 2 indicating the position along that scaffold)

First, we'll do some setup (loading packages and data):

```julia
# load required packages
using MultivariateStats, CairoMakie, DataFrames, CSV, DelimitedFiles

# set data file names
genotype_file_name = "SparrowDemo_genotypes.012"
individuals_file_name = genotype_file_name * ".indv"
position_file_name = genotype_file_name * ".pos"
metadataFile = "SparrowDemo.Fst_groups.txt"
```

Now we load the data:

```julia
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
```

By examining the size of the `genosOnly` matrix, we can see that the data matrix includes 47 individuals and 45,986 loci. 

We choose the groups to include in our analysis, and the colors to represent them. 
These names refer to entries in the `Fst_group` column of the metadata file. 

```julia
groups_to_plot = ["GCSP","PSWS","GWCS"]
group_colors = ["gold","red","blue"]
```

### PCA

If we want to use Principal Components Analysis (PCA) to visualize a good 2-dimensional representation
of genomic relationships among individuals, we first need to impute missing genotypes.
We will modify the type of matrix so that missing values are represented by `missing` rather than `-1`, 
and then we'll impute:

```julia
# change missing data format
genosOnly_with_missing = Matrix{Union{Missing, Int16}}(genosOnly)
genosOnly_with_missing[genosOnly_with_missing .== -1] .= missing;
genosOnly_with_missing = Matrix{Union{Missing, Float32}}(genosOnly_with_missing) 

# impute the missing genotypes
imputed_genos = Impute.svd(genosOnly_with_missing)
```

Now we are ready to call the `GenomicDiversity.plotPCA()` function to do the PCA and make the plot:

```julia
PCAmodel = GenomicDiversity.plotPCA(imputed_genos, ind_with_metadata,
        groups_to_plot, group_colors;
        sampleSet="Zonotrichia sparrows", regionText="whole_genome",
        flip1=false, flip2=false, showPlot=true)
PCAmodel.PCAfig  # shows the figure
```

The resulting plot:

![](demoData/SparrowDemo_data_McCallumetal2024/SparrowDemo_PCA.jpg)

The above shows genomic variation illustrated in 2 dimensions, showing Golden-crowned sparrows in yellow, _pugetensis_ White-crowned Sparrows in red, and _gambelii_ White-crowned Sparrows in blue.

### Genotype-by-individual plot of highly differentiated SNPs

Our overall dataset has close to 46,000 SNPs, far too many to visualize meaningfully in a single plot.
We can come up with a meaningful subset by examining variation on a single scaffold (i.e., chromosome in this case)
and only showing the SNPs that strongly differ in allele frequency between chosen groups.

We'll calculate allele frequencies and sample sizes for each group, 
and then genetic differentiation (known as F<sub>ST</sub>) between the groups:

```julia
freqs, sampleSizes = GenomicDiversity.getFreqsAndSampleSizes(genosOnly_with_missing,
                            ind_with_metadata.Fst_group, groups_to_plot)

Fst, FstNumerator, FstDenominator, pairwiseNamesFst = GenomicDiversity.getFst(freqs,
                                            sampleSizes, groups_to_plot)
```

Now we will choose a scaffold and specify some other parameters for the plotting algorithm:

```julia
chr = "CM018231.2" # the name of a scaffold in the reference genome
regionInfo = GenomicDiversity.chooseChrRegion(pos_whole_genome, chr) # this gets the maximum position for the chromosome
group1 = "GCSP"   # the alleles most common in this  group will be assigned the same color in the graph
groupsToCompare = "GCSP_PSWS" # The groups to compare for the Fst filter below
Fst_cutoff =  0.8
missingFractionAllowed = 0.2
```

Finally, we can actually make the plot:

```julia
plotInfo = GenomicDiversity.plotGenotypeByIndividualWithFst(groupsToCompare, 
    Fst_cutoff, missingFractionAllowed, regionInfo, 
    pos_whole_genome, Fst, pairwiseNamesFst, 
    genosOnly_with_missing, ind_with_metadata, freqs, 
    groups_to_plot, group_colors)
plotInfo[1] # this outputs the plot
```

The resulting plot:

![](demoData/SparrowDemo_data_McCallumetal2024/SparrowDemo_GBI.jpg)

The above shows genotypic variation among _Zonotrichia_ sparrows along a single chromosome. Each row represents an individual, with leftmost and rightmost colors indicating  Golden-crowned sparrows in yellow, _pugetensis_ White-crowned Sparrows in red, and _gambelii_ White-crowned Sparrows in blue. Dark purple and light purple boxes represent homozygotes for alternate loci, and boxes with triangles represent heterozygotes. Missing data is indicated by a horizontal line. Only highly differentiated (F<sub>ST</sub> > 0.8) loci are shown.