# This file serves as an example script that I am using during development of GenomicDiversity.jl V0.3.0

# load required packages
# using GenomicDiversity, MultivariateStats, CairoMakie, DataFrames, CSV, DelimitedFiles, Impute

using MultivariateStats, CairoMakie, DataFrames, CSV, DelimitedFiles, Impute

import Pkg; Pkg.add(path="https://github.com/darreni/GenomicDiversity.jl")
using GenomicDiversity

# set data file names
pathName = "demoData/SparrowDemo_data_McCallumetal2024/"
metadataFile = pathName * "SparrowDemo.Fst_groups.txt"
genotypesFile = pathName * "SparrowDemo_genotypes.012"
individualsFile = genotypesFile * ".indv"
positionsFile = genotypesFile * ".pos"


#= """
    GenoData(indInfo::DataFrame, genotypes::Matrix, positions::DataFrame)

Constructs a `GenoData` struct, which stores metadata of individuals, a genotype matrix of individuals and SNPS, and genomic position info.

​# Arguments / keywords
- `indInfo`: A DataFrame containing metadata.
- `genotypes`: A Matrix containing genotypes for individuals (rows) and loci (columns).
- `positions`: A 2-column DataFrame containing genomic location info of SNPs.

# Notes
You can access each component of this data structure by appending ".indInfo" or ".genotypes" or "positions" to the name of the GenoData object.
"""
mutable struct GenoData
    indInfo::DataFrame
    genotypes::Matrix
    positions::DataFrame
    function GenoData(indInfo::DataFrame, genotypes::Matrix, positions::DataFrame)
        if size(indInfo, 1) != size(genotypes, 1)
            error("number of individuals (rows) differs between metadata dataframe and genotype matrix")
        end
        if size(genotypes, 2) != size(positions, 1)
            error("number of positions differs between genotype matrix and positions dataframe")
        end
        new(indInfo, genotypes, positions)
    end
end =#

# Now we load the data:

#= """
    loadGenoData(metadataFile, individualsFile, genotypesFile, positionsFile)

Imports data from four files and builds a `GenoData` object. 

​# Arguments
- `metadataFile`: The path and name of the metadata file.
- `individualsFile`: The path and name of the file providing names of the individuals (corresponding to rows of the genotype matrix).
- `genotypesFile`: The path and name of the file containing the matrix of genotypes.
- `positionsFile`: The path and name of the file containing the genomic position info of the loci corresponding to columns of the genotype matrix.

# Notes
Returns a `GenoData`` object, containing:
- `indInfo`: A DataFrame containing metadata.
- `genotypes`: A Matrix containing genotypes for individuals (rows) and loci (columns).
- `positions`: A 2-column DataFrame containing genomic location info of SNPs.
"""
function loadGenoData(metadataFile, individualsFile, genotypesFile, positionsFile)
    # load the metadata
    indMetadata = DataFrame(CSV.File(metadataFile))  # uses CSV, DataFrames packages
    # load the list of individuals
    indNames = DataFrame(CSV.File(individualsFile; header=["ind"], types=[String]))
    # check that rows in each imported object are the same:
    if nrow(indMetadata) != size(indNames, 1)
        error("Number of rows in metadataFile different than number of individuals in individualsFile")
    end
    # put individual names and metadata together, to enable easy confirmation that they match
    indNames_with_indMetadata = hcat(indNames, indMetadata)
    if hasproperty(indNames_with_indMetadata, "ID")
        if indNames_with_indMetadata.ind != indNames_with_indMetadata.ID
            println("WARNING: Individual names differ between $metadataFile and $individualsFile")
        end
    end
    # load the genomic positions of loci
    positions = DataFrame(CSV.File(positionsFile; header=["chrom", "position"], types=[String, Int]))
    # load the genotype matrix
    geno = readdlm(genotypesFile, '\t', Int16, '\n')  # uses DelimitedFiles package
    genosOnly = geno[:, Not(1)] #remove first column, which was just a row index
    if nrow(positions) != size(genosOnly, 2)
        error("Number of loci differs between $genotypesFile and $positionsFile")
    end
    println(string("Have read in genotypic data at ", size(genosOnly, 2)," loci for ", size(genosOnly, 1), " individuals. \n"))
    return GenoData(indNames_with_indMetadata, genosOnly, positions)
end =#

sparrowGenoData = loadGenoData(metadataFile, individualsFile, genotypesFile, positionsFile)


#= 
# load the metadata
sparrowMetadata = DataFrame(CSV.File(metadataFile))

# load the list of individuals
ind = DataFrame(CSV.File(individuals_file_name; header=["ind"], types=[String]))

# check that rows in each imported object are the same:
indNum = size(ind, 1) # number of individuals
if nrow(sparrowMetadata) != indNum
    println("WARNING: number of rows in metadata file different than number of individuals in .indv file")
else println("Good news: number of rows in metadata file matches the number of individuals in .indv file")
end

# put individual names and metadata together, to enable easy confirmation that they match
ind_with_metadata = hcat(ind, sparrowMetadata)

# load the genomic positions of loci
pos_whole_genome = DataFrame(CSV.File(position_file_name; header=["chrom", "position"], types=[String, Int]))

# load the genotype matrix
geno = readdlm(genotype_file_name, '\t', Int16, '\n')
loci_count = size(geno, 2) - 1   # because the first column is not a SNP (just a count from zero)
print(string("Read in genotypic data at ", loci_count," loci for ", indNum, " individuals. \n"))
genosOnly = geno[:, Not(1)] #remove first column, which was just a row index

sparrowGenoData = GenoData(ind_with_metadata, genosOnly, pos_whole_genome)
 =#

groupsToPlot = ["GCSP","PSWS","GWCS"]

groupColorKey = Dict("GCSP" => "gold",
                    "PSWS" => "red",
                    "GWCS" => "blue")

groupColors = [groupColorKey[i] for i in groupsToPlot]



### PCA2

#= If we want to use Principal Components Analysis (PCA) to visualize a good 2-dimensional representation
of genomic relationships among individuals, we first need to impute missing genotypes.
We will modify the type of matrix so that missing values are represented by `missing` rather than `-1`, 
and then we'll impute: =#


#= """
    imputeMissingGenotypes(genotypeMatrix::Matrix; method = "SVD")

Imputes missing values in the genotype matrix.

​# Arguments
- `genotypeMatrix`: Matrix of genotypes, where missing data are indicated with either `-1` or `missing`.
- `method`: The method used for imputing. Presently "SVD" is the only option (and default).  

# Notes
Returns the genotype matrix with filled-in missing values, of type Float32. 
"""
function imputeMissingGenoData(genotypeMatrix::Matrix; method = "SVD")
    T = eltype(genotypeMatrix)  # get the types of the cells, to use in next line
    genotypeMatrixWithMissing = Matrix{Union{Missing, T}}(genotypeMatrix)  # allow missing type
    replace!(genotypeMatrixWithMissing, -1 => missing)
    # convert to Float32 numbers
    genotypeMatrixWithMissing = Matrix{Union{Missing, Float32}}(genotypeMatrixWithMissing) 
    # impute the missing genotypes
    println("Imputing now. Thanks for your patience, as this may take some time.")
    return Impute.svd(genotypeMatrixWithMissing)  # uses Impute package
end

"""
    imputeMissingGenotypes(myGenoData::GenoData; method = "SVD")

When applied to a GenoData object, imputes the data in the genotypes object, and returns a new GenoData object.
"""
function imputeMissingGenoData(myGenoData::GenoData; method = "SVD")
    genotypesImputed = imputeMissingGenoData(myGenoData.genotypes; method = "SVD")
    return GenoData(myGenoData.indInfo, genotypesImputed, myGenoData.positions)
end =#

sparrowGenoDataImputed = imputeMissingGenoData(sparrowGenoData)

# Now we are ready to call the `GenomicDiversity.plotPCA()` function to do the PCA and make the plot:

#= PCAmodel = plotPCA(sparrowGenoDataImputed.genotypes, sparrowGenoDataImputed.indInfo,
        groupsToPlot, groupColors;
        sampleSet="Zonotrichia sparrows", regionText="whole_genome",
        flip1=false, flip2=false, showPlot=true)
# PCAmodel.PCAfig  # shows the figure =#


#= """
    plotPCA(genotypes, ind_with_metadata, groups_to_plot_PCA, group_colors_PCA; sampleSet = "", regionText="", flip1 = false, flip2 = false)

Do a principal components analysis based on genotypes of individuals (colored by group), and display PC1 vs. PC2.

​# Arguments
- `genotypes`: Genotype matrix (individuals in rows, loci in columns, with no missing data--can impute prior to this, to ensure no missing).
- `indMetadata`: Matrix containing metadata for individuals (make sure there is an `Fst_group`` column).
- `groups_to_plot_PCA`: Vector of names of groups to include.
- `group_colors_PCA`: Vector listing plotting color for each each group.
Optional arguments:
- `sampleSet`: Name of the sample set to appear in the plot title.
- `regionText`: Name of the genomic region to appear in the plot title.
- `flip1`: Set to `true` if wanting to flip PC1 (i.e., multiply by -1).
- `flip2`: Same but for PC2.
- `showPlot`: Set to `false` if not wanting the function to draw display.
- `autolimitaspect_setting`: For unconstrained axes, set to `nothing` (default is `1` for 1:1 axes).
- `lineOpacity`: Opacity of the line color of symbols.
- `fillOpacity`: Opacity of the fill color of symbols.
- `symbolSize`: Size of symbols.
- `plotTitle`: Specify if a specific title is wanted; otherwise will fill in a title.
- `showTitle`: Set to `false` for no title. 
- `xLabelText`: Can choose an x axis label.
- `yLabelText`: Can choose a y axis label.
- `labelSize`: Size of x and y axis labels.

# Notes

Returns a tuple containing:
1) `model`: the PCA model.
2) `metadata`: the metadata of individuals included in the PCA.
2) `values`: PC values of individuals (before any flipping of axes).
3) `PC1`: PC1 values of individuals(after any flipping).
4) `PC2`: PC2 values of individuals (after any flipping).
5) `PCAfig`: The PCA figure.
"""
function TESTplotPCA(genotypes, indMetadata, groups_to_plot_PCA, group_colors_PCA; 
                    sampleSet = "", regionText="",
                    flip1 = false, flip2 = false,
                    showPlot = true, autolimitaspect_setting = 1,
                    lineOpacity = 0.8, fillOpacity = 0.2,
                    symbolSize = 10, plotTitle = nothing, showTitle = true,
                    xLabelText = "Genomic PC1", yLabelText = "Genomic PC2",
                    labelSize = 24)
    
    selection = map(in(groups_to_plot_PCA), indMetadata.Fst_group)
    matrixForPCA = Matrix{Float32}(transpose(genotypes[selection, :]))
    indMetadata_groupSelected = indMetadata[selection, :]

    if showTitle && isnothing(plotTitle)
        if sampleSet == "" && regionText == ""
            plotTitle = "PCA"
        elseif sampleSet == ""
            plotTitle = string("PCA: ", regionText)
        elseif regionText == ""
            plotTitle = string("PCA of ", sampleSet)
        else 
            plotTitle = string("PCA of ", sampleSet, ": ", regionText)
        end
    elseif !showTitle  # if showTitle is false, then set title to ""
        plotTitle = ""
    end

    # this now works well--appears to NOT be scaling the variables to variance 1, which is good
    PCA_indGenos = fit(PCA, matrixForPCA; method = :svd, maxoutdim=3); # good to suppress output of this--otherwise many lines
    PCA_values = predict(PCA_indGenos, matrixForPCA)

    # eigvecs(PCA_indGenos)
    # loadings(PCA_indGenos)
    # eigvals(PCA_indGenos)
    # projection(PCA_indGenos)
    # tvar(PCA_indGenos)

    # if chosen to, flip axes:
    if flip1
        PC1 = -PCA_values[1,:]
    else
        PC1 = PCA_values[1,:]
    end
    if flip2
        PC2 = -PCA_values[2,:]
    else
        PC2 = PCA_values[2,:]
    end

    try
        global f = CairoMakie.Figure()
        global ax = Axis(f[1, 1],
            title = plotTitle,
            xlabel = xLabelText, xlabelsize = labelSize,
            ylabel = yLabelText, ylabelsize = labelSize,
            autolimitaspect = autolimitaspect_setting)
        hidedecorations!(ax, label = false, ticklabels = false, ticks = false) # hide background lattice
        for i in eachindex(groups_to_plot_PCA) 
            selection = indMetadata_groupSelected.Fst_group .== groups_to_plot_PCA[i]
            CairoMakie.scatter!(ax, PC1[selection], PC2[selection], marker = :diamond, color = (group_colors_PCA[i], fillOpacity), markersize = symbolSize, strokewidth=0.5, strokecolor = ("black", lineOpacity))
        end
        showPlot && display(f)
    catch # Makie sometimes has an error due to the "autolimitaspect = 1" above, so adding this "try-catch" structure to keep program going:
        println("Did not success in drawing PCA with 1:1 axes for ", regionText,
                ", so drawing with non-proportional axes.")
        global f = CairoMakie.Figure()
        global ax = Axis(f[1, 1],
            title = plotTitle,
            xlabel = xLabelText, xlabelsize = labelSize,
            ylabel = yLabelText, ylabelsize = labelSize)
        hidedecorations!(ax, label = false, ticklabels = false, ticks = false) # hide background lattice
        for i in eachindex(groups_to_plot_PCA) 
            selection = indMetadata_groupSelected.Fst_group .== groups_to_plot_PCA[i]
            CairoMakie.scatter!(ax, PC1[selection], PC2[selection], marker = :diamond, color = (group_colors_PCA[i], fillOpacity), markersize = symbolSize, strokewidth = 0.5, strokecolor = ("black", lineOpacity))
        end
        showPlot && display(f) 
    end 

    return (model = PCA_indGenos, metadata = indMetadata_groupSelected, values = PCA_values, PC1 = PC1, PC2 = PC2, PCAfig = f)
end

"""
    plotPCA(myGenoData::GenoData, groups_to_plot_PCA, group_colors_PCA; sampleSet = "", regionText="", flip1 = false, flip2 = false)

Apply the plotPCA function to a `GenoData`` object.
"""
function TESTplotPCA(myGenoData::GenoData, groups_to_plot_PCA, group_colors_PCA; 
                    sampleSet = "", regionText="",
                    flip1 = false, flip2 = false,
                    showPlot = true, autolimitaspect_setting = 1,
                    lineOpacity = 0.8, fillOpacity = 0.2,
                    symbolSize = 10, plotTitle = nothing, showTitle = true,
                    xLabelText = "Genomic PC1", yLabelText = "Genomic PC2",
                    labelSize = 24)
    return TESTplotPCA(myGenoData.genotypes, myGenoData.indInfo, groups_to_plot_PCA, group_colors_PCA; 
                        sampleSet = sampleSet, regionText = regionText,
                        flip1 = flip1, flip2 = flip2,
                        showPlot = showPlot, autolimitaspect_setting = autolimitaspect_setting,
                        lineOpacity = lineOpacity, fillOpacity = fillOpacity,
                        symbolSize = symbolSize, plotTitle = plotTitle, showTitle = showTitle,
                        xLabelText = xLabelText, yLabelText = yLabelText,
                        labelSize = labelSize)
end =#

#= PCAmodel = plotPCA(sparrowGenoDataImputed.genotypes, sparrowGenoDataImputed.indInfo,
        groupsToPlot, groupColors;
        sampleSet="Zonotrichia sparrows", regionText="whole_genome",
        flip1=false, flip2=false, showPlot=true)
PCAmodel.PCAfig  =#

PCAmodel = plotPCA(sparrowGenoDataImputed,
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

#= """
    getFreqsAndSampleSizes(genoData, indGroup, groupsToCalc)

Calculate allele frequencies and sample sizes for each group and SNP.

​# Arguments
- `genoData`: The genotype matrix, where rows are individuals and columns are loci, with genotype codes 0,1,2 meaning homozygous reference, heterozygote, homozygous alternate, and missing genotypes can be either -1 or `missing`.
- `indGroup`: A vector providing the group name each individual belongs to.
- `groupsToCalc`: A list of group names to include in calculations.

# Notes
Returns a tuple containing 1) a matrix of frequencies, and 2) a matrix of samples sizes (in both, rows are groups and columns are loci). 
"""
function TESTgetFreqsAndSampleSizes(genoData, indGroup, groupsToCalc)
    groupCount = length(groupsToCalc)
    freqs = Array{Float32,2}(undef, groupCount, size(genoData, 2))
    sampleSizes = Array{Int16,2}(undef, groupCount, size(genoData, 2))
    for i in 1:groupCount
        selection = (indGroup .== groupsToCalc[i]) # gets the correct rows for individuals in the group 
        geno0counts = sum(isequal.(genoData[selection, :], 0), dims=1) # count by column the number of 0 genotypes (homozygous ref)
        geno1counts = sum(isequal.(genoData[selection, :], 1), dims=1) # same for 1 genotypes (heterozygous)
        geno2counts = sum(isequal.(genoData[selection, :], 2), dims=1) # same for 2 genotypes (homozygous alternate) 
        sumGenoCounts = geno0counts .+ geno1counts .+ geno2counts
        sampleSizes[i, :] = sumGenoCounts
        freqs[i, :] = ((2 .* geno2counts) .+ geno1counts) ./ (2 * sumGenoCounts)
    end
    return freqs, sampleSizes
end 

"""
    getFreqsAndSampleSizes(myGenoData::GenoData, groupsToCalc)

Call getFreqsAndSampleSizes() using a `GenoData` object as input.
"""
function TESTgetFreqsAndSampleSizes(myGenoData::GenoData, groupsToCalc)
    return TESTgetFreqsAndSampleSizes(myGenoData.genotypes,
                                    myGenoData.indInfo.Fst_group, groupsToCalc)
end =#

# freqs, sampleSizes = getFreqsAndSampleSizes(sparrowGenoDataImputed.genotypes,
#                    sparrowGenoDataImputed.indInfo.Fst_group, groupsToPlot)

freqs, sampleSizes = getFreqsAndSampleSizes(sparrowGenoDataImputed, groupsToPlot)


Fst, FstNumerator, FstDenominator, pairwiseNamesFst = getFst(freqs,
                                            sampleSizes, groupsToPlot)


# Now we will choose a scaffold and specify some other parameters for the plotting algorithm:

#= """
    chooseChrRegion(pos, chr; positionMin=1, positionMax=NaN)

Designate a specific region of a chromosome for subsequent analysis.

​# Arguments
- `pos`: A matrix containing genomic location of each locus; must have a `chrom` column (name of scaffold) and a `position` column (numerical position on scaffold).
- `chr`: The name of the chosen scaffold.
- `positionMin`: Optional (default 1); the starting location of the chosen region.
- `positionMax`: Optional (default the end of the scaffold); the ending location of the chosen region.

# Notes
Returns a tuple containing:
- the scaffold name
- the minimum position
- the maximum position
- a string describing the region (e.g. `chr Z: 20000 to 1000000`)

This function does not actually filter the genotype matrix--that can be done by providing the output of this to functions such as `plotGenotypeByIndividual()`.
"""
function TESTchooseChrRegion(pos, chr; positionMin=1, positionMax=NaN)
    if positionMin == 1 && isnan(positionMax) # then include whole chromosome
        positionMax = last(pos.position[pos.chrom .== chr]) # get highest position value on chromosome
        regionText = string(chr,"_whole")
    elseif positionMin > 1 && isnan(positionMax)
        positionMax = last(pos.position[pos.chrom .== chr]) # get highest position value on chromosome
        regionText = string(chr,"_from",positionMin,"to",positionMax)
    else
        regionText = string("chr ", chr, " ",positionMin," to ",positionMax)
    end
    return chr, positionMin, positionMax, regionText
end

"""
    chooseChrRegion(myGenoData::GenoData, chr; positionMin=1, positionMax=NaN)

Call chooseChrRegion using a GenoData object.
"""
function TESTchooseChrRegion(myGenoData::GenoData, chr; positionMin=1, positionMax=NaN)
    return chooseChrRegion(myGenoData.positions, chr; positionMin=1, positionMax=NaN)
end =#

chr = "CM018231.2" # the name of a scaffold in the reference genome
regionInfo = chooseChrRegion(sparrowGenoData, chr) # this gets the maximum position for the chromosome
group1 = "GCSP"   # the alleles most common in this  group will be assigned the same color in the graph
groupsToCompare = "GCSP_PSWS" # The groups to compare for the Fst filter below
Fst_cutoff =  0.8
missingFractionAllowed = 0.2

# Finally, we can actually make the plot:

#= """
    plotGenotypeByIndividualWithFst(groupsToCompare, Fst_cutoff, missingFractionAllowed,
                            regionInfo, pos, Fst, pairwiseNamesFst,
                            genoData, indMetadata, freqs, plotGroups, plotGroupColors;
                            colorAllelesByGroup = true, group1 = plotGroups[1],
                            indFontSize=10, figureSize=(1200, 1200),
                            show_SNP_density = true, densityPlotColor = :steelblue1,
                            plotTitle = nothing, titleFontSize = 20,
                            highlightRegionStarts = [],
                            highlightRegionEnds = [],
                            highlightRegionColor = "red")

Construct a genotype-by-individual plot, including only loci that pass thresholds for Fst and missing data. 

Under the default setting, alleles are colored (dark purple vs. light purple) according to whichever allele is designated as `group1`. 

​# Arguments
- `groupsToCompare`: The two groups (in format `name1_name2`) for filtering loci based on Fst. This can be a single string, or a vector of strings, each with a population comparison to include if Fst is above Fst_cutoff.
- `Fst_cutoff`: The minimum Fst for a locus to be included in the plot.
- `missingFractionAllowed`: The maximum missing genotype fraction for a locus to be included.
- `regionInfo`: Information regarding the scaffold and region of focus; a tuple provided by `chooseChrRegion()`.
- `pos`: Matrix providing genomic location of each locus; there must be a `chrom` column and a `position` column.
- `Fst`: Matrix of Fst values; can be produced by `getFst()`.
- `pairwiseNamesFst`: Vector of pairwise names corresponding to rows in the `Fst` matrix.
- `genoData`: Matrix containing genotype data (individuals in rows, loci in columns).
- `indMetadata`: Matrix of metadata for individuals; must contain `Fst_group` and `plot_order` columns.
- `freqs`: Matrix of alternate allele frequencies for each group (row) and locus (column).
- `plotGroups`: Vector of group names to include in plot.
- `plotGroupColors`: Vector of plotting colors corresponding to the groups.
- `colorAllelesByGroup`: Optional; set to `false` to color alleles according to reference and alternate.
- `group1`: Optional (default is `plotGroups[1]`); when `colorAllelesByGroup` is `true`, this is the group that determine which allele is dark purple. 
- `indFontSize`: Optional; the font size of the individual ID labels.
- `figureSize`: Optional; the size of the figure; default is `(1200, 1200)`.  
- `show_SNP_density`: Optional; default is `true` to show a density plot. 
- `densityPlotColor`: Optional; default is `:steelblue1`
- `plotTitle`: Optional; default will make a title. For no title, set to `""`.
- `titleFontSize`; Optional; default is `20`.
- `highlightRegionStarts`: Optional; the left locations of regions to highlight on the scaffold.
- `highlightRegionEnds`: Optional; the right locations of regions to highlight on the scaffold. 
- `highlightRegionColor`: Optional; default red; the color of the highlight bar along the scaffold.

# Notes
Returns a tuple containing:
- the figure
- the plotted genotypes
- the numerical positions (in the chosen scaffold) of the plotted loci
- the sorted metadata matrix for the plotted individuals
"""
function TESTplotGenotypeByIndividualWithFst(groupsToCompare, Fst_cutoff, missingFractionAllowed,
                            regionInfo, pos, Fst, pairwiseNamesFst,
                            genoData, indMetadata, freqs, plotGroups, plotGroupColors;
                            colorAllelesByGroup = true, group1 = plotGroups[1],
                            indFontSize=10, figureSize=(1200, 1200),
                            show_SNP_density = true, densityPlotColor = :steelblue1,
                            plotTitle = nothing, titleFontSize = 20,
                            highlightRegionStarts = [],
                            highlightRegionEnds = [],
                            highlightRegionColor = "red")

    chr, positionMin, positionMax, regionText = regionInfo
    if isnothing(plotTitle)
        plotTitle = string(regionText, ": genotypes Fst>", Fst_cutoff, " loci between ", groupsToCompare)
    end

    # write ≥ using \ge tab, and ≤ using \le tab.
    # note the operator ".&" (bitwise &) is much more efficient than ".&&"
    # Select based on the genome region:
    selection_position = (pos.chrom .== chr) .&
                (pos.position .≥ positionMin) .&
                (pos.position .≤ positionMax)

    # Filter SNPs according to Fst of the groupsToCompare
    if typeof(groupsToCompare) == String  # If one scaffold name was provided
        rowChoice = findfirst(pairwiseNamesFst .== groupsToCompare)
        selection_Fst = (.!isnan.(Fst[rowChoice, :])) .& (Fst[rowChoice, :] .≥ Fst_cutoff)
    elseif typeof(groupsToCompare) == Vector{String}
        rowChoice = findfirst(pairwiseNamesFst .== groupsToCompare[1])
        selection_Fst = (.!isnan.(Fst[rowChoice, :])) .& (Fst[rowChoice, :] .≥ Fst_cutoff)
        if size(groupsToCompare, 1) ≥ 2
            for i in eachindex(groupsToCompare)[2:end]
                rowChoice = findfirst(pairwiseNamesFst .== groupsToCompare[i])
                selection_Fst = selection_Fst .| (.!isnan.(Fst[rowChoice, :])) .& (Fst[rowChoice, :] .≥ Fst_cutoff)
            end
        end
    end
    selection = selection_position .& selection_Fst
    SNP_genotypes = genoData[:, selection]
    SNP_freqs = freqs[:, selection]
    SNP_positions = pos.position[selection]
    SNP_genotypes_subset = SNP_genotypes[indMetadata.Fst_group .∈ Ref(plotGroups), :]
    # convert missing to -99
    SNP_genotypes_subset[ismissing.(SNP_genotypes_subset)] .= -99

    indMetadata_subset = indMetadata[indMetadata.Fst_group .∈ Ref(plotGroups), :]
    # The section below does the default of coloring the allele most common in group1 as the dark purple.
    # Otherwise, set colorAllelesByGroup=false to have the alleles colored according to ref vs. alternate 
    if colorAllelesByGroup
        altAlleleHiInGroup1 = SNP_freqs[findfirst(plotGroups .== group1), :] .> 0.5
        SNP_genotypes_subset[:, altAlleleHiInGroup1] = 2 .- SNP_genotypes_subset[:, altAlleleHiInGroup1]
    end
    # Choose sorting order by plot_order column in input metadata file
    #sorted.SNP.genotypes.subset = SNP.genotypes.subset[order(SNP.genotypes.subset$group, SNP.genotypes.subset$ID),]
    sorted_SNP_genotypes_subset = SNP_genotypes_subset[sortperm(indMetadata_subset.plot_order, rev=false), :]
    numInds = size(sorted_SNP_genotypes_subset, 1) # fast way to get number of rows
    sorted_indMetadata_subset = indMetadata_subset[sortperm(indMetadata_subset.plot_order, rev=false), :]

    # filter out the SNPs that have too much missing data:
    numberMissing = sum(isequal.(sorted_SNP_genotypes_subset, -1) .| 
                        isequal.(sorted_SNP_genotypes_subset, 3) .|
                        ismissing.(sorted_SNP_genotypes_subset), dims=1)
    fractionMissing = numberMissing / numInds
    selection = vec(fractionMissing .≤ missingFractionAllowed)
    sorted_SNP_genotypes_subset2 = sorted_SNP_genotypes_subset[:, selection]
    SNP_positions_subset2 = SNP_positions[selection]
    num_SNPs_to_plot = length(SNP_positions_subset2)

    # Set up the plot window:
    f = CairoMakie.Figure(size=figureSize)

    # Set up the main Axis: 
    ax = Axis(f[1, 1],
        title = plotTitle,
        # xlabel = "location",
        # ylabel = "individual",
        titlesize = titleFontSize,
        limits=(0.5 - 0.09 * (num_SNPs_to_plot), 0.5 + 1.09 * (num_SNPs_to_plot),
            0.5 - 0.3 * numInds, 0.5 + numInds)
    )
    hidedecorations!(ax) # hide background lattice and axis labels
    hidespines!(ax) # hide box around plot

    genotypeColors = ["#3f007d", "#807dba", "#dadaeb", "grey50"]  # purple shades from colorbrewer

    # plot evenly spaced by SNP order along chromosome:
    # make top part of fig (genotypes for individuals)
    labelCushion = num_SNPs_to_plot / 100
    label_x_left = 0.5 - labelCushion
    label_x_right = 0.5 + num_SNPs_to_plot + labelCushion
    colorBoxCushion = 0.07 * num_SNPs_to_plot
    groupColorBox_x_left = 0.5 - colorBoxCushion
    groupColorBox_x_right = 0.5 + num_SNPs_to_plot + colorBoxCushion
    boxWidth = 0.005 * num_SNPs_to_plot * 2
    groupColorBox_x_left = [-boxWidth, -boxWidth, boxWidth, boxWidth, -boxWidth] .+ groupColorBox_x_left
    groupColorBox_x_right = [-boxWidth, -boxWidth, boxWidth, boxWidth, -boxWidth] .+ groupColorBox_x_right
    groupColorBox_y = [0.4, -0.4, -0.4, 0.4, 0.4]

    for i in 1:numInds
        y = numInds + 1 - i  # y is location for plotting; this reverses order of plot top-bottom
        labelText = last(split(sorted_indMetadata_subset.ID[i], "_"))  # this gets the last part of the sample ID (usually the main ID part)
        # put sample label on left side:
        CairoMakie.text!(label_x_left, y; text=labelText, align=(:right, :center), fontsize=indFontSize)
        # put sample label on left side:
        CairoMakie.text!(label_x_right, y; text=labelText, align=(:left, :center), fontsize=indFontSize)
        boxColor = plotGroupColors[findfirst(plotGroups .== sorted_indMetadata_subset.Fst_group[i])]
        CairoMakie.poly!(Point2f.(groupColorBox_x_left, (y .+ groupColorBox_y)), color=boxColor)
        CairoMakie.poly!(Point2f.(groupColorBox_x_right, (y .+ groupColorBox_y)), color=boxColor)
    end

    # generate my own plotting symbol (a rectangle)
    box_x = [-0.5, -0.5, 0.5, 0.5, -0.5]
    box_y = [0.4, -0.4, -0.4, 0.4, 0.4]
    # generate triangles for plotting heterozygotes
    triangle1_x = [-0.5, -0.5, 0.5, -0.5]
    triangle1_y = [0.4, -0.4, 0.4, 0.4]
    triangle2_x = [-0.5, 0.5, 0.5, -0.5]
    triangle2_y = [-0.4, -0.4, 0.4, -0.4]
    # cycle through individuals, graphing each type of genotype:
    for i in 1:numInds
        y = numInds + 1 - i  # y is location for plotting; this reverses order of plot top-bottom
        CairoMakie.lines!([0.5, num_SNPs_to_plot + 0.5], [y, y], color="grey40")
        genotypes = sorted_SNP_genotypes_subset2[i, :]
        hom_ref_locs = findall(isequal.(genotypes, 0))
        if length(hom_ref_locs) > 0
            for j in eachindex(hom_ref_locs)
                CairoMakie.poly!(Point2f.((hom_ref_locs[j] .+ box_x), (y .+ box_y)), color=genotypeColors[1])
                #polygon(hom.ref.locs[j]+symbol.x, y+symbol.y, border=NA, col=genotype.colors[1])
            end
        end
        het_locs = findall(isequal.(genotypes, 1))
        if length(het_locs) > 0
            for j in eachindex(het_locs)
                CairoMakie.poly!(Point2f.((het_locs[j] .+ triangle1_x), (y .+ triangle1_y)), color=genotypeColors[1])
                CairoMakie.poly!(Point2f.((het_locs[j] .+ triangle2_x), (y .+ triangle2_y)), color=genotypeColors[3])
            end
        end
        hom_alt_locs = findall(isequal.(genotypes, 2))
        if length(hom_alt_locs) > 0
            for j in eachindex(hom_alt_locs)
                CairoMakie.poly!(Point2f.((hom_alt_locs[j] .+ box_x), (y .+ box_y)), color=genotypeColors[3])
            end
        end
    end

    # make lower part of figure (indicating position along chromosome)
    chrLine_y = 0.5 - 0.2numInds
    topHatchLine_y1 = 0.6
    topHatchLine_y2 = 0.5 - 0.02numInds
    lowHatchLine_y1 = 0.5 - 0.18numInds
    lowHatchLine_y2 = 0.5 - 0.2numInds
    # draw chromosome line and text labels
    CairoMakie.lines!([0.5, num_SNPs_to_plot + 0.5], [chrLine_y, chrLine_y], color="black", linewidth=5)
    CairoMakie.text!(0.5 + (num_SNPs_to_plot / 2), chrLine_y - 0.04numInds; text=string("Location along chromosome ", chr), align=(:center, :center), fontsize=30)
    CairoMakie.text!(0.5, chrLine_y - 0.025numInds; text=string(positionMin), align=(:center, :center), fontsize=20)
    CairoMakie.text!(0.5 + num_SNPs_to_plot, chrLine_y - 0.025numInds; text=string(positionMax), align=(:center, :center), fontsize=20)
    # draw lines from SNPs to chromosome line
    chrPlotRatio = num_SNPs_to_plot / (positionMax - positionMin)
    for i in eachindex(SNP_positions_subset2)
        CairoMakie.lines!([i, i], [topHatchLine_y1, topHatchLine_y2], color="grey20", linewidth=1)
        CairoMakie.lines!([i, 0.5 + chrPlotRatio * (SNP_positions_subset2[i] - positionMin)], [topHatchLine_y2, lowHatchLine_y1], color="grey20", linewidth=1)
        CairoMakie.lines!([0.5 + chrPlotRatio * (SNP_positions_subset2[i] - positionMin), 0.5 + chrPlotRatio * (SNP_positions_subset2[i] - positionMin)], [lowHatchLine_y1, lowHatchLine_y2], color="grey20", linewidth=1)
    end
    # draw highlight regions on chromosome
    highlightRegionLine_y = 0.5 - 0.21numInds
    yValues = [highlightRegionLine_y, highlightRegionLine_y]
    if length(highlightRegionStarts) > 0
        for i in eachindex(highlightRegionStarts)
            xValues = [0.5 + chrPlotRatio * (highlightRegionStarts[i] - positionMin), 0.5 + chrPlotRatio * (highlightRegionEnds[i] - positionMin)]
            CairoMakie.lines!(xValues, yValues, color=highlightRegionColor, linewidth=18)
        end
    end

    if show_SNP_density
        # show SNP density plot, using an inset axis
        selection = (pos.chrom .== chr) .&
                    (pos.position .≥ positionMin) .&
                    (pos.position .≤ positionMax)
        pos_region = pos[selection, :] 
        CairoMakie.text!(0.5 + (num_SNPs_to_plot / 2), chrLine_y - 0.075numInds; text=string("Showing ", num_SNPs_to_plot, " SNPs out of ", size(pos_region, 1), " (distribution in lower plot)"), align=(:center, :center), fontsize=20)
        density_plot_height = 0.05
        inset_ax = Axis(f[1, 1],
                    width = Relative(1 / 1.18), # calculated these from the above settings
                    height = Relative(density_plot_height), 
                    halign = :center,
                    valign = 0.1057 / 1.3,  # I had to sort of titrate this--the calculations (based on the above) were slightly off
                    backgroundcolor=:white,
                    limits=(positionMin, positionMax, 0, nothing))
        hidedecorations!(inset_ax) # hide background lattice and axis labels
        hidespines!(inset_ax)
        density!(pos_region.position, color = (densityPlotColor, 0.5))   #color = (:lightsteelblue1, 0.5)     
    end

    display(f)

    return f, sorted_SNP_genotypes_subset2, SNP_positions_subset2, sorted_indMetadata_subset
end

"""
    plotGenotypeByIndividualWithFst(myGenoData::GenoData, groupsToCompare, 
                            Fst_cutoff, missingFractionAllowed,
                            regionInfo, Fst, pairwiseNamesFst,
                            freqs, plotGroups, plotGroupColors;
                            colorAllelesByGroup = true, group1 = plotGroups[1],
                            indFontSize=10, figureSize=(1200, 1200),
                            show_SNP_density = true, densityPlotColor = :steelblue1,
                            plotTitle = nothing, titleFontSize = 20,
                            highlightRegionStarts = [],
                            highlightRegionEnds = [],
                            highlightRegionColor = "red")

Call the plotGenotypeByIndividualWithFst() function using a GenoData object.
"""
function TESTplotGenotypeByIndividualWithFst(myGenoData::GenoData, groupsToCompare, 
                            Fst_cutoff, missingFractionAllowed,
                            regionInfo, Fst, pairwiseNamesFst,
                            freqs, plotGroups, plotGroupColors;
                            colorAllelesByGroup = true, group1 = plotGroups[1],
                            indFontSize = 10, figureSize = (1200, 1200),
                            show_SNP_density = true, densityPlotColor = :steelblue1,
                            plotTitle = nothing, titleFontSize = 20,
                            highlightRegionStarts = [],
                            highlightRegionEnds = [],
                            highlightRegionColor = "red")
    return plotGenotypeByIndividualWithFst(groupsToCompare, Fst_cutoff, missingFractionAllowed,
                            regionInfo, myGenoData.positions, Fst, pairwiseNamesFst,
                            myGenoData.genotypes, myGenoData.indInfo, freqs, plotGroups, plotGroupColors;
                            colorAllelesByGroup = colorAllelesByGroup, group1 = group1,
                            indFontSize = indFontSize, figureSize = figureSize,
                            show_SNP_density = show_SNP_density, densityPlotColor = densityPlotColor,
                            plotTitle = plotTitle, titleFontSize = titleFontSize,
                            highlightRegionStarts = highlightRegionStarts,
                            highlightRegionEnds = highlightRegionEnds,
                            highlightRegionColor = highlightRegionColor)
end =#

plotInfo = plotGenotypeByIndividualWithFst(groupsToCompare, 
    Fst_cutoff, missingFractionAllowed, regionInfo, 
    sparrowGenoDataImputed.positions, Fst, pairwiseNamesFst, 
    sparrowGenoDataImputed.genotypes, sparrowGenoDataImputed.indInfo, freqs, 
    groupsToPlot, groupColors)
plotInfo[1] # this outputs the plot

# GenomicDiversity.plotGenotypeByIndividualWithFst()

plotInfo = plotGenotypeByIndividualWithFst(sparrowGenoDataImputed, groupsToCompare, 
    Fst_cutoff, missingFractionAllowed, 
    regionInfo, Fst, pairwiseNamesFst, 
    freqs, groupsToPlot, groupColors)
