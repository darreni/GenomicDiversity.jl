#= 
GenomicDiversity.jl
This file started by Darren Irwin on 24 June 2023
as a set of functions for processing and displaying SNP data.
    
Converted to an actual package (private) with a Github repo on 11 November 2023.

Changed name from SNPlots to GenoPlots to GenomicDiversity on 15Nov2024
    =#

# Useful site for developing packages: https://julialang.org/contribute/developing_package/

# Notes for self, regarding the development version:
# After making changes to this package, run: 
# julia> ] # Go to the package mode
# (GenomicDiversity) pkg> activate .
# (GenomicDiversity) pkg> test # to run basic test
# Then, where being used, need to run:
# import Pkg; Pkg.add(path="/Users/darrenirwin/.julia/dev/GenomicDiversity.jl")
# using GenomicDiversity



module GenomicDiversity

# make functions available in calling namespace through "using GenomicDiversity" command, without needing to prefix module name: 
export greet_GenomicDiversity,
    getFreqsAndSampleSizes,
    getSitePi,
    getRegionPi,
    getPairwiseNames,
    getDxy,
    getRegionDxy,
    getFst,
    limitIndsToPlot,
    plotPCA,
    chooseChrRegion,
    plotGenotypeByIndividual,
    plotGenotypeByIndividualWithFst,
    getRollingMean,
    getWindowedFst,
    plotGenomeFst,
    getWindowedIndHet,
    standardizeIndHet,
    getVarWindowedIndHet,
    getWindowBoundaries

using MultivariateStats
using CairoMakie

"""
    greet_GenomicDiversity()

A simple function for testing that the module has loaded.
"""
function greet_GenomicDiversity()
    greeting = "Hello I am the GenomicDiversity package!"
    println(greeting)
    return greeting
end


# function copied from IrwinLabGenomicsAnalysisScriptV2.jl :
# calculate sample size (of individuals) and frequency of alternate allele. 
# genoData is Matrix where rows are individuals and columns are loci.
# indGroup is a vector of group names indicating the group each individual belongs to.
# groupsToCalc is a vector of group names for which freqs and sample sizes will be calculated.
# (homozygote alt coded as 2; heterozygote as 1, ref homozygote as 0, missing as -1 or missing)
"""
    getFreqsAndSampleSizes(genoData, indGroup, groupsToCalc)

Calculate allele frequencies and sample sizes for each group and SNP.

​# Arguments
- `genoData`: The genotype matrix, where rows are individuals and columns are loci, with genotype codes 0,1,2 meaning homozygous reference, heterozygote, homozygous alternate, and missing genotypes can be either -1 or `missing`.
- `indGroup`: A vector providing the group name each individual belongs to.
- `groupsToCalc`: A list of group names to include in calculations.

# Notes
Returns a tuple containing 1) a matrix of frequencies, and 2) a matrix of samples sizes (in both, rows are groups and columns are loci). 
"""
function getFreqsAndSampleSizes(genoData, indGroup, groupsToCalc)
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
    getSitePi(freqs, sampleSizes)

Calculate non-biased nucleotide diversity (pi) at each site for each population, with correction for sample size, using allele freqs and sample size as input.

​# Arguments
- `freqs`: Matrix containing alternate allele frequencies for each group (row) and SNP (column).
- `sampleSizes`: Matrix containing sample sizes for each group (row) and SNP (column).
"""
function getSitePi(freqs, sampleSizes)
    sitePi = @. 2freqs*(1 - freqs) * (2sampleSizes / (2sampleSizes - 1)) # the @. macro broadcasts the dot to each math operator in the formula
    # change infinites (when sample_size is somehow 1) to NaN:
    replace(sitePi, Inf => NaN)
    return sitePi
end


"""
    getRegionPi(sitePiMatrix)

Given a matrix representing per-group (rows) per-site (columns) Pi values (within-group nucleotide diversity), produce a vector containing the mean pi for each group.
"""
function getRegionPi(sitePiMatrix)
    meanPi = fill(-9.0, size(sitePiMatrix, 1))
    for i in axes(sitePiMatrix, 1)
        # get real values (ignore NaN):
        nums = sitePiMatrix[i, .!isnan.(sitePiMatrix[i,:])]
        # get mean of those real values:
        meanPi[i] = sum(nums) / length(nums)
    end
    return meanPi
end


# function copied from IrwinLabGenomicsAnalysisScriptV2.jl :
"""
    getPairwiseNames(groupsToCalc)::Vector{String}

Using a list of names, return a list of paired names in format "name1_name2".

​# Arguments
- `groupsToCalc`: The list of names.

# Notes
Returns a vector of paired names.
"""
function getPairwiseNames(groupsToCalc)::Vector{String}
    groupCount = length(groupsToCalc)
    pairwiseNames = []
    for i in 1:(groupCount-1)
        for j in (i+1):groupCount
            push!(pairwiseNames, string(groupsToCalc[i], "_", groupsToCalc[j]))
        end
    end
    return pairwiseNames
end


"""
    getDxy(freqs, groupsToCalc)

Calculate per-locus Dxy (pi_between) for each pair of populations.

# Arguments
- `freqs`: Matrix containing alternate allele frequencies for each group (row) and SNP (column).
- `groupsToCalc`: Vector containing the names of the groups corresponding to each row of `freqs`.

# Notes
Returns a tuple containing: 
- Matrix of Dxy values for each comparison (row) and locus (column).
- Vector of comparison names.
"""
function getDxy(freqs, groupsToCalc)
    pairwiseNames = getPairwiseNames(groupsToCalc)
    groupCount = length(groupsToCalc)
    Dxy = Array{Float32, 2}(undef, length(pairwiseNames), size(freqs, 2))
    rowCount = 1
    for i in 1:(groupCount-1)
        for j in (i+1):groupCount
            Dxy[rowCount,:] = @. freqs[i,:]*(1-freqs[j,:]) + freqs[j,:]*(1-freqs[i,:])
            rowCount += 1 # adds one to rowCount
        end
    end
    return Dxy, pairwiseNames
end


"""
    getRegionDxy(siteDxyMatrix)

Given a matrix representing per-population-pair (rows) per-site (columns) Dxy values (within-group nucleotide diversity), produce a vector containing the mean Dxy for each population pair.
"""
function getRegionDxy(siteDxyMatrix)
    meanDxy = fill(-9.0, size(siteDxyMatrix, 1))
    for i in axes(siteDxyMatrix, 1)
        # get real values (ignore NaN):
        nums = siteDxyMatrix[i, .!isnan.(siteDxyMatrix[i,:])]
        # get mean of those real values:
        meanDxy[i] = sum(nums) / length(nums)
    end
    return meanDxy
end


# function copied from IrwinLabGenomicsAnalysisScriptV2.jl :
"""
    getFst(freqs, sampleSizes, groupsToCalc; among=false)

Calculate Fst for each locus, between pairs of groups and (optionally) among all groups. 
    
Fst is calculated according to Weir & Cockerham (1984), including their correction for sample size and number of groups.

​# Arguments
- `freqs`: Matrix containing alternate allele frequencies for each group (row) and SNP (column).
- `sampleSizes`: Matrix containing sample sizes for each group (row) and SNP (column).
- `groupsToCalc`: The list of group names.
- `among`: Optional argument--set to "true" to also calculate Fst among all groups.

# Notes
Returns a tuple containing: 
- Matrix of Fst values for each comparison (row) and locus (column).
- Matrix of numerator values from the Fst calculation for each comparison (row) and locus (column).
- Matrix of denominator values from the Fst calculation for each comparison (row) and locus (column).
- Vector of comparison names.
"""
function getFst(freqs, sampleSizes, groupsToCalc; among=false)
    pairwiseNames = getPairwiseNames(groupsToCalc)
    groupCount = length(groupsToCalc)
    Fst = Array{Float32,2}(undef, length(pairwiseNames), size(freqs, 2))
    FstNum = Array{Float32,2}(undef, length(pairwiseNames), size(freqs, 2))
    FstDen = Array{Float32,2}(undef, length(pairwiseNames), size(freqs, 2))
    rowCount = 1
    for i in 1:(groupCount-1)
        for j in (i+1):groupCount
            r = 2 # number of populations in comparison (from Weir & Cockerham 1984 p. 1363)
            # the "@." macro adds the dot operator (broadcasting) to all function calls to the right
            n_bar = @. (sampleSizes[i, :] + sampleSizes[j, :]) / 2 # mean sample size 
            p_bar = @. (sampleSizes[i, :] * freqs[i, :] + sampleSizes[j, :] * freqs[j, :]) / (2 * n_bar) # mean frequency
            C = @. sqrt(((sampleSizes[i, :] - n_bar)^2 + (sampleSizes[j, :] - n_bar)^2) / (r - 1)) / n_bar  # The CV of sample sizes, based on the equation for sample SD--Mike & I solved the n_c equation in Weir&Cockerham1984 and confirmed this. (see the division by 1 in the numerator, prior to square root) 
            s2 = @. (sampleSizes[i, :] * ((freqs[i, :] - p_bar)^2) + sampleSizes[j, :] * ((freqs[j, :] - p_bar)^2)) / ((r - 1) * n_bar)   # allele frequency variance, using sample variance, as per W&C84
            FstNum[rowCount, :] = @. s2 - ((p_bar * (1 - p_bar) - ((r - 1) / r) * s2) / (2 * n_bar - 1))
            FstDenTerm1 = @. (1 - ((2 * n_bar * C^2) / ((2 * n_bar - 1) * r))) * p_bar * (1 - p_bar)
            FstDenTerm2 = @. (1 + ((2 * n_bar * (r - 1) * C^2) / ((2 * n_bar - 1) * r))) * s2 / r
            FstDen[rowCount, :] = FstDenTerm1 .+ FstDenTerm2
            Fst[rowCount, :] = FstNum[rowCount, :] ./ FstDen[rowCount, :]
            rowCount += 1 # adds one to rowCount
        end
    end
    if among == true
        # Add a line for the "Fst_among" for each site, 
        # among all the groups listed in "groups" variable above (so r can be greater than 2), 
        # using the Weir&Cockerham 1984 approach to correct for sample size.
        Fst_among = Array{Float32,2}(undef, 1, size(freqs, 2))
        FstNum_among = Array{Float32,2}(undef, 1, size(freqs, 2))
        FstDen_among = Array{Float32,2}(undef, 1, size(freqs, 2))
        for i in axes(freqs, 2)  # axes(A, dim) function gets the indices of A along dimension dim (recommended over 1:size(freqs,2))
          popFreqs = freqs[:,i]
          popFreqs_NaN_removed = popFreqs[.!isnan.(popFreqs)] # finds only the real (non-NaN) values
          r = length(popFreqs_NaN_removed)    # "r" in Weir&Cockerham; the number of pops with non-NaN allele freqs
          popSampleSizes = sampleSizes[:,i]
          popSampleSizesPositive = popSampleSizes[.!isnan.(popFreqs)]
          n_bar = sum(popSampleSizesPositive) / length(popSampleSizesPositive) # mean sample size
          p_bar = sum(popFreqs_NaN_removed) / length(popFreqs_NaN_removed) # mean frequency
          CV_sample_size = sqrt(sum((popSampleSizesPositive .- n_bar).^2) / (r - 1)) / n_bar     # Note that used the equation for sample SD here--Mike Whitlock & I solved the n_c equation in Weir&Cockerham1984 and confirmed this.
          allele_freq_variance = sum(popSampleSizesPositive .* ((popFreqs_NaN_removed .- p_bar).^2)) / ((r - 1)*n_bar)    # using sample variance, as per W&C84
          FstNum_among[1,i] = allele_freq_variance - ((p_bar*(1-p_bar) - ((r-1)/r)*allele_freq_variance)/(2*n_bar - 1))
          FstDen_among_term1 = (1 - (2*n_bar*((CV_sample_size)^2)/((2*n_bar - 1)*r))) * p_bar*(1-p_bar) 
          FstDen_among_term2 = (1 + ((2*n_bar*(r-1)*(CV_sample_size^2))/((2*n_bar-1)*r))) * allele_freq_variance / r
          FstDen_among[1,i] = FstDen_among_term1 + FstDen_among_term2
          Fst_among[1,i] = FstNum_among[1,i] / FstDen_among[1,i]		
        end
        Fst = vcat(Fst, Fst_among)
        FstNum = vcat(FstNum, FstNum_among)
        FstDen = vcat(FstDen, FstDen_among)
        push!(pairwiseNames, "Fst_among")
    end
    return Fst, FstNum, FstDen, pairwiseNames
end


"""
    limitIndsToPlot(plotGroups, numIndsToPlot, genoData, indMetadata)

Select a subset of individuals, based on numbers within each group, for subsequent analysis.

​# Arguments
- `plotGroups`: Vector of names of groups to include.
- `numIndsToPlot`: Vector of maximum number of individuals to include from each group (order as in `plotGroups`).
- `genoData`: Matrix of genotypes (rows are individuals, columns are loci).
- `indMetadata`: Matrix containing metadata for individuals (make sure there is an `Fst_group` column).
- `sortByMissing`: Set to true to choose individuals with the least missing data among these loci.

# Notes
Returns a tuple containing:
- Matrix of genotypes for included individuals.
- Matrix of metadata for included individuals.
"""
function limitIndsToPlot(plotGroups, numIndsToPlot, genoData, indMetadata;
                            sortByMissing = false)
    if sortByMissing
        missingGenos = (genoData .== -1) .| ismissing.(genoData)
        missingGenosPerInd = vec(sum(missingGenos, dims=2))
        rankedOrderOfMissingPerInd = sortperm(missingGenosPerInd)
        genoData = genoData[rankedOrderOfMissingPerInd, :]
        indMetadata = indMetadata[rankedOrderOfMissingPerInd, :]
    end
    cumulativeRowSelection = []
    for i in eachindex(plotGroups)
        rowSelection = findall(indMetadata.Fst_group .== plotGroups[i])
        if length(rowSelection) > numIndsToPlot[i]
            rowSelection = rowSelection[1:numIndsToPlot[i]]
        end
        cumulativeRowSelection = vcat(cumulativeRowSelection, rowSelection)
    end
    genoData_included = genoData[cumulativeRowSelection, :]
    indMetadata_included = indMetadata[cumulativeRowSelection, :]
    return genoData_included, indMetadata_included
end


"""
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
function plotPCA(genotypes, indMetadata, groups_to_plot_PCA, group_colors_PCA; 
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



# Option to focus on a region of chromosome.
# If not specified, this will set appropriate positionMin and positionMax for a chromosome
# and will make a regionText string.
# This will NOT actually filter SNPs to this range
# (that is done within the plotGenotypeByIndividual function below).
"""
    chooseChrRegion(pos, chr; positionMin=1, positionMax=NaN)

Designate a specific region of a chromosome for subsequent analysis.

​# Arguments
- `pos`: A matrix containing genomic location of each locus; must have a `chrom` column (name of scaffold) and a `position` column (numerical position on scaffold).
- `chr`: The name of the chosen scaffold.
- `positionMin`: Optional (default 1); the starting location of the chosen region.
- `positionMax`: Optional (default the end of the scaffold); the ending location of the chosen region.

# Notes
Returns a tuple containing:
- the scaffold name.
- the minimum position.
- the maximum position.
- a string describing the region (e.g. `chr Z: 20000 to 1000000`)

This function does not actually filter the genotype matrix--that can be done by providing the output of this to functions such as `plotGenotypeByIndividual()`.
"""
function chooseChrRegion(pos, chr; positionMin=1, positionMax=NaN)
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
    plotGenotypeByIndividual(regionInfo, pos, genoData, 
                                indMetadata, freqs,
                                plotGroups, plotGroupColors;
                                missingFractionAllowed = 1.0,
                                colorAllelesByGroup = true, group1 = plotGroups[1],
                                indFontSize=10, figureSize=(1200, 1200),
                                show_SNP_density = true, densityPlotColor = :steelblue1,
                                plotTitle = nothing,
                                indColorLeftProvided = false,
                                indColorRightProvided = false)

Construct a genotype-by-individual plot, with option to filter out SNPs with too much missing data. 

Under the default setting, alleles are colored (dark purple vs. light purple) according to whichever allele is designated as `group1`. 

​# Arguments

- `regionInfo`: Information regarding the scaffold and region of focus; a tuple provided by `chooseChrRegion()`.
- `pos`: Matrix providing genomic location of each locus; there must be a `chrom` column and a `position` column.
- `genoData`: Matrix containing genotype data (individuals in rows, loci in columns).
- `indMetadata`: Matrix of metadata for individuals; must contain `Fst_group` and `plot_order` columns.
- `freqs`: Matrix of alternate allele frequencies for each group (row) and locus (column).
- `plotGroups`: Vector of group names to include in plot.
- `plotGroupColors`: Vector of plotting colors corresponding to the groups.
- `missingFractionAllowed`: Optional; The maximum missing genotype fraction for a locus to be included (default is 1.0, so no filtering).
- `colorAllelesByGroup`: Optional; set to `false` to color alleles according to reference and alternate.
- `group1`: Optional (default is `plotGroups[1]`); when `colorAllelesByGroup` is `true`, this is the group that determine which allele is dark purple. 
- `indFontSize`: Optional; the font size of the individual ID labels.
- `figureSize`: Optional; the size of the figure; default is `(1200, 1200)`.  
- `show_SNP_density`: Optional; default is `true` to show a density plot. 
- `densityPlotColor`: Optional; default is `:steelblue1`
- `plotTitle`: Optional; default will make a title. For no title, set to `""`.
- `indColorLeftProvided`: Optional; Default is `false`. Set to `true` if there is a column labeled `indColorLeft` in the metadata providing color of each individual for plotting on left side.
- `indColorRightProvided`: Optional; same as above but for right side (requires `indColorRight` column in metadata).

# Notes
Returns a tuple containing:
- the figure
- the plotted genotypes
- the numerical positions (in the chosen scaffold) of the plotted loci
- the sorted metadata matrix for the plotted individuals
"""
function plotGenotypeByIndividual(regionInfo, pos, genoData, 
                                    indMetadata, freqs,
                                    plotGroups, plotGroupColors;
                                    missingFractionAllowed = 1.0,
                                    colorAllelesByGroup = true, group1 = plotGroups[1],
                                    indFontSize=10, figureSize=(1200, 1200),
                                    show_SNP_density = true, densityPlotColor = :steelblue1,
                                    plotTitle = nothing,
                                    indColorLeftProvided = false,
                                    indColorRightProvided = false)

    chr, positionMin, positionMax, regionText = regionInfo
    if isnothing(plotTitle)
        plotTitle = string(regionText, ": genotypes of loci with < ", missingFractionAllowed, " missing")
    end

    # write ≥ using \ge tab, and ≤ using \le tab.
    # note the operator ".&" (bitwise &) is much more efficient than ".&&"
    # Select based on the genome region:
    selection_position = (pos.chrom .== chr) .&
                (pos.position .≥ positionMin) .&
                (pos.position .≤ positionMax)

    SNP_genotypes = genoData[:, selection_position]
    SNP_freqs = freqs[:, selection_position]
    SNP_positions = pos.position[selection_position]
    SNP_genotypes_subset = SNP_genotypes[indMetadata.Fst_group .∈ Ref(plotGroups), :]
    indMetadata_subset = indMetadata[indMetadata.Fst_group .∈ Ref(plotGroups), :]
    # The section below does the default of coloring the allele most common in group1 as the dark purple.
    # Otherwise, set colorAllelesByGroup=false to have the alleles colored according to ref vs. alternate 
    if colorAllelesByGroup
        altAlleleHiInGroup1 = SNP_freqs[findfirst(plotGroups .== group1), :] .> 0.5
        SNP_genotypes_subset[:, altAlleleHiInGroup1] = 2 .- SNP_genotypes_subset[:, altAlleleHiInGroup1]
        SNP_genotypes_subset[SNP_genotypes_subset .== 3] .= -1  # convert original -1 values back to -1
    end
    # Choose sorting order by plot_order column in input metadata file
    #sorted.SNP.genotypes.subset = SNP.genotypes.subset[order(SNP.genotypes.subset$group, SNP.genotypes.subset$ID),]
    sorted_SNP_genotypes_subset = SNP_genotypes_subset[sortperm(indMetadata_subset.plot_order, rev=false), :]
    numInds = size(sorted_SNP_genotypes_subset, 1) # fast way to get number of rows
    sorted_indMetadata_subset = indMetadata_subset[sortperm(indMetadata_subset.plot_order, rev=false), :]

    # filter out the SNPs that have too much missing data:
    numberMissing = sum(isequal.(sorted_SNP_genotypes_subset, -1) .| ismissing.(sorted_SNP_genotypes_subset), dims=1)
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
        titlesize=30,
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
        if indColorLeftProvided
            boxColorLeft = sorted_indMetadata_subset.indColorLeft[i]
        else
            boxColorLeft = plotGroupColors[findfirst(plotGroups .== sorted_indMetadata_subset.Fst_group[i])]
        end
        if indColorRightProvided
            boxColorRight = sorted_indMetadata_subset.indColorRight[i]
        else
            boxColorRight = plotGroupColors[findfirst(plotGroups .== sorted_indMetadata_subset.Fst_group[i])]
        end
        CairoMakie.poly!(Point2f.(groupColorBox_x_left, (y .+ groupColorBox_y)), color=boxColorLeft)
        CairoMakie.poly!(Point2f.(groupColorBox_x_right, (y .+ groupColorBox_y)), color=boxColorRight)
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
        hom_ref_locs = findall(genotypes .== 0)
        if length(hom_ref_locs) > 0
            for j in eachindex(hom_ref_locs)
                CairoMakie.poly!(Point2f.((hom_ref_locs[j] .+ box_x), (y .+ box_y)), color=genotypeColors[1])
                #polygon(hom.ref.locs[j]+symbol.x, y+symbol.y, border=NA, col=genotype.colors[1])
            end
        end
        het_locs = findall(genotypes .== 1)
        if length(het_locs) > 0
            for j in eachindex(het_locs)
                CairoMakie.poly!(Point2f.((het_locs[j] .+ triangle1_x), (y .+ triangle1_y)), color=genotypeColors[1])
                CairoMakie.poly!(Point2f.((het_locs[j] .+ triangle2_x), (y .+ triangle2_y)), color=genotypeColors[3])
            end
        end
        hom_alt_locs = findall(genotypes .== 2)
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
    plotGenotypeByIndividualWithFst(groupsToCompare, Fst_cutoff, missingFractionAllowed,
                            regionInfo, pos, Fst, pairwiseNamesFst,
                            genoData, indMetadata, freqs, plotGroups, plotGroupColors;
                            colorAllelesByGroup=true, group1=plotGroups[1])

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
- `highlightRegionStarts`: Optional; the left locations of regions to highlight on the scaffold.
- `highlightRegionEnds`: Optional; the right locations of regions to highlight on the scaffold. 
- `highlightRegionColor`: Optional; default red; the color of the highlight bar along the scaffold

# Notes
Returns a tuple containing:
- the figure
- the plotted genotypes
- the numerical positions (in the chosen scaffold) of the plotted loci
- the sorted metadata matrix for the plotted individuals
"""
function plotGenotypeByIndividualWithFst(groupsToCompare, Fst_cutoff, missingFractionAllowed,
                            regionInfo, pos, Fst, pairwiseNamesFst,
                            genoData, indMetadata, freqs, plotGroups, plotGroupColors;
                            colorAllelesByGroup = true, group1 = plotGroups[1],
                            indFontSize=10, figureSize=(1200, 1200),
                            show_SNP_density = true, densityPlotColor = :steelblue1,
                            plotTitle = nothing,
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
    indMetadata_subset = indMetadata[indMetadata.Fst_group .∈ Ref(plotGroups), :]
    # The section below does the default of coloring the allele most common in group1 as the dark purple.
    # Otherwise, set colorAllelesByGroup=false to have the alleles colored according to ref vs. alternate 
    if colorAllelesByGroup
        altAlleleHiInGroup1 = SNP_freqs[findfirst(plotGroups .== group1), :] .> 0.5
        SNP_genotypes_subset[:, altAlleleHiInGroup1] = 2 .- SNP_genotypes_subset[:, altAlleleHiInGroup1]
        SNP_genotypes_subset[SNP_genotypes_subset.==3] .= -1  # convert original -1 values back to -1
    end
    # Choose sorting order by plot_order column in input metadata file
    #sorted.SNP.genotypes.subset = SNP.genotypes.subset[order(SNP.genotypes.subset$group, SNP.genotypes.subset$ID),]
    sorted_SNP_genotypes_subset = SNP_genotypes_subset[sortperm(indMetadata_subset.plot_order, rev=false), :]
    numInds = size(sorted_SNP_genotypes_subset, 1) # fast way to get number of rows
    sorted_indMetadata_subset = indMetadata_subset[sortperm(indMetadata_subset.plot_order, rev=false), :]

    # filter out the SNPs that have too much missing data:
    numberMissing = sum(isequal.(sorted_SNP_genotypes_subset, -1) .| ismissing.(sorted_SNP_genotypes_subset), dims=1)
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
        titlesize=30,
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
        hom_ref_locs = findall(genotypes .== 0)
        if length(hom_ref_locs) > 0
            for j in eachindex(hom_ref_locs)
                CairoMakie.poly!(Point2f.((hom_ref_locs[j] .+ box_x), (y .+ box_y)), color=genotypeColors[1])
                #polygon(hom.ref.locs[j]+symbol.x, y+symbol.y, border=NA, col=genotype.colors[1])
            end
        end
        het_locs = findall(genotypes .== 1)
        if length(het_locs) > 0
            for j in eachindex(het_locs)
                CairoMakie.poly!(Point2f.((het_locs[j] .+ triangle1_x), (y .+ triangle1_y)), color=genotypeColors[1])
                CairoMakie.poly!(Point2f.((het_locs[j] .+ triangle2_x), (y .+ triangle2_y)), color=genotypeColors[3])
            end
        end
        hom_alt_locs = findall(genotypes .== 2)
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
    getRollingMean(inputVector, windowSize)

Calculate non-overlapping windowed arithmetic average (i.e., means) from a vector of numeric values.

​# Arguments
- `inputVector`: A vector of numeric values (can include NaN).
- `windowSize`: An integer indicating the window number of values to include in each mean.

# Notes
- Returns a vector of windowed means. 
- In the calculations, NaN values are ignored.
- If all values in a window are NaN, then the mean is NaN.
"""
function getRollingMean(inputVector, windowSize)
    numWindows = length(inputVector) ÷ windowSize  # this is the quotient (rounds to integer)
    rollingMeans = Vector{Float64}(undef, numWindows) # set up place to store data
    for i in 1:numWindows
        start = 1 + (i-1)*windowSize
        windowValues = inputVector[start : (start + windowSize - 1)]  # this ignores the NaN values
        windowValues_real = windowValues[.!isnan.(windowValues)]
        if length(windowValues_real) > 0
            rollingMeans[i] = sum(windowValues_real) / length(windowValues_real)
        else
            rollingMeans[i] = NaN # if there are no non-NaN values in the inputVector
        end
    end
    return rollingMeans
end


"""
    getWindowedFst(FstNumerator, FstDenominator, pos, windowSize)

Calculate windowed Fst along chromosomal regions, uses output from each site's Fst (produced by the getFst() function). 

This is calculated according to Weir&Cockerham1984 (with sample size and pop number correction), calculated as windowed numerator over windowed denominator, in whole windows starting on left side of the region.

​# Arguments
- `FstNumerator`: Produced by the getFst() function, this is a vector of the Fst numerator values for each locus. 
- `FstDenominator`: Produced by the getFst() function, this is a vector of the Fst denominator values for each locus. 
- `pos`: Matrix providing genomic location of each locus; there must be a `position` column.
- `windowSize`: An integer indicating the window number of values to include in each mean.

# Notes
- returns a tuple of two vectors: rolling mean position, and rolling mean Fst
- the calculations ignore NaN values (which would be in both numerator and denominator)
- returns NaN for a window if there are no real values in the window 
"""
function getWindowedFst(FstNumerator, FstDenominator, pos, windowSize)
    rollingMeanPos = getRollingMean(pos.position, windowSize) # get per-window mean SNP position along scaffold
    rollingMeanFst = Array{Float64, 2}(undef, size(FstNumerator, 1), length(rollingMeanPos))
    for i in axes(FstNumerator, 1)
        rollingMeanFstNumerator = getRollingMean(FstNumerator[i,:], windowSize)
        rollingMeanFstDenominator = getRollingMean(FstDenominator[i,:], windowSize)
        rollingMeanFst[i,:] = rollingMeanFstNumerator ./ rollingMeanFstDenominator
    end
    return rollingMeanPos, rollingMeanFst
end


"""
    plotGenomeFst(scaffolds_to_plot, 
                    windowed_Fst_all,
                    pairwiseNamesFst,
                    windowed_pos_all,
                    groupsToPlotFst,
                    groupColorsFst;
                    lineTransparency = 0.8,
                    fillTransparency = 0.2,
                    figureSize = (1200, 1200))

Plot windowed Fst along chromosomal regions, for multiple population comparisons.

​# Arguments
- `scaffolds_to_plot`: a vector of the scaffold names to include in the plot.
- `windowed_Fst_all`: a matrix of Fst values (rows representing different population comparisons; columns different genomic regions).
- `pairwiseNamesFst`: a vector of names indicating the pairs of populations in rows of `windowed_Fst_all`.
- `windowed_pos_all`: dataframe containing genomic locations of windows (must have columns `chrom` and `mean_position`).
- `groupsToPlotFst`: a vector of strings indicating the population comparisons to include.
- `groupColorsFst`: the colors to use to plot the comparisons in `groupsToPlotFst`
- `lineTransparency` (optional): transparency of the lines
- `fillTransparency` (optional): transparency of the fill color
- `figureSize` (optional): size of the figure (hight, width)

# Notes
- Returns the figure.
"""
function plotGenomeFst(scaffolds_to_plot,
    windowed_Fst_all,
    pairwiseNamesFst,
    windowed_pos_all,
    groupsToPlotFst,
    groupColorsFst;
    lineTransparency=0.8,
    fillTransparency=0.2,
    figureSize=(1200, 1200))

    # get sizes of chromosomes (note this actually isn't the true length, just the mean position of the rightmost window):
    scaffold_lengths = repeat([-1], length(scaffolds_to_plot))# vector of Int64
    for i in eachindex(scaffolds_to_plot)
        selection = windowed_pos_all.chrom .== scaffolds_to_plot[i]
        if sum(selection) == 0
            scaffold_lengths[i] = 0
        else
            scaffold_lengths[i] = round(Int, maximum(windowed_pos_all.mean_position[selection]))
        end
    end

    scaffoldPlottedAlready = falses(length(scaffolds_to_plot))
    bpPerRow = 200_000_000
    bpGapBetweenScaffolds = 10_000_000
    scaffoldInRow = Vector{Vector{String}}() # declare a vector of vectors to store the chromosome names per row
    scaffoldBpStartInRow = Vector{Vector{Int}}()
    row = 1

    while sum(scaffoldPlottedAlready .== false) > 0  # repeat until all chr have a row to be plotted in
        bpStart = 0
        remainingRowLength = bpPerRow
        push!(scaffoldInRow, []) # initialize an empty vector for the next row
        push!(scaffoldBpStartInRow, []) # initialize an empty vector for the next row
        placeInRow = 1
        # look through chromosomes and find some to go in this row
        for i in eachindex(scaffolds_to_plot)
            # if not plotted and short enough, add to row:
            if scaffoldPlottedAlready[i] .== false .&& scaffold_lengths[i] .<= remainingRowLength
                push!(scaffoldInRow[row], scaffolds_to_plot[i])
                push!(scaffoldBpStartInRow[row], bpStart)
                remainingRowLength = remainingRowLength - scaffold_lengths[i] - bpGapBetweenScaffolds
                scaffoldPlottedAlready[i] = true
                placeInRow += 1
                bpStart = bpStart + scaffold_lengths[i] + bpGapBetweenScaffolds
            end
        end
        row += 1
    end
    println(scaffoldInRow)

    # Plot in the order defined above
    f = CairoMakie.Figure(size=figureSize)
    numRows = length(scaffoldInRow)
    # make one axis object per row:
    axs = [Axis(f[i, 1], limits=(0, bpPerRow,
            0, 1.15), ylabel="Fst") for i in 1:numRows]
    for i in 1:numRows
        hidexdecorations!(axs[i]) #hide background lattice and axis labels
        hidespines!(axs[i], :t, :r, :b) # hide box around plot
    end

    for rowNum in 1:numRows
        numScaffoldsInRow = length(scaffoldInRow[rowNum])
        nextPlotShiftBP = 0
        for orderNum in 1:numScaffoldsInRow
            scaffoldName = scaffoldInRow[rowNum][orderNum]
            plotBpLength = scaffold_lengths[findfirst(scaffolds_to_plot .== scaffoldName)]
            # select one scaffold from the whole genome:
            selection = windowed_pos_all.chrom .== scaffoldName
            windowed_pos_selected = view(windowed_pos_all.mean_position, selection)
            windowed_Fst_selected = view(windowed_Fst_all, :, selection)
            if length(windowed_pos_selected) > 0  # only plot if there is data for that scaffold        
                text!(axs[rowNum], scaffoldBpStartInRow[rowNum][orderNum] + 0.005bpPerRow, 1.05; text=scaffoldName, align=(:left, :center), fontsize=20)
                for j in eachindex(groupsToPlotFst)
                    xx = vcat(first(windowed_pos_selected), windowed_pos_selected, last(windowed_pos_selected)) .+ scaffoldBpStartInRow[rowNum][orderNum]
                    FstRow = findfirst(pairwiseNamesFst .== groupsToPlotFst[j])
                    yy = vcat(0, windowed_Fst_selected[FstRow, :], 0)
                    poly!(axs[rowNum], Point2f.(xx, yy), color=(groupColorsFst[j], fillTransparency), strokecolor=(groupColorsFst[j], lineTransparency), strokewidth=1)
                end
            end
        end
    end
    display(f)
    return f
end


"""
    getWindowedIndHet(genotypes, pos, windowSize)

Calculate windowed heterozygosity for each individual, across a genomic region.

​# Arguments
- `genotypes`: The genotype matrix, where rows are individuals and columns are loci, with genotype codes 0,1,2 meaning homozygous reference, heterozygote, homozygous alternate, and missing genotypes can be either -1 or `missing`.
- `pos`: A matrix containing genomic location of each locus; must have a `position` column (numerical position on scaffold). 
- `windowSize`: An integer indicating the number of loci to include in each mean.

# Notes
Returns a tuple containing 1) a vector of windowed mean positions, and 2) a matrix of windowed heterozgyosity values (rows are individuals and columns are loci). 
"""
function getWindowedIndHet(genotypes, pos, windowSize)
    rollingMeanPos = getRollingMean(pos.position, windowSize) # get per-window mean SNP position along scaffold
    # convert both homozygotes to value zero, leaving heterozgyotes value one;
    # also covert missing and -1 values to NaN, which getRollingMean function deals with
    heterozygotes = replace(genotypes, 2 => 0, missing => NaN, -1 => NaN)

    rollingMeanHet = Array{Float64, 2}(undef, size(heterozygotes, 1), length(rollingMeanPos))
    # cycle through individuals
    for i in axes(heterozygotes, 1)
        rollingMeanHet[i,:] = getRollingMean(heterozygotes[i,:], windowSize)
    end
    return rollingMeanPos, rollingMeanHet
end


"""
    standardizeIndHet(windowedIndHet)

Divide each individual's windowed heterozgyosity by the average heterozygosity of all individuals (at each window).

​# Arguments
- `windowedIndHet`: The matrix of individual windowed heterozygosity values (rows are individuals; columns represent windowed genomic regions)

# Notes
Returns a matrix of standardized individual windowed heterozygosities.  
"""
function standardizeIndHet(windowedIndHet)
    meanPerWindow_windowedIndHet = sum(windowedIndHet, dims=1) ./ size(windowedIndHet, 1)
    windowedIndHet_standardized = Array{Float64,2}(undef, size(windowedIndHet, 1), size(windowedIndHet, 2))
    for row in eachindex(windowedIndHet[:,1])
        windowedIndHet_standardized[row, :] = windowedIndHet[row, :] ./ vec(meanPerWindow_windowedIndHet)
    end
    return windowedIndHet_standardized
end


"""
    getVarWindowedIndHet(windowedIndHet_standardized)

Calculate the windowed variance in standardized individual heterozygosity. 

​# Arguments
- `windowedIndHet_standardized`: The matrix of standardized individual windowed heterozygosity values (rows are individuals; columns represent windowed genomic regions)

# Notes
Returns a vector of variances in standardized individual heterozygosities (each value corresponds to a different windowed genomic region).  
"""
function getVarWindowedIndHet(windowedIndHet_standardized)
    variancePerWindow_windowedIndHet = vec(var(windowedIndHet_standardized; dims = 1))
end


"""
    getWindowBoundaries(lociVector, windowSize)

Determine the start and end scaffold positions for each window along a single scaffold. 

​# Arguments
- `lociVector`: Vector of included locus positions along a scaffold.
- `windowSize`: Size of each window (in number of loci).

# Notes
- returns a Matrix, where rows represent windows and columns are starting location and ending location of each window along the scaffold. 
"""
function getWindowBoundaries(lociVector, windowSize)
    numWindows = length(lociVector) ÷ windowSize  # this is the quotient (rounds to integer)
    windowBoundaries = Matrix{Int}(undef, numWindows, 2) # set up place to store data
    for i in 1:numWindows
        startIndex = 1 + (i-1)*windowSize
        endIndex = startIndex + windowSize - 1
        windowBoundaries[i, 1] = lociVector[startIndex] # this is the position of the first locus in the window
        windowBoundaries[i, 2] = lociVector[endIndex] # this is the position of the last locus in the window
    end
    return windowBoundaries
end


end # of module GenomicDiversity
