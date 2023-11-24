#= 
SNPlots.jl
This file started by Darren Irwin on 24 June 2023
as a set of functions for processing and displaying SNP data.
    
Converted to an actual package (private) with a Github repo on 11 November 2023.
    =#

module SNPlots

# make functions available in calling namespace through "using SNPlots" command, without needing to prefix module name: 
export greet_SNPlots,
    getFreqsAndSampleSizes, 
    getPairwiseNames,
    getFst,
    limitIndsToPlot,
    plotPCA,
    chooseChrRegion,
    plotGenotypeByIndividual

using MultivariateStats
using CairoMakie

"""
    greet_SNPlots()

A simple function for testing that the module has loaded.
"""
function greet_SNPlots()
    greeting = "Hello I am the SNPlots package!"
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
    genoData[ismissing.(genoData)] .= -1 # if "missing" datatype is use, convert to -1
    groupCount = length(groupsToCalc)
    freqs = Array{Float32,2}(undef, groupCount, size(genoData, 2))
    sampleSizes = Array{Int16,2}(undef, groupCount, size(genoData, 2))
    for i in 1:groupCount
        selection = (indGroup .== groupsToCalc[i]) # gets the correct rows for individuals in the group 
        geno0counts = sum(genoData[selection, :] .== 0, dims=1) # count by column the number of 0 genotypes (homozygous ref)
        geno1counts = sum(genoData[selection, :] .== 1, dims=1) # same for 1 genotypes (heterozygous)
        geno2counts = sum(genoData[selection, :] .== 2, dims=1) # same for 2 genotypes (homozygous alternate) 
        sumGenoCounts = geno0counts .+ geno1counts .+ geno2counts
        sampleSizes[i, :] = sumGenoCounts
        freqs[i, :] = ((2 .* geno2counts) .+ geno1counts) ./ (2 * sumGenoCounts)
    end
    return freqs, sampleSizes
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

# Notes
Returns a tuple containing:
- Matrix of genotypes for included individuals.
- Matrix of metadata for included individuals.
"""
function limitIndsToPlot(plotGroups, numIndsToPlot, genoData, indMetadata)
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
- `sampleSet`: Optional, for specifying a name of the sample set to appear in the plot title.
- `regionText`: Optional, for specifying a name of the genomic region to appear in the plot title.
- `flip1`: Optional, set to `true` if wanting to flip PC1 (i.e., multiply by -1).
- `flip2`: Same but for PC2.
- `showPlot`: Set to `false` if not wanting the function to draw display
- `autolimitaspect_setting`: for unconstrained axes, set to `nothing` (default is `1` for 1:1 axes)

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
                    showPlot = true, autolimitaspect_setting = 1)
    
    selection = map(in(groups_to_plot_PCA), indMetadata.Fst_group)
    matrixForPCA = Matrix{Float32}(transpose(genotypes[selection, :]))
    indMetadata_groupSelected = indMetadata[selection, :]

    if sampleSet == "" && regionText == ""
        plotTitle = "PCA"
    elseif sampleSet == ""
        plotTitle = string("PCA: ", regionText)
    elseif regionText == ""
        plotTitle = string("PCA of ", sampleSet)
    else 
        plotTitle = string("PCA of ", sampleSet, ": ", regionText)
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
            xlabel = "PC1",
            ylabel = "PC2",
            autolimitaspect = autolimitaspect_setting
        )
        for i in eachindex(groups_to_plot_PCA) 
            selection = indMetadata_groupSelected.Fst_group .== groups_to_plot_PCA[i]
            CairoMakie.scatter!(ax, PC1[selection], PC2[selection], marker = :diamond, color=group_colors_PCA[i], markersize=10, strokewidth=0.5)
        end
        showPlot && display(f)
    catch # Makie sometimes has an error due to the "autolimitaspect = 1" above, so adding this "try-catch" structure to keep program going:
        println("Did not success in drawing PCA with 1:1 axes for ", regionText,
                ", so drawing with non-proportional axes.")
        global f = CairoMakie.Figure()
        global ax = Axis(f[1, 1],
            title = plotTitle,
            xlabel = "PC1",
            ylabel = "PC2",
        )
        for i in eachindex(groups_to_plot_PCA) 
            selection = indMetadata_groupSelected.Fst_group .== groups_to_plot_PCA[i]
            CairoMakie.scatter!(ax, PC1[selection], PC2[selection], marker = :diamond, color=group_colors_PCA[i], markersize=10, strokewidth=0.5)
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
    plotGenotypeByIndividual(groupsToCompare, Fst_cutoff, missingFractionAllowed,
                            regionInfo, pos, Fst, pairwiseNamesFst,
                            genoData, indMetadata, freqs, plotGroups, plotGroupColors;
                            colorAllelesByGroup=true, group1=plotGroups[1])

Construct a genotype-by-individual plot, including only loci that pass thresholds for Fst and missing data. 

Under the default setting, alleles are colored (dark purple vs. light purple) according to whichever allele is designated as `group1`. 

​# Arguments
- `groupsToCompare`: The two groups (in format `name1_name2`) for filtering loci based on Fst.
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


# Notes
Returns a tuple containing:
- the figure
- the plotted genotypes
- the numerical positions (in the chosen scaffold) of the plotted loci
- the sorted metadata matrix for the plotted individuals
"""
function plotGenotypeByIndividual(groupsToCompare, Fst_cutoff, missingFractionAllowed,
                            regionInfo, pos, Fst, pairwiseNamesFst,
                            genoData, indMetadata, freqs, plotGroups, plotGroupColors;
                            colorAllelesByGroup=true, group1=plotGroups[1],
                            indFontSize=10, figureSize=(1200, 1200))   
    chr, positionMin, positionMax, regionText = regionInfo
    # if the genoData has missing values, then convert to -1:
    genoData[ismissing.(genoData)] .= -1
    # Filter SNPs according to Fst of the groupsToCompare
    rowChoice = findfirst(pairwiseNamesFst .== groupsToCompare)
    # write ≥ using \ge tab, and ≤ using \le tab.
    # note the operator ".&" (bitwise &) is much more efficient than ".&&"
    selection = (pos.chrom .== chr) .&
                (Fst[rowChoice, :] .≥ Fst_cutoff) .&
                (.!isnan.(Fst[rowChoice, :])) .&
                (pos.position .≥ positionMin) .&
                (pos.position .≤ positionMax)
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
        SNP_genotypes_subset[SNP_genotypes_subset.==3] .= -1
    end
    # Choose sorting order by plot_order column in input metadata file
    #sorted.SNP.genotypes.subset = SNP.genotypes.subset[order(SNP.genotypes.subset$group, SNP.genotypes.subset$ID),]
    sorted_SNP_genotypes_subset = SNP_genotypes_subset[sortperm(indMetadata_subset.plot_order, rev=false), :]
    numInds = size(sorted_SNP_genotypes_subset, 1) # fast way to get number of rows
    sorted_indMetadata_subset = indMetadata_subset[sortperm(indMetadata_subset.plot_order, rev=false), :]

    # filter out the SNPs that have too much missing data:
    numberMissing = sum(sorted_SNP_genotypes_subset .== -1 .| ismissing.(sorted_SNP_genotypes_subset), dims=1)
    fractionMissing = numberMissing / numInds
    selection = vec(fractionMissing .≤ missingFractionAllowed)
    sorted_SNP_genotypes_subset2 = sorted_SNP_genotypes_subset[:, selection]
    SNP_positions_subset2 = SNP_positions[selection]
    num_SNPs_to_plot = length(SNP_positions_subset2)

    # Set up the plot window:
    f = CairoMakie.Figure(resolution=figureSize)
    ax = Axis(f[1, 1],
        title=string(regionText, ": genotypes Fst>", Fst_cutoff, " loci between ", groupsToCompare),
        # xlabel = "location",
        # ylabel = "individual",
        titlesize=30,
        limits=(0.5 - 0.09 * (num_SNPs_to_plot + 0.5), 1.09 * (num_SNPs_to_plot + 0.5),
            0.5 - 0.3 * numInds, numInds + 1)
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

    # make lower part of figure (indicating position along chromosome) #NOTE THAT CHROMOSOME LENGTH NOT QUITE THE TRUE LENGTH
    chrLine_y = 0.5 - 0.2numInds
    topHatchLine_y1 = 0.6
    topHatchLine_y2 = 0.5 - 0.02numInds
    lowHatchLine_y1 = 0.5 - 0.18numInds
    lowHatchLine_y2 = 0.5 - 0.2numInds
    # draw chromosome line and text labels
    CairoMakie.lines!([0.5, num_SNPs_to_plot + 0.5], [chrLine_y, chrLine_y], color="black", linewidth=5)
    CairoMakie.text!(0.5 + (num_SNPs_to_plot / 2), chrLine_y - 0.05numInds; text=string("Location along chromosome ", chr), align=(:center, :center), fontsize=30)
    CairoMakie.text!(0.5, chrLine_y - 0.025numInds; text=string(positionMin), align=(:center, :center), fontsize=20)
    CairoMakie.text!(0.5 + num_SNPs_to_plot, chrLine_y - 0.025numInds; text=string(positionMax), align=(:center, :center), fontsize=20)
    # draw lines from SNPs to chromosome line
    chrPlotRatio = num_SNPs_to_plot / (positionMax - positionMin)
    for i in eachindex(SNP_positions_subset2)
        CairoMakie.lines!([i, i], [topHatchLine_y1, topHatchLine_y2], color="grey20", linewidth=1)
        CairoMakie.lines!([i, 0.5 + chrPlotRatio * (SNP_positions_subset2[i] - positionMin)], [topHatchLine_y2, lowHatchLine_y1], color="grey20", linewidth=1)
        CairoMakie.lines!([0.5 + chrPlotRatio * (SNP_positions_subset2[i] - positionMin), 0.5 + chrPlotRatio * (SNP_positions_subset2[i] - positionMin)], [lowHatchLine_y1, lowHatchLine_y2], color="grey20", linewidth=1)
    end
    display(f)
    return f, sorted_SNP_genotypes_subset2, SNP_positions_subset2, sorted_indMetadata_subset
end

end # of module SNPlots
