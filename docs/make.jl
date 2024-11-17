using Documenter
using GenomicDiversity

makedocs(
    sitename = "GenomicDiversity",
    format = Documenter.HTML(),
    modules = [GenomicDiversity]
)

deploydocs(
    repo = "https://github.com/darreni/GenomicDiversity.jl"
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
