using Documenter
using Patter

makedocs(
    sitename = "Patter",
    format = Documenter.HTML(),
    modules = [Patter]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/edwardlavender/Patter.jl.git"
)
