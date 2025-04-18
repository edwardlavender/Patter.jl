using Documenter
using Patter

makedocs(
    sitename = "Patter",
    # Boost permitted html size
    # 300 * 2^10 = 300 kiB (default = 200)
    format = Documenter.HTML(size_threshold = 300 * 2^10, 
                             size_threshold_warn = 300 * 2^10), 
    modules = [Patter]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/edwardlavender/Patter.jl.git"
)
