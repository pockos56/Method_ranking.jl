module Method_ranking

    # Packages
    using Conda
    using PyCall
    using ScikitLearn
    using CSV
    using Statistics
    using DataFrames
    using ProgressBars
    using FreqTables
    cat = pyimport("catboost")
    pcp = pyimport("pubchempy")
    pd = pyimport("padelpy")
 
    include("app.jl")

    export rank_methods
end
