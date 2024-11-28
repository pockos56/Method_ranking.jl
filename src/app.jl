using CSV, DataFrames, ProgressBars, Statistics, FreqTables, ScikitLearn, PyCall, Conda

"""
A function that ranks and outputs the score of chromatography methods based on their suitability for the analysis of specific compounds given their SMILES or Inchikeys.

# Examples
```julia-repl
julia> rank_methods(["CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CC(C)(C1=CC=C(C=C1)O)C2=CC=C(C=C2)O", "C1=CC=C(C=C1)C[C@@H](C(=O)O)N"])
```
"""
function rank_methods(mol_identifiers::Union{String, Vector{String}}; w_homo=100, w_first=-10, w_last=-5, w_No=-40, windows=5, sorted=true, presence_model_version="latest", window_model_version="latest", identifier = "auto")
    # Packages
    cat = pyimport("catboost")
    pcp = pyimport("pubchempy")
    pd = pyimport("padelpy")
    cd(@__DIR__)
    path_models = (pwd() * "\\data\\models\\")
    
    function from_INCHIKEY_to_fp(mol_identifiers)
        # Auto assignment to identifier argument
        if identifier == "auto"
            if all(occursin.(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]{1}$", mol_identifiers))
                identifier = "INCHIKEY"
            else identifier = "SMILES" end
        else  identifier = uppercase(identifier)
        end

        if identifier != "INCHIKEY" && identifier != "SMILES"
            error("Please use one of the following argmuents as the identifier argument: inchikey, SMILES, auto") end

        # Load Inchikeys and list of method names
        method_names_list = CSV.read(path_models*"method_names.csv", DataFrame)
        if typeof(mol_identifiers) == String        # e.g. mol_identifiers = "CCCCO"
            mol_identifiers = [mol_identifiers]
        end
        if typeof(mol_identifiers) == Vector{String}
            mol_identifiers = Vector(string.(unique(strip.(mol_identifiers)))) 
        else error("rank_methods function requires the compound info as argument. E.g. rank_methods([\"C1=CC=CC=C1\",\"CCCC\",\"CC(CC)O\"])")
        end

        # Create dataset dataframe
        dataset_df = DataFrame()
        for i = 1:size(method_names_list,1)
            dataset_df = append!(dataset_df, DataFrame("Method Name"=>method_names_list[i,"Method Name"], "ID"=>mol_identifiers)) end
        insertcols!(dataset_df, "Presence"=>"To be predicted")
        insertcols!(dataset_df, "RT Window"=>"To be predicted")

        # Create FP dictionary
        fp_dict = DataFrame()
        for comp_i in ProgressBar(mol_identifiers)    # comp_i = mol_identifiers[2]
            fp_temp = DataFrame()
            shortest_cid = []
            try
                fp_temp = []
                if identifier == "INCHIKEY"
                    shortest_cid = pcp.get_compounds(comp_i, "inchikey")[1]
                    fp_temp = DataFrame(pd.from_smiles(shortest_cid.isomeric_smiles, fingerprints=true, descriptors=false))
                    fp_temp = convert.(Int16, parse.(Float16, fp_temp))
                elseif identifier == "SMILES"
                    fp_temp = DataFrame(pd.from_smiles(comp_i, fingerprints=true, descriptors=false))
                    fp_temp = convert.(Int16, parse.(Float16, fp_temp))
                end
            catch
                error("Failure to calculate fingerprint for compound: $comp_i. \nPossible solutions: \n1. Check spelling. \n2. Visit PubChem database and update to the latest InChIKey.\n3. Remove from the list of molecular identifiers and rerun.")
            end
            if size(fp_temp,2) != 1661
                error("An incorrect type of fingerprints has been calculated. This error should never occur. Check if the file descriptors.xml is in the correct directory.") end

            fp_dict_to_add = hcat(DataFrame("ID"=>comp_i), fp_temp)
            fp_dict = append!(fp_dict, fp_dict_to_add)
        end
        #println("Fingerprints calculated for all compounds.")
        
        # Create PubChem compressed FP dictionary
        function compressPubChemFPs(ACfp::DataFrame, PCfp::DataFrame)
            FP1tr = convert.(Float64, ACfp)
            pubinfo = convert.(Int, Matrix(PCfp))
            findidx(FP_number) = findfirst(x -> x .== "PubchemFP$FP_number", names(PCfp))

            # Ring counts
            FP1tr[!, "PCFP-r3"] = pubinfo[:, findidx(115)]
            FP1tr[!, "PCFP-r3"][pubinfo[:, findidx(122)].==1] .= 2
            FP1tr[!, "PCFP-r4"] = pubinfo[:, findidx(129)]
            FP1tr[!, "PCFP-r4"][pubinfo[:, findidx(136)].==1] .= 2
            FP1tr[!, "PCFP-r5"] = pubinfo[:, findidx(143)]
            FP1tr[!, "PCFP-r5"][pubinfo[:, findidx(150)].==1] .= 2
            FP1tr[!, "PCFP-r5"][pubinfo[:, findidx(157)].==1] .= 3
            FP1tr[!, "PCFP-r5"][pubinfo[:, findidx(164)].==1] .= 4
            FP1tr[!, "PCFP-r5"][pubinfo[:, findidx(171)].==1] .= 5

            FP1tr[!, "PCFP-r6"] = pubinfo[:, findidx(178)]
            FP1tr[!, "PCFP-r6"][pubinfo[:, findidx(185)].==1] .= 2
            FP1tr[!, "PCFP-r6"][pubinfo[:, findidx(192)].==1] .= 3
            FP1tr[!, "PCFP-r6"][pubinfo[:, findidx(199)].==1] .= 4
            FP1tr[!, "PCFP-r6"][pubinfo[:, findidx(206)].==1] .= 5
            FP1tr[!, "PCFP-r7"] = pubinfo[:, findidx(213)]
            FP1tr[!, "PCFP-r7"][pubinfo[:, findidx(220)].==1] .= 2
            FP1tr[!, "PCFP-r8"] = pubinfo[:, findidx(227)]
            FP1tr[!, "PCFP-r8"][pubinfo[:, findidx(234)].==1] .= 2
            FP1tr[!, "PCFP-r9"] = pubinfo[:, findidx(241)]     #160, 127
            FP1tr[!, "PCFP-r10"] = pubinfo[:, findidx(248)]    #167, 134

            #minimum number of type of rings
            arom = zeros(size(pubinfo, 1))
            arom[(arom.==0).&(pubinfo[:, findidx(261)].==1)] .= 4
            arom[(arom.==0).&(pubinfo[:, findidx(259)].==1)] .= 3
            arom[(arom.==0).&(pubinfo[:, findidx(257)].==1)] .= 2
            arom[(arom.==0).&(pubinfo[:, findidx(255)].==1)] .= 1
            FP1tr[!, "minAromCount"] = arom
            het = zeros(size(pubinfo, 1))
            het[(het.==0).&(pubinfo[:, findidx(262)].==1)] .= 4
            het[(het.==0).&(pubinfo[:, findidx(260)].==1)] .= 3
            het[(het.==0).&(pubinfo[:, findidx(258)].==1)] .= 2
            het[(het.==0).&(pubinfo[:, findidx(256)].==1)] .= 1
            FP1tr[!, "minHetrCount"] = het

            custom_FPs_df = DataFrame(convert.(Int16, Matrix(FP1tr)), names(FP1tr))
            println("Compressed fingerprints calculation: Completed.")

            return custom_FPs_df
        end
        compressed = compressPubChemFPs(fp_dict[:,2:781], fp_dict[:,782:end])
        fp_dict = hcat(fp_dict[:,1], compressed)
        rename!(fp_dict, Dict(:x1 => "ID"))

        ## Associate the dataset_df to the fp_dict
        fp_dataset = DataFrame()
        for compound_i in (1:size(dataset_df,1))   # Delete after debugging: compound_i = 3
            fp_dict_comp = fp_dict[findfirst(x -> x.== dataset_df[compound_i,"ID"], fp_dict[:,"ID"]),:]
            fp_dataset_to_add = hcat(DataFrame(dataset_df[compound_i,:]),DataFrame(fp_dict_comp[2:end]))
            fp_dataset = append!(fp_dataset, fp_dataset_to_add)
        end
        #println("FP dataset ready.")
        return fp_dataset
    end
    function predict_presence_and_windows(fp_dataset)
        # Making a new dataframe
        dataset = deepcopy(fp_dataset)

        # Load models from file
        presence_model_available = readdir(path_models)[contains.(readdir(path_models),"presence_model_v")]
        if presence_model_version == "latest"
            presence_model_file = path_models*presence_model_available[end]
        else
            presence_model_file = path_models * "presence_model_v" * presence_model_version * ".bin"
        end

        try
            global presence_model = cat.CatBoostClassifier().load_model(presence_model_file)
        catch
            error("Available presence model versions: $([versions for versions in presence_model_available]). Enter version as string, e.g. '0.1' or 'latest'")
        end

        window_model_available = readdir(path_models)[contains.(readdir(path_models),"window_model_v")]
        if window_model_version == "latest"
            window_model_file = path_models*window_model_available[end]
        else
            window_model_file = path_models * "window_model_v" * window_model_version * ".bin"
        end

        try
            global window_model = cat.CatBoostClassifier().load_model(window_model_file)
        catch
            error("Available window model versions: $([versions for versions in window_model_available]). Enter version as string, e.g. '0.1' or 'latest'")
        end
        # Checking impossible error
        if presence_model.is_fitted() == false || window_model.is_fitted() == false 
            error("Failed loading models. This error should never occur.") end

        # Making presence predictions
        X_set = hcat(dataset[:,"Method Name"], dataset[:,end-789:end])
        dataset[:,"Presence"] = presence_model.predict(Matrix(X_set))
        dataset[:,"RT Window"] = window_model.predict(Matrix(X_set))[:]

        # The presence prediction overlaps the RT window prediction (i.e. If a compound does not appear in a specific method, there is no RT window with it.)
        dataset[findall(x->x .== "NaN", dataset[:,"Presence"]), "RT Window"] .= "No"
        println("Retention information prediction: Completed.")
        return dataset[:,1:5]
    end
    function set_score_legacy(dataset)
        scores_list = DataFrame()
        no_of_methods = length(unique(dataset[:, "Method Name"]))
    
        for method_i in (1:no_of_methods)    # method_i=2
            method_name_i = unique(dataset[:, "Method Name"])[method_i]
            method_i_idx = findall(x -> x .== method_name_i, dataset[:, "Method Name"])
            method_set = dataset[method_i_idx, :]
            RT_windows = dataset[method_i_idx, "RT Window"]

            no_of_theoretical_comps = length(RT_windows)
            no_of_present_comps = length(RT_windows[RT_windows.!="No"])
            ratios = freqtable(RT_windows)[1:end-1] / no_of_present_comps
            ratios_df = DataFrame("Window" => string.((collect('A':'Z'))[1:windows]), "Ratio" => zeros(windows))
            try
                for k = 1:length(ratios)
                    ratios_df[findfirst(x -> x .== names(ratios)[1][k], ratios_df[:, "Window"]), "Ratio"] = ratios[k]
                end    
            catch
                error("Make sure that the number of windows is set correctly.")
            end
            
            # Score variables
            z_homo = w_homo * std(ratios_df[:, "Ratio"])
            z_first = w_first * ratios_df[ratios_df.Window .== ratios_df.Window[1], "Ratio"][1]
            z_last = w_last * ratios_df[ratios_df.Window .== ratios_df.Window[end], "Ratio"][1]
            z_No = w_No * ((no_of_theoretical_comps - no_of_present_comps) / no_of_theoretical_comps)

            # Score
            score_val = z_homo + z_first + z_last + z_No
            scores_list_to_add = DataFrame("Method number"=>method_i, "Method Name"=>method_name_i, "Score"=>score_val)
            scores_list = append!(scores_list, scores_list_to_add)
        end
        scores_list_sorted = sort(scores_list, "Score", rev=true)

        println("The optimal method for the selected compounds is $(scores_list_sorted[1,"Method Name"])")
        println("The 2nd optimal method for the selected compounds is $(scores_list_sorted[2,"Method Name"])")
        println("The 3rd optimal method for the selected compounds is $(scores_list_sorted[3,"Method Name"])")
        if sorted
            return scores_list_sorted
        else return scores_list
        end
    end
    function set_score(dataset)
        scores_list = DataFrame()
        no_of_methods = length(unique(dataset[:, "Method Name"]))
    
        for method_i in ProgressBar(1:no_of_methods)    # method_i=4
            method_name_i = unique(dataset[:, "Method Name"])[method_i]
            method_i_idx = findall(x -> x .== method_name_i, dataset[:, "Method Name"])
            method_set = dataset[method_i_idx, :]
            RT_windows = dataset[method_i_idx, "RT Window"]
    
            no_of_theoretical_comps = length(RT_windows)
            no_of_present_comps = length(RT_windows[RT_windows.!="No"])
    
            ratios = freqtable(RT_windows[RT_windows .!= "No"]) / no_of_present_comps
            ratios_df = DataFrame("Window" => string.((collect('A':'Z'))[1:windows]), "Ratio" => zeros(windows))
            try
                for k = 1:length(ratios)
                    ratios_df[findfirst(x -> x .== names(ratios)[1][k], ratios_df[:, "Window"]), "Ratio"] = ratios[k]
                end    
            catch
                error("Make sure that the number of windows is set correctly.")
            end
            
            # Score variables
            z_homo = w_homo * exp(-2 * std(ratios_df[:, "Ratio"]))
            z_first = w_first * ratios_df[ratios_df.Window .== ratios_df.Window[1], "Ratio"][1]
            z_last = w_last * ratios_df[ratios_df.Window .== ratios_df.Window[end], "Ratio"][1]
            z_No = w_No * ((no_of_theoretical_comps - no_of_present_comps) / no_of_theoretical_comps)
    
            # Score
            function score_sigmoid_norm(score, n_comps = no_of_theoretical_comps)
                function sigmoid(score, steepness::Float64=0.05)
                    midpoint = 65 + 2 * log(n_comps)
                    return 100 / (1 + exp(-steepness * (score - midpoint)))
                end
                #return clamp(transformed_score, 0, 100)
                return (sigmoid(score) - sigmoid(0)) * (100 / (sigmoid(100) - sigmoid(0)))
            end
            score_val_raw = z_homo + z_first + z_last + z_No
            score_val = score_sigmoid_norm(score_val_raw)
            if isempty(ratios)
                score_val = 0.0 end
            scores_list_to_add = DataFrame("Method number"=>method_i, "Method Name"=>method_name_i, "Score"=>Float64(score_val))
            scores_list = append!(scores_list, scores_list_to_add)
        end
        scores_list_sorted = sort(scores_list, "Score", rev=true)
    
        println("The optimal method for the selected compounds is $(scores_list_sorted[1,"Method Name"])")
        println("The 2nd optimal method for the selected compounds is $(scores_list_sorted[2,"Method Name"])")
        println("The 3rd optimal method for the selected compounds is $(scores_list_sorted[3,"Method Name"])")

        if sorted
            return scores_list_sorted
        else return scores_list
        end
    end
    
    fp_dataset = from_INCHIKEY_to_fp(mol_identifiers)
    dataset = predict_presence_and_windows(fp_dataset)
    scores = set_score(dataset)
    return scores
end

# Example 1
#mol_identifiers = ["RYYVLZVUVIJVGH-UHFFFAOYSA-N", "IISBACLAFKSPIT-UHFFFAOYSA-N", "COLNVLDHVKWLRT-QMMMGPOBSA-N", "ZYGHJZDHTFUPRJ-UHFFFAOYSA-N", "CMPQUABWPXYYSH-UHFFFAOYSA-N", "MUMGGOZAMZWBJJ-DYKIIFRCSA-N", "CBCKQZAAMUWICA-UHFFFAOYSA-N", "OYGQVDSRYXATEL-UHFFFAOYSA-N", "BSYNRYMUTXBXSQ-UHFFFAOYSA-N", "UFWIBTONFRDIAS-UHFFFAOYSA-N"]
#scores = rank_methods(mol_identifiers);

# Example 2
#mol_identifiers = ["CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CC(C)(C1=CC=C(C=C1)O)C2=CC=C(C=C2)O", "C1=CC=C(C=C1)C[C@@H](C(=O)O)N "]
#scores = rank_methods(mol_identifiers);
