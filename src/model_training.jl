## In this file we split the pest_mix into individual files that contain the masses.
## We split the measurements (different pest mix) of each method
## We assign compound names to the peaks by correlating the assigned mass to the pest_mix info.
using CSV, DataFrames, Dates, Statistics, ProgressBars, FreqTables
project_path = "C://Users//alex_//Documents//Github//Method_ranking//"
## This script splits the raw pesticide mix info files into 10 individual files. One for each mixture number
function split_into_10_individual_mixtures(full_file_dir="C://Users//alex_//Documents//Github//Method_ranking//data//AZ_raw//Pesticide mix info//", full_file_name="Pesticide RESTEK multiresidue information.csv")
    mix_list = CSV.read(full_file_dir * full_file_name, DataFrame, header=1)
    mix_list = deleteat!(mix_list, findall(ismissing.(mix_list[:, "Chemical Name"])))        # Delete the empty spaces inbetween different Pest mixes
    header_indices = findall(x -> occursin("mix", x), mix_list[:, "Chemical Name"])          # Find where the headers are (e.g. Pest mix #1)
    for i = 1:length(header_indices)
        list_start = header_indices[i] + 1
        if i != length(header_indices)
            list_end = header_indices[i+1] - 1
        else
            list_end = size(mix_list, 1)
        end
        mix_temp_list = mix_list[list_start:list_end, :]
        CSV.write(full_file_dir * "pest_mix_$i.csv", mix_temp_list)
    end
end
# split_into_10_individual_mixtures()

## This script splits the measurements (different pest mix) of each method
function identification(; samples_csv="C://Users//alex_//Documents//Github//Method_ranking//data//AZ_raw//AZ_raw_csv//", pest_dir="C://Users//alex_//Documents//Github//Method_ranking//data//AZ_raw//Pesticide mix info//", mass_col="Base Peak (Combined)", min_area_count=100000, H_adduct_mass=1.0078, mass_tolerance=100e-3)

    method_files = readdir(samples_csv)
    methods_names = [i[1:end-4] for i in method_files]
    all_methods = DataFrame()
    for method_no in ProgressBar(1:length(method_files))               # 21 methods  method_no=1
        method_file = CSV.read(samples_csv * method_files[method_no], DataFrame, header=2)
        if names(method_file)[1] != "Name"
            method_file = CSV.read(samples_csv * method_files[method_no], DataFrame, header=1)
        end
        if any(contains.(strip.(names(method_file)), "Pesticide Mix"))        # Some method_files have an extra Pesticide Mix column. It is removed here
            select!(method_file, Not("Pesticide Mix "))
        end
        if any(contains.(strip.(names(method_file)), "Injection")) == false     # Some method_files miss the Injection column. It is added here
            method_file = insertcols(method_file, 13, "Injection" => 1)
        end

        for date_i = 1:size(method_file, 1)       # Sort the measurements by date to ensure the 1-10 injections correspond to pest mix #1-10
            method_file[date_i, "Date Acquired"] = Dates.format(DateTime(method_file[date_i, "Date Acquired"][1:end-4], "EEEE, U dd, yyyy H:MM:SS p"), "yyyy-mm-dd HH:MM:SS")
        end
        sort!(method_file, "Date Acquired")
        method_file = insertcols(method_file[:, 2:end], 1, "Name" => "NaN")                                 # Changed the Name column type from Missing to String
        method_file = insertcols(method_file, 2, "InChIKey" => "NaN")                                      # Adding InChIKey column
        method_file = insertcols(method_file, 6, "Method" => methods_names[method_no])                     # Adding Method column
        # Check for missing values and removing the peaks with missing information
        while any(ismissing.(Matrix(method_file[:, 1:end-2])))
            missing_idx = findall(x -> ismissing.(x), Matrix(method_file[:, 1:end-2]))[1][1]
            deleteat!(method_file, missing_idx)
        end
        method_file = method_file[method_file[:, "Area"].>=min_area_count, :]                             # Filter on a minimum Area count of 100,000
        uniq_vials = unique(method_file[:, "Vial"])

        for pest_mix = 1:length(uniq_vials)                                                                                   # Injection 1-10, Pest mix 1-10
            sample_file = method_file[findall(x -> x .== uniq_vials[pest_mix], method_file[:, "Vial"]), :]
            sample_file = insertcols(sample_file, 6, "Mix" => pest_mix)
            ref_file = CSV.read(pest_dir * "pest_mix_$(pest_mix).csv", DataFrame)
            ref_file = sort(insertcols(ref_file, 5, "H_MonoisotopicMass" => (ref_file[:, "MonoisotopicMass"] .+ H_adduct_mass)), "H_MonoisotopicMass", rev=true)
            if names(ref_file)[1] == "Chemical Name"                                                           # If true, "Chemical Name" should be changed to "Name"
                ref_file = insertcols(ref_file[:, 2:end], 1, "Name" => ref_file[:, "Chemical Name"])
            end

            # We associate the Masses from the sample to the reference file and put the Name and INCHIKEY
            for peak_no = 1:size(sample_file, 1)
                peak_mass_string = strip(string(sample_file[peak_no, mass_col]))
                if contains(peak_mass_string, '(')      # In some cases the Assigned mass contains extra info that must be removed e.g. "81.5 (2)"
                    peak_mass = parse(Float64, peak_mass_string[1:(findfirst(x -> x .== '(', peak_mass_string)-2)])
                else
                    peak_mass = parse(Float64, peak_mass_string[1:(length(peak_mass_string))])
                end
                difference = abs.(peak_mass .- ref_file[:, "H_MonoisotopicMass"])
                if any(difference .<= mass_tolerance)                                                     # If true, we have found a match
                    if length(findall(x -> x .== findmin(difference)[1], difference)) == 1                # If true, we have exactly one match
                        matched_comp = ref_file[argmin(difference), :]
                        sample_file[peak_no, "Name"] = matched_comp["Name"]
                        sample_file[peak_no, "InChIKey"] = matched_comp["InChIKey"]
                    elseif length(findall(x -> x .== findmin(difference)[1], difference)) > 1             # If true, we have several matches (e.g. isomer)
                        matched_comps = ref_file[findall(x -> x .== findmin(difference)[1], difference), :]
                        name = matched_comps[1, "Name"]
                        inchik = matched_comps[1, "InChIKey"]
                        for matched_compounds_i = 2:size(matched_comps, 1)
                            name = name * "/" * matched_comps[matched_compounds_i, "Name"]
                            inchik = inchik * "/" * matched_comps[matched_compounds_i, "InChIKey"]
                        end
                        sample_file[peak_no, "Name"] = name
                        sample_file[peak_no, "InChIKey"] = inchik
                    end
                end
            end
            all_methods = append!(all_methods, sample_file)
        end
    end
    println("Identification successful for all $(length(method_files)) methods")
    all_methods[all_methods[:, "Name"].!="NaN", "Name"] .= lowercase.(all_methods[all_methods[:, "Name"].!="NaN", "Name"])
    dropnan_df = all_methods[findall(x -> x != "NaN", all_methods[:, "Name"]), :]
    return (all_methods, dropnan_df)
end
ident_df, dropnan_df = identification() # Unidentified peaks are removed to dropnan_df

## Retention times
## For each method and mix combination, make a dataframe with the relevant data and make an extra column with whether it's present or not.
function set_RT(dropnan_df=dropnan_df, pest_dir="C://Users//alex_//Documents//Github//Method_ranking//data//AZ_raw//Pesticide mix info//", rt_tolerance=0.2)
    # Prepare a dataframe to insert the values Present? and RT Window
    reference_list = DataFrame()
    for ref_file_no = 1:(sum(contains.(readdir(pest_dir), "pest_mix")))
        mix_to_add = CSV.read(pest_dir * "pest_mix_$(ref_file_no).csv", DataFrame)
        mix_to_add[!, 1] = convert.(String, mix_to_add[:, 1])
        mix_to_add[!, 2] = convert.(String, mix_to_add[:, 2])
        if names(mix_to_add)[1] == "Chemical Name"                                                           # If true, "Chemical Name" should be change to "Name"
            mix_to_add = insertcols(mix_to_add[:, 2:end], 1, "Name" => mix_to_add[:, "Chemical Name"])
        end
        insertcols!(mix_to_add, 1, "Reference Mix" => ref_file_no)
        reference_list = append!(reference_list, mix_to_add, promote=true)
    end
    reference_list[!, "Name"] = lowercase.(reference_list[!, "Name"])
    reference_list = insertcols(reference_list, 3, "Matched Peaks" => -1)
    reference_list = insertcols(reference_list, 4, "RT Window" => "?")
    reference_list = insertcols(reference_list, 4, "RT" => -1.0)

    # For each method
    Method_names = unique(dropnan_df[:, "Method"])
    ref_method_all = DataFrame()
    for method_no in ProgressBar(1:length(unique(dropnan_df[:, "Method"])))
        method_i = dropnan_df[dropnan_df[:, "Method"].==Method_names[method_no], :]
        ref_method_i = deepcopy(reference_list)
        ref_method_i = insertcols(ref_method_i, 1, "Method Name" => Method_names[method_no])
        for ref_comp_idx = 1:length(ref_method_i[:, "Name"])
            comp_name = ref_method_i[ref_comp_idx, "Name"]
            comp_idx = findall(x -> occursin.(comp_name, x), method_i[:, "Name"])
            comp_subdf = method_i[comp_idx, :]
            matches_no = size(comp_subdf, 1)
            if (matches_no != 0) && (length(unique(comp_subdf[:, "Mix"])) != 1) && (unique(comp_subdf[:, "Mix"]) != ref_method_i[ref_comp_idx, "Reference Mix"])
                error("There is a match with the compound '$comp_name' that exists in more than one mix. This error should never occur")
            end
            ref_method_i[ref_comp_idx, "Matched Peaks"] = matches_no

            # These are the compounds that are shown zero times (i.e. 0 peaks corresponds to 1 compound)
            if matches_no == 0
                ref_method_i[ref_comp_idx, "RT"] = 0.0
                ref_method_i[ref_comp_idx, "RT Window"] = "NaN"
                # These are the compounds that are shown once (i.e. 1 peak corresponds to 1 compound)
            elseif matches_no == 1
                ref_method_i[ref_comp_idx, "RT"] = comp_subdf[:, "RT"][1]
                # These are the compounds that are shown twice (i.e. 2 peaks correspond to 1 compound) 
                # (1) If the two peaks are close to each other
                # (2) If the two peaks are far from each other
            elseif matches_no == 2
                # (1)      
                if abs(maximum(comp_subdf[:, "RT"]) - minimum(comp_subdf[:, "RT"])) <= rt_tolerance
                    ref_method_i[ref_comp_idx, "RT"] = round(mean(comp_subdf[:, "RT"]), digits=5)
                    # (2)
                else
                    ref_method_i[ref_comp_idx, "RT"] = comp_subdf[:, "RT"][argmax(comp_subdf[:, "Area"])]
                end
                # These are the compounds that are shown more than twice (i.e. >2 peaks correspond to 1 compound)
                # (1) If all the peaks are close to each other. All peaks should be kept
                # (2) At least one of the peaks must be removed. Should be in a while true break loop to remove all outliers. 
            elseif matches_no > 2
                # (1) 
                if abs(maximum(comp_subdf[:, "RT"]) - minimum(comp_subdf[:, "RT"])) <= rt_tolerance
                    ref_method_i[ref_comp_idx, "RT"] = round(mean(comp_subdf[:, "RT"]), digits=5)
                    # (2)
                else
                    set = comp_subdf[:, "RT"]
                    while true
                        set_normalized = set .- mean(set)
                        distances_sum = [sum(abs.(set_normalized_i .- set_normalized)) for set_normalized_i in set_normalized]
                        # These are the compounds that are shown twice (i.e. 2 peaks correspond to 1 compound) 
                        # (1) If the two peaks are close to each other
                        # (2) If the two peaks are far from each other
                        if length(distances_sum) == 2
                            set = comp_subdf[findall(x -> x in set, comp_subdf[:, "RT"]), :]
                            # (1)
                            if abs(maximum(set[:, "RT"]) - minimum(set[:, "RT"])) <= rt_tolerance
                                ref_method_i[ref_comp_idx, "RT"] = round(mean(set[:, "RT"]), digits=5)
                                # (2)
                            else
                                ref_method_i[ref_comp_idx, "RT"] = set[argmax(comp_subdf[:, "Area"]), "RT"]
                            end
                            break
                            # If more than 2 peaks and they are close to each other, take the average.
                        elseif abs(maximum(set) - minimum(set)) <= rt_tolerance
                            ref_method_i[ref_comp_idx, "RT"] = round(mean(set), digits=5)
                            break
                            # Otherwise, remove the peak with the most RT distance to the rest and repeat
                        else
                            set = set[Not(argmax(distances_sum))]
                        end
                    end
                end
            else
                error("Error with $matches_no number of matches for compound $comp_name. This error should never occur")
            end
        end
        append!(ref_method_all, ref_method_i)
    end
    println("Retention times completed")
    return ref_method_all
end
training_dataset = set_RT();

## Making RT windows
function set_RT_windows(training_dataset=training_dataset, windows=5, normalise_to_first_peak=false)
    dataset = deepcopy(training_dataset)
    # Normalising the retention times
    insertcols!(dataset, "RT Window", "RT_norm" => -1.0)
    min_RT_method_i = [minimum(ident_df[ident_df[:, "Method"].==(i), "RT"]) for i in unique(ident_df[:, "Method"])]
    max_RT_method_i = [maximum(ident_df[ident_df[:, "Method"].==(i), "RT"]) for i in unique(ident_df[:, "Method"])]
    for method_no in ProgressBar(1:length(unique(dataset[:, "Method Name"])))
        method_name_i = unique(dataset[:, "Method Name"])[method_no]
        method_i_idx = findall(x -> x .== method_name_i, dataset[:, "Method Name"])
        if normalise_to_first_peak
            RT_min = min_RT_method_i[method_no]
        elseif normalise_to_first_peak == false
            RT_min = 0
        else
            error("Set the normalization parameter. Should the normalization be done on the first peak?")
        end
        RT_max = max_RT_method_i[method_no]
        dataset[method_i_idx, "RT_norm"] = (dataset[method_i_idx, "RT"] .- RT_min) / (RT_max - RT_min)
    end

    # Checking for errors
    NaN_idx = findall(x -> x .== "NaN", dataset[:, "RT Window"])
    if any(dataset[NaN_idx, "RT_norm"] .> 0)
        error("A compound that appeared as a peak is assigned to the Not_present window. This error should never occur.")
    end
    if length(findall(x -> x > 1, dataset[:, "RT_norm"])) != 0
        error("A compound has a normalised RT higher than 1. This error should never occur.")
    end

    # Defining RT windows
    if normalise_to_first_peak == true || windows != 5
        error("To be implemented. The windows need to change for that setting.")
    end

    dataset[findall(x -> 0 < x <= 0.2, dataset[:, "RT_norm"]), "RT Window"] .= "A"
    dataset[findall(x -> 0.2 < x <= 0.4, dataset[:, "RT_norm"]), "RT Window"] .= "B"
    dataset[findall(x -> 0.4 < x <= 0.6, dataset[:, "RT_norm"]), "RT Window"] .= "C"
    dataset[findall(x -> 0.6 < x <= 0.8, dataset[:, "RT_norm"]), "RT Window"] .= "D"
    dataset[findall(x -> 0.8 < x <= 1.0, dataset[:, "RT_norm"]), "RT Window"] .= "E"
    println("RT windows completed")
    return dataset
end
dataset = set_RT_windows();

## Scoring
## The ideal score for three compounds would be for them to be in windows B,C,D.
## So we have positive score for homogeneity in the 5 windows. 
## We have negative score for compounds that are not present (NaN).
## We have medium negative score for window A and small negative score for window E.
function set_score(dataset; w_homo=100, w_first=-5, w_last=-2, w_NaN=-20, windows=5, sorted=true)
    scores_list = DataFrame()
    no_of_methods = length(unique(dataset[:, "Method Name"]))

    for method_i in ProgressBar(1:no_of_methods)    # method_i=2
        method_name_i = unique(dataset[:, "Method Name"])[method_i]
        method_i_idx = findall(x -> x .== method_name_i, dataset[:, "Method Name"])
        method_set = dataset[method_i_idx, :]
        RT_windows = dataset[method_i_idx, "RT Window"]

        no_of_theoretical_comps = length(RT_windows)
        no_of_present_comps = length(RT_windows[RT_windows.!="NaN"])
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
        z_first = w_first * ratios_df[ratios_df.Window.==ratios_df.Window[1], "Ratio"][1]
        z_last = w_last * ratios_df[ratios_df.Window.==ratios_df.Window[end], "Ratio"][1]
        z_NaN = w_NaN * ((no_of_theoretical_comps - no_of_present_comps) / no_of_theoretical_comps)

        # Score
        score_val = z_homo + z_first + z_last + z_NaN
        scores_list_to_add = DataFrame("Method number" => method_i, "Method Name" => method_name_i, "Score" => score_val)
        scores_list = append!(scores_list, scores_list_to_add)
    end
    scores_list_sorted = sort(scores_list, "Score", rev=true)

    println("The optimal method for the selected compounds is $(scores_list_sorted[1,"Method Name"])")
    println("The 2nd optimal method for the selected compounds is $(scores_list_sorted[2,"Method Name"])")
    println("The 3rd optimal method for the selected compounds is $(scores_list_sorted[3,"Method Name"])")
    if sorted
        return scores_list_sorted
    else
        return scores_list
    end
end
scores = set_score(dataset);
#= Exploratory analysis
## PROBLEM with duplicate measurement files.
duplicate_score_vals = names(sort(freqtable(scores_list))[end-1:end])[1]
duplicate_idx = findall(x -> x in duplicate_score_vals, scores_list)
duplicates = hcat(duplicate_idx, scores_list[duplicate_idx], unique(dataset[:, "Method Name"])[duplicate_idx])

# plots
windows = ["A","B","C","D","E"]
comps = [2,0,2,15,1]
comps = [3,5,4,5,3]
sum(comps)
bar(windows,comps,xlabel="RT window", ylabel="Compounds", legend=false, c=:darkgrey,ylims=(0,15))
=#

# Model
using ScikitLearn, PyCall, Conda
using ScikitLearn.CrossValidation: train_test_split
#using ScikitLearn.GridSearch: RandomizedSearchCV

cat = pyimport("catboost")
cat_utils = pyimport("catboost.utils")
plt = pyimport("matplotlib")
jblb = pyimport("joblib")
pcp = pyimport("pubchempy")
pd = pyimport("padelpy")
#ipywidgets = pyimport("ipywidgets")
#IPython = pyimport("IPython")

############################################################################################
# Creating fingerprints
function from_CID_to_fp_fast(dataset::DataFrame; dictionary_file_dir="C:\\Users\\alex_\\Documents\\GitHub\\Method_ranking\\data\\Fingerprints\\full_fingerprints_dictionary.csv", create_dict=false, pubchem_compression=true)
    ## PubChem compression FPs
    function compressPubChemFPs(ACfp::DataFrame, PCfp::DataFrame)
        FP1tr = convert.(Float64, ACfp)
        pubinfo = convert.(Int, Matrix(PCfp))
        findidx(FP_number) = findfirst(x -> x .== "PubchemFP$FP_number", names(PCfp))

        #ring counts
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
        println("Compressed fingerprints calculated")

        return custom_FPs_df
    end

    rep = deepcopy(dataset)
    if isempty(dictionary_file_dir) || create_dict
        unique_dataset_CIDs = unique(rep[:, "PubChem CID"])

        shortest_cid = pcp.get_compounds(unique_dataset_CIDs[1], "cid")[1]
        fp_temp = DataFrame(pd.from_smiles(shortest_cid.isomeric_smiles, fingerprints=true, descriptors=false))
        fp_temp = convert.(Int16, parse.(Float16, fp_temp))
        fp_padel_dict = hcat(DataFrame("Name" => unique(rep[:, "Name"])[1], "InChIKey" => unique(rep[:, "InChIKey"])[1], "CID" => shortest_cid.cid), fp_temp)

        for i in ProgressBar(2:length(unique_dataset_CIDs))
            try
                shortest_cid = pcp.get_compounds(unique_dataset_CIDs[i], "cid")[1]
                fp_temp = DataFrame(pd.from_smiles(shortest_cid.isomeric_smiles, fingerprints=true, descriptors=false))
                fp_temp = convert.(Int16, parse.(Float16, fp_temp))
                fp_padel_dict_to_add = hcat(DataFrame("Name" => unique(rep[:, "Name"])[i], "InChIKey" => unique(rep[:, "InChIKey"])[i], "CID" => shortest_cid.cid), fp_temp)
                fp_padel_dict = append!(fp_padel_dict, fp_padel_dict_to_add)
            catch
                println("Problem with compound $i. Continuing...")
                continue
            end
        end
        # PubChem Compression
        if pubchem_compression
            ACfp = fp_padel_dict[:, findall(x -> occursin.("APC2D", x), names(fp_padel_dict))]
            PCfp = fp_padel_dict[:, findall(x -> occursin.("Pubchem", x), names(fp_padel_dict))]
            fp_padel_dict = hcat(fp_padel_dict[:, 1:3], compressPubChemFPs(ACfp, PCfp))
        end

        CSV.write(dictionary_file_dir, fp_padel_dict)
        return fp_padel_dict
    else
        fp_padel_dict = CSV.read(dictionary_file_dir, DataFrame)
        fp_padel_dict[:, 3:end] = Int.(fp_padel_dict[:, 3:end])
    end

    ## Associate the rep to the fp_padel_dict
    rep_FP = DataFrame()
    for compound_i in ProgressBar(1:size(rep, 1))   # Delete after debugging: compound_i = 2015
        fp_padel_dict_comp_idx = findfirst(x -> x .== rep[compound_i, "Name"], fp_padel_dict[:, "Name"])
        fp_padel_dict_comp = fp_padel_dict[fp_padel_dict_comp_idx, :]
        # Checking if the provided PubChem CID is the same as the calculated CID. This part can be removed later!
        if fp_padel_dict_comp["CID"] != rep[compound_i, "PubChem CID"]
            println("There is a CID mismatch at compound $compound_i")
        end
        #
        fp_padel_dict_comp
        rep_FP_to_add = hcat(DataFrame(rep[compound_i, :]), DataFrame(fp_padel_dict_comp[4:end]))
        rep_FP = append!(rep_FP, rep_FP_to_add)
    end
    println("FP calculation complete")
    return rep_FP
end
fp_dataset = from_CID_to_fp_fast(dataset)

############################################################################################
## Modeling
function confusion_matrix_df(y_true, y_pred)
    classes = sort(unique(vcat(y_true, y_pred)))  # Determine classes and sort them alphabetically
    class_to_index = Dict(c => i for (i, c) in enumerate(classes))
    conf_matrix = zeros(Int, length(classes), length(classes))

    for (t, p) in zip(y_true, y_pred)
        conf_matrix[class_to_index[t], class_to_index[p]] += 1
    end

    df = hcat(DataFrame(Class=classes), DataFrame(conf_matrix, classes))  # Add class names as the first column
    return df
end

## Presence model => IS THE COMP PRESENT?
fp_dataset_presence = deepcopy(fp_dataset)
fp_dataset_presence[findall(x -> x .!= "NaN", fp_dataset_presence[:, "RT Window"]), "RT Window"] .= "Yes"
seed = 1
# Presence model - Train test split
last_info_col_idx = findfirst(x -> x .== "SMILES", names(fp_dataset_presence)) # SMILES should be the last informational column.
train_inchikeys, test_inchikeys = train_test_split(unique(fp_dataset_presence[:, "InChIKey"]), test_size=0.2, random_state=seed)
train_idx = findall(x -> x in train_inchikeys, fp_dataset_presence[:, "InChIKey"])
test_idx = findall(x -> x in test_inchikeys, fp_dataset_presence[:, "InChIKey"])

X_train_df = hcat(fp_dataset_presence[train_idx, 1], fp_dataset_presence[train_idx, (last_info_col_idx+1):end])
X_train = Matrix(X_train_df)
X_test_df = hcat(fp_dataset_presence[test_idx, 1], fp_dataset_presence[test_idx, (last_info_col_idx+1):end])
X_test = Matrix(X_test_df)
Y_train = fp_dataset_presence[train_idx, "RT Window"]
Y_test = fp_dataset_presence[test_idx, "RT Window"]
train_pool = cat.Pool(X_train, label=Y_train, cat_features=[0])
test_pool = cat.Pool(X_test, label=Y_test, cat_features=[0])

# Presence model - Optimization v0.1
function optimise_presence(n_iter; csv_write::Bool=false, seed=2)
    ## Presence model => IS THE COMP PRESENT?
    fp_dataset_presence = deepcopy(fp_dataset)
    fp_dataset_presence[findall(x -> x .!= "NaN", fp_dataset_presence[:, "RT Window"]), "RT Window"] .= "Yes"

    # Presence model - Train test split
    last_info_col_idx = findfirst(x -> x .== "SMILES", names(fp_dataset_presence)) # SMILES should be the last informational column.
    train_inchikeys, test_inchikeys = train_test_split(unique(fp_dataset_presence[:, "InChIKey"]), test_size=0.2, random_state=seed)
    train_idx = findall(x -> x in train_inchikeys, fp_dataset_presence[:, "InChIKey"])
    test_idx = findall(x -> x in test_inchikeys, fp_dataset_presence[:, "InChIKey"])

    X_train_df = hcat(fp_dataset_presence[train_idx, 1], fp_dataset_presence[train_idx, (last_info_col_idx+1):end])
    X_train = Matrix(X_train_df)
    X_test_df = hcat(fp_dataset_presence[test_idx, 1], fp_dataset_presence[test_idx, (last_info_col_idx+1):end])
    X_test = Matrix(X_test_df)
    Y_train = fp_dataset_presence[train_idx, "RT Window"]
    Y_test = fp_dataset_presence[test_idx, "RT Window"]
    train_pool = cat.Pool(X_train, label=Y_train, cat_features=[0])
    test_pool = cat.Pool(X_test, label=Y_test, cat_features=[0])

    grid = Dict("learning_rate" => collect(0.02:0.0025:0.05), "iterations" => collect(1500:250:2000), "depth" => collect(6:1:10), "l2_leaf_reg" => collect(1:2:9), "random_seed" => collect(1:3))
    results_df = DataFrame()
    presence_model = cat.CatBoostClassifier(early_stopping_rounds=50, thread_count=-1, verbose=false)
    grid_search_results = presence_model.randomized_search(grid, train_pool, n_iter=n_iter, cv=3, search_by_train_test_split=false, shuffle=true, verbose=false)

    # Reporting
    results = hcat(DataFrame(grid_search_results["params"]), rename(DataFrame(DataFrame(grid_search_results["cv_results"])[end, 1:2]), "iterations" => "actual iterations"))
    results = hcat(results, DataFrame("Train accuracy" => score(presence_model, X_train, Y_train), "Test accuracy" => score(presence_model, X_test, Y_test)))
    if csv_write
        CSV.write("params_presence_model_v0.1.csv", results, append=true)  #end    
    end
    results_df = append!(results_df, results)
    return results_df, presence_model
end
#results_df, presence_model = optimise_presence(400, csv_write=true)    # 4-Sep
results_df_1, presence_model_1 = optimise_presence(150, seed=1 , csv_write=false)
CSV.write("20240905_Optimization_presence_seed1.csv",results_df_1)
results_df_2, presence_model_2 = optimise_presence(150, seed=2 , csv_write=false)
CSV.write("20240905_Optimization_presence_seed2.csv",results_df_2)
results_df_3, presence_model_3 = optimise_presence(150, seed=3 , csv_write=false)
CSV.write("20240905_Optimization_presence_seed3.csv",results_df_3)

presence_model = presence_model_1
# Loading model 
presence_model = cat.CatBoostClassifier()
presence_model.load_model(project_path*"models//"*"presence_model_v0.2.bin")
c_train = confusion_matrix_df(Y_train, convert.(String,presence_model.predict(X_train)))
acc_train = (c_train[1,2] + c_train[2,3])/ sum(Matrix(c_train[:,2:3]))
c_test = confusion_matrix_df(Y_test, convert.(String,presence_model.predict(X_test)))
acc_test = (c_test[1,2] + c_test[2,3])/ sum(Matrix(c_test[:,2:3]))

# Presence model - Changing prediction threshold
using Plots
fpr, tpr, thresholds = cat_utils.get_roc_curve(presence_model, train_pool)
fpr_test, tpr_test, thresholds = cat_utils.get_roc_curve(presence_model, test_pool)
plot(fpr, tpr, xlabel="FPR", ylabel="TPR", label="Train")
plot!(fpr_test, tpr_test, xlabel="FPR", ylabel="TPR", label="Test")
threshold = cat_utils.select_threshold(presence_model, curve=cat_utils.get_roc_curve(presence_model, train_pool))
cat_utils.select_threshold(presence_model, curve=cat_utils.get_roc_curve(presence_model, test_pool))
presence_model.set_probability_threshold(binclass_probability_threshold=threshold)

z2 = score(presence_model, X_train, Y_train)                      # Train set accuracy
z3 = score(presence_model, X_test, Y_test)                        # Test set accuracy
z4 = presence_model.predict_proba(X_train)                        # y_hat_train
z4 = presence_model.predict(X_train)                              # y_hat_train
z5 = presence_model.predict(X_test)                               # y_hat_test

# Presence model - Save model (optimized hyperparameters + optimized prediction threshold)
if presence_model.is_fitted()
    presence_model.save_model("presence_model_v0.2.bin") 
end

# Presence model - Plots
thresholds_train, fpr_train = cat_utils.get_fpr_curve(presence_model, train_pool)
thresholds_train, fnr_train = cat_utils.get_fnr_curve(presence_model, train_pool)
thresholds_test, fpr_test = cat_utils.get_fpr_curve(presence_model, test_pool)
thresholds_test, fnr_test = cat_utils.get_fnr_curve(presence_model, test_pool)

plot(thresholds_train, fpr_train, xlabel="Threshold", label="Train FPR", legend=:right)
plot!(thresholds_train, fnr_train, label="Train FNR")
plot!(thresholds_train, fnr_train + fpr_train, label="Sum(FPR,FNR)")
plot!([threshold, threshold], [0, 1], label="Chosen Threshold")

plot(thresholds_test, fpr_test, xlabel="Threshold", label="Test FPR", legend=:right)
plot!(thresholds_test, fnr_test, label="Test FNR")
plot!(thresholds_test, fnr_test + fpr_test, label="Sum(FPR,FNR)")
plot!([threshold, threshold], [0, 1], label="Chosen Threshold")

############################################################################################
## Window model => WHAT IS THE WINDOW?
fp_dataset_window = deepcopy(fp_dataset[findall(x -> x .!= "NaN", fp_dataset[:, "RT Window"]), :])
seed=2
#method_names_list = DataFrame("Method Name"=>unique(fp_dataset[:,"Method Name"]))
#CSV.write("method_names.csv", method_names_list)
# Train test split
last_info_col_idx = findfirst(x -> x .== "SMILES", names(fp_dataset_window)) # SMILES should be the last informational column.
train_inchikeys, test_inchikeys = train_test_split(unique(fp_dataset_window[:, "InChIKey"]), test_size=0.2, random_state=seed)
train_idx = findall(x -> x in train_inchikeys, fp_dataset_window[:, "InChIKey"])
test_idx = findall(x -> x in test_inchikeys, fp_dataset_window[:, "InChIKey"])

X_train_df = hcat(fp_dataset_window[train_idx, 1], fp_dataset_window[train_idx, (last_info_col_idx+1):end])
X_train = Matrix(X_train_df)
X_test_df = hcat(fp_dataset_window[test_idx, 1], fp_dataset_window[test_idx, (last_info_col_idx+1):end])
X_test = Matrix(X_test_df)
Y_train = fp_dataset_window[train_idx, "RT Window"]
Y_test = fp_dataset_window[test_idx, "RT Window"]

# Window model - Optimization v0.1
function optimise_window_v01(n_iter; csv_write::Bool=false, seed=2)
    ## Window model => WHAT IS THE WINDOW?
    fp_dataset_window = deepcopy(fp_dataset[findall(x -> x .!= "NaN", fp_dataset[:, "RT Window"]), :])

    # Train test split
    last_info_col_idx = findfirst(x -> x .== "SMILES", names(fp_dataset_window)) # SMILES should be the last informational column.
    train_inchikeys, test_inchikeys = train_test_split(unique(fp_dataset_window[:, "InChIKey"]), test_size=0.2, random_state=seed)
    train_idx = findall(x -> x in train_inchikeys, fp_dataset_window[:, "InChIKey"])
    test_idx = findall(x -> x in test_inchikeys, fp_dataset_window[:, "InChIKey"])

    X_train_df = hcat(fp_dataset_window[train_idx, 1], fp_dataset_window[train_idx, (last_info_col_idx+1):end])
    X_train = Matrix(X_train_df)
    X_test_df = hcat(fp_dataset_window[test_idx, 1], fp_dataset_window[test_idx, (last_info_col_idx+1):end])
    X_test = Matrix(X_test_df)
    Y_train = fp_dataset_window[train_idx, "RT Window"]
    Y_test = fp_dataset_window[test_idx, "RT Window"]
    train_pool = cat.Pool(X_train, label=Y_train, cat_features=[0])
    test_pool = cat.Pool(X_test, label=Y_test, cat_features=[0])

    grid = Dict("learning_rate" => collect(0.02:0.0025:0.05), "iterations" => collect(1500:250:2000), "depth" => collect(6:1:10), "l2_leaf_reg" => collect(1:2:9), "random_seed" => collect(1:3))
    results_df = DataFrame()
    window_model = cat.CatBoostClassifier(early_stopping_rounds=50, thread_count=-1, verbose=false)
    grid_search_results = window_model.randomized_search(grid, train_pool, n_iter=n_iter, cv=3, search_by_train_test_split=false, shuffle=true, verbose=false)

    # Reporting
    results = hcat(DataFrame(grid_search_results["params"]), rename(DataFrame(DataFrame(grid_search_results["cv_results"])[end, 1:2]), "iterations" => "actual iterations"))
    results = hcat(results, DataFrame("Train accuracy" => score(window_model, X_train, Y_train), "Test accuracy" => score(window_model, X_test, Y_test)))
    if csv_write
        file_name = "params_window_model_v0.1.csv"
        if any(readdir(project_path*"models//") .== file_name)
        CSV.write(project_path*"models//"*file_name, results, append=false) 
        else
        CSV.write("params_window_model_v0.1.csv", results, append=true)  end    
    end
    results_df = append!(results_df, results)
    return results_df, window_model
end
# Window model - Optimization v0.2
function optimise_window_v02(n_iter; csv_write::Bool=false, seed=2)
    ## Window model => WHAT IS THE WINDOW?
    fp_dataset_window = deepcopy(fp_dataset[findall(x -> x .!= "NaN", fp_dataset[:, "RT Window"]), :])

    # Train test split
    last_info_col_idx = findfirst(x -> x .== "SMILES", names(fp_dataset_window)) # SMILES should be the last informational column.
    train_inchikeys, test_inchikeys = train_test_split(unique(fp_dataset_window[:, "InChIKey"]), test_size=0.2, random_state=seed)
    train_idx = findall(x -> x in train_inchikeys, fp_dataset_window[:, "InChIKey"])
    test_idx = findall(x -> x in test_inchikeys, fp_dataset_window[:, "InChIKey"])

    X_train_df = hcat(fp_dataset_window[train_idx, 1], fp_dataset_window[train_idx, (last_info_col_idx+1):end])
    X_train = Matrix(X_train_df)
    X_test_df = hcat(fp_dataset_window[test_idx, 1], fp_dataset_window[test_idx, (last_info_col_idx+1):end])
    X_test = Matrix(X_test_df)
    Y_train = fp_dataset_window[train_idx, "RT Window"]
    Y_test = fp_dataset_window[test_idx, "RT Window"]
    train_pool = cat.Pool(X_train, label=Y_train, cat_features=[0])
    test_pool = cat.Pool(X_test, label=Y_test, cat_features=[0])

    grid = Dict("learning_rate" => collect(0.02:0.0025:0.05), "iterations" => collect(1500:250:2000), "depth" => collect(6:1:10), "l2_leaf_reg" => collect(1:2:9), "random_seed" => collect(1:3), "loss_function"=>"AUC")
    results_df = DataFrame()
    window_model = cat.CatBoostClassifier(early_stopping_rounds=50, thread_count=-1, verbose=false)
    grid_search_results = window_model.randomized_search(grid, train_pool, n_iter=n_iter, cv=3, search_by_train_test_split=false, shuffle=true, verbose=false)
    cv_search_results = cat.cv(pool=train_pool, params=grid, shuffle=true, verbose=false)

    # Reporting
    results = hcat(DataFrame(grid_search_results["params"]), rename(DataFrame(DataFrame(grid_search_results["cv_results"])[end, 1:2]), "iterations" => "actual iterations"))
    results = hcat(results, DataFrame("Train accuracy" => score(window_model, X_train, Y_train), "Test accuracy" => score(window_model, X_test, Y_test)))
    if csv_write
        file_name = "params_window_model_v0.1.csv"
        if any(readdir(project_path*"models//") .== file_name)
        CSV.write(project_path*"models//"*file_name, results, append=false) 
        else
        CSV.write("params_window_model_v0.1.csv", results, append=true)  end    
    end
    results_df = append!(results_df, results)
    return results_df, window_model
end

#results_df, presence_model = optimise_presence(400, csv_write=true)    # 4-Sep
results_df_1, window_model_1 = optimise_window_v01(150, seed=1 , csv_write=false)
CSV.write(project_path*"models//"*"20240917_Optimization_window_seed1.csv",results_df_1)
results_df_2, window_model_2 = optimise_window_v01(150, seed=2 , csv_write=false)
CSV.write(project_path*"models//"*"20240917_Optimization_window_seed2.csv",results_df_2)
results_df_3, window_model_3 = optimise_window_v01(225, seed=3 , csv_write=false)
CSV.write(project_path*"models//"*"20240917_Optimization_window_seed3.csv",results_df_3)

window_model = cat.CatBoostClassifier(early_stopping_rounds=50, thread_count=-1, verbose=100, depth=7, iterations=2000, l2_leaf_reg=1, learning_rate=0.0325, random_state=3)
window_model.fit(X_train, Y_train, cat_features=[0], eval_set=(X_test, Y_test))

window_model = cat.CatBoostClassifier(early_stopping_rounds=50, thread_count=-1, verbose=100, depth=7, iterations=2000, l2_leaf_reg=1, learning_rate=0.0325, random_state=3)
window_model.fit(X_train, Y_train, cat_features=[0])

z2 = score(window_model, X_train, Y_train)                      # Train set accuracy
z3 = score(window_model, X_test, Y_test)                        # Test set accuracy
z4 = window_model.predict_proba(X_train)                        # y_hat_train
z4 = window_model.predict(X_train)                              # y_hat_train
z5 = window_model.predict(X_test)                               # y_hat_test

# Optimization
train_pool = cat.Pool(X_train, label=Y_train, cat_features=[0])
test_pool = cat.Pool(X_test, label=Y_test, cat_features=[0])
window_model = cat.CatBoostClassifier(early_stopping_rounds=50, thread_count=-1, verbose=100)
grid = Dict("learning_rate" => [0.01, 0.02, 0.03], "iterations" => [700, 1000, 1300], "depth" => [4, 6, 8, 10], "l2_leaf_reg" => [1, 3, 5, 7, 9], "random_seed" => [1, 2, 3, 4, 5])
grid_fast = Dict("learning_rate" => [0.02], "random_seed" => [1, 2])
grid_search_results_window = window_model.grid_search(grid_fast, train_pool, cv=3, shuffle=true, verbose=100, plot_file="window_model_params.html")
opt_params_window = DataFrame(grid_search_results_window["params"])
cv_results_window = DataFrame(grid_search_results_window["cv_results"])
if window_model.is_fitted() == false
    error("Presence model has not been fitted.")
end
window_model.get_params()

# Save model (optimized hyperparameters)
window_model.save_model(project_path*"models//"*"window_model_v0.1.bin")

# Plots
# Load model (optimized hyperparameters)
window_model = cat.CatBoostClassifier()
window_model.load_model(project_path*"models//"*"window_model_v0.1.bin")

thresholds_train, fpr_train = cat_utils.get_fpr_curve(window_model, train_pool)
thresholds_train, fnr_train = cat_utils.get_fnr_curve(window_model, train_pool)
thresholds_test, fpr_test = cat_utils.get_fpr_curve(window_model, test_pool)
thresholds_test, fnr_test = cat_utils.get_fnr_curve(window_model, test_pool)
plot(thresholds_train, fpr_train, xlabel="Threshold", label="Train FPR", legend=:right)
plot!(thresholds_train, fnr_train, label="Train FNR")
plot!(thresholds_train, fnr_train + fpr_train, label="Sum(FPR,FNR)")
plot!([threshold, threshold], [0, 1], label="Chosen Threshold")

# Confusion Matrix
using LinearAlgebra
c_train = confusion_matrix_df(Y_train, convert.(String, window_model.predict(X_train)[:]))
acc_train = (sum(diag(Matrix(c_train[:,2:end]))))/ sum(Matrix(c_train[:,2:end]))
c_test = confusion_matrix_df(Y_test, convert.(String,window_model.predict(X_test)))
acc_test = (sum(diag(Matrix(c_test[:,2:end]))))/ sum(Matrix(c_test[:,2:end]))

using Plots
heatmap(Matrix(c_test)[:, 2:end])
