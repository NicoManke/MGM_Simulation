#Macrophyte Growth Model (MGM)
#
#Anne Lewerentz <anne.lewerentz@uni-wuerzburg.de>
#(c) 2021-2022, licensed under the terms of the MIT license
#
#Contains all different output functions that can be used to generate output from MGM
#e.g. from within R
#
#Functions:
#* CHARISMA_biomass() : Returns: Mean summer biomass for all lakes, species, and multiple depths
#* CHARISMA_biomass_parallel() : same output, but runs parallalised
#* CHARISMA_biomass_parallel_lastNyears(): same output, but for multiple years
#* CHARISMA_biomass_onedepth() : same output, but just for one depth
#
#Functions that where used during development of the code to test connection with R
#* CHARISMA_test_15lakes_4depths()
#* CHARISMA_parallel_test_15lakes_4depths()
#
#Out-dated functions
#* CHARISMA_VE()


#using
#    HCubature, #for Integration
#    DelimitedFiles, # for function writedlm, used to write output files
#    Dates, #to create output folder
#    Distributions, Random, #for killWithProbability
#    CSV, #For virtual Ecologist
#    DataFrames, #For virtual Ecologist
#    StatsBase #For virtual Ecologist

# Include functions
#include("structs.jl")
#include("defaults.jl")
#include("input.jl")
#include("functions.jl")
#include("run_simulation.jl")
#include("output.jl")



"""
    CHARISMA_biomass()

Function to run Charisma without saving output files

Arguments used from settings: none

Returns: Mean summer biomass for all lakes, species, and multiple depths

"""
function CHARISMA_biomass()

        # Get Settings for selection of lakes, species & depth
        cd(dirname(@__DIR__))
        GeneralSettings = parseconfigGeneral("./input/general.config.txt")
        depths = parse.(Float64, GeneralSettings["depths"])
        nyears = parse.(Int64, GeneralSettings["years"])
        nlakes = length(GeneralSettings["lakes"]) #VE
        nspecies = length(GeneralSettings["species"]) #VE

        # Define output structure
        Macroph = zeros(Float64, (nspecies*nlakes), (length(GeneralSettings["depths"])+2))
        j=0 # counter

        # Loop for model run for selected lakes, species and depths
        for l in 1:length(GeneralSettings["lakes"])
            println(GeneralSettings["lakes"][l])

            for s in 1:length(GeneralSettings["species"])
                println(GeneralSettings["species"][s])
                j=j+1 #counter

                #Get settings
                settings = getsettings(GeneralSettings["lakes"][l], GeneralSettings["species"][s])

                #Test if setting are logic; if not break
                if testSettings(settings)!=0
                      break
                end


                push!(settings, "years" => parse.(Int64,GeneralSettings["years"])[1]) #add "years" from GeneralSettings
                push!(settings, "yearsoutput" => parse.(Int64,GeneralSettings["yearsoutput"])[1]) #add "years" from GeneralSettings
                push!(settings, "modelrun" => GeneralSettings["modelrun"][1]) #add "modelrun" from GeneralSettings

                # Simulate environment
                dynamicData = Dict{Int16, DayData}()
                environment = simulateEnvironment(settings, dynamicData)

                # Get macrophytes in multiple depths
                result = simulateMultipleDepth(depths,settings,dynamicData) #Biomass, Number, indWeight, Height,
                # [depths][1=superInd][day*year,parameter]]

                # (select total (seed+tuber) biomass [or all values] of fieldday of last year of simulation for all depths)
                junefirst=nyears[1]*365-(365-152) #N of output day
                augustlast=nyears[1]*365-(365-243) #N of output day

                # Export summer mean of total Biomass for each depth
                for i = 1:4 #depths
                        Macroph[j,i] = mean(result[i][1][junefirst:augustlast,1]) #result[i][1][day,1:4]
                end

                Spec_number_as_string=split(split(GeneralSettings["species"][s],".")[2],"_")[end]
                Spec_number=parse(Int, Spec_number_as_string)
                Lake_number_as_string=split(split(GeneralSettings["lakes"][l],".")[2],"_")[end]
                Lake_number=parse(Int, Lake_number_as_string)
                Macroph[j,5]=Spec_number #species Number GeneralSettings["species"][s]
                Macroph[j,6]=Lake_number #lake Number

            end
        end
    return (Macroph) #Table For all lakes (&species) together
end

#CHARISMA_biomass()




"""
    CHARISMA_biomass_parallel()

Function to run Charisma without saving output files

Arguments used from settings: none

Returns: Mean summer biomass for all lakes, species, and multiple depths

"""
function CHARISMA_biomass_parallel()
        # Get Settings for selection of lakes, species & depth
        cd(dirname(@__DIR__))
        GeneralSettings = parseconfigGeneral("./input/general.config.txt")
        depths = parse.(Float64, GeneralSettings["depths"])
        nyears = parse.(Int64, GeneralSettings["years"])
        nlakes = length(GeneralSettings["lakes"]) #VE
        nspecies = length(GeneralSettings["species"]) #VE

        # Define output structure
        Macroph = zeros(Float64, (nspecies*nlakes), (length(GeneralSettings["depths"])+2))
        lak = zeros(nlakes)

        # Loop for model run for selected lakes, species and depths
        Threads.@threads for l in 1:length(GeneralSettings["lakes"])

            #println(GeneralSettings["lakes"][l])
            lak[l]=l
            for s in 1:length(GeneralSettings["species"])
                #j=j+1 #counter
                j=Int.((lak[l]-1)*length(GeneralSettings["species"]) + s)
                #println(GeneralSettings["species"][s])

                #Get settings
                settings = getsettings(GeneralSettings["lakes"][l], GeneralSettings["species"][s])
                push!(settings, "years" => parse.(Int64,GeneralSettings["years"])[1]) #add "years" from GeneralSettings
                push!(settings, "yearsoutput" => parse.(Int64,GeneralSettings["yearsoutput"])[1]) #add "years" from GeneralSettings
                push!(settings, "modelrun" => GeneralSettings["modelrun"][1]) #add "modelrun" from GeneralSettings

                #Add Spec & Lake number to output
                Spec_number_as_string=split(split(GeneralSettings["species"][s],".")[2],"_")[end]
                Spec_number=parse(Int, Spec_number_as_string)
                Lake_number_as_string=split(split(GeneralSettings["lakes"][l],".")[2],"_")[end]
                Lake_number=parse(Int, Lake_number_as_string)
                Macroph[j,5]=Spec_number #species Number
                Macroph[j,6]=Lake_number #lake Number

                #Test if setting are logic; if not break
                if testSettings(settings)!=0
                      break
                end

                # Simulate environment
                dynamicData = Dict{Int16, DayData}()
                environment = simulateEnvironment(settings, dynamicData)

                # Get macrophytes in multiple depths
                result = simulateMultipleDepth_parallel(depths,settings,dynamicData) #Biomass, Number, indWeight, Height,
                # [depths][1=superInd][day*year,parameter]]

                # Virtual Ecologist
                # (select total (seed+tuber) biomass [or all values] of fieldday of last year of simulation for all depths)
                junefirst=nyears[1]*365-(365-152) #N of output day
                augustlast=nyears[1]*365-(365-243) #N of output day

                # Export summer mean of total Biomass for each depth
                for i = 1:4 #depths
                        Macroph[j,i] = mean(result[i][1][junefirst:augustlast,1]) #result[i][1][day,1:4]
                end

            end
        end
    return (Macroph) #Table For all lakes (&species) together
end

#CHARISMA_biomass_parallel()

#using Profile
#@profile CHARISMA_biomass_parallel()
#@time CHARISMA_biomass_parallel()
#@time CHARISMA_biomass()



"""
    CHARISMA_biomass_parallel_lastNyears()

Function to run Charisma without saving output files

Arguments used from settings: none

Returns: Mean summer biomass for all lakes, species, and multiple depths

"""
function CHARISMA_biomass_parallel_lastNyears()
        # Get Settings for selection of lakes, species & depth
        cd(dirname(@__DIR__))
        GeneralSettings = parseconfigGeneral("./input/general.config.txt")
        depths = parse.(Float64, GeneralSettings["depths"])
        nyears = parse.(Int64, GeneralSettings["years"])
        nlakes = length(GeneralSettings["lakes"]) #VE
        nspecies = length(GeneralSettings["species"]) #VE
        nyearsoutput = parse.(Int64, GeneralSettings["yearsoutput"])[1]


        # Define output structure
        Macroph = zeros(Float64, (nspecies*nlakes*nyearsoutput), (length(GeneralSettings["depths"])+3))
        lak = zeros(nlakes)


        # Loop for model run for selected lakes, species and depths
        Threads.@threads for l in 1:length(GeneralSettings["lakes"])

            #println(GeneralSettings["lakes"][l])
            lak[l]=l
            for s in 1:length(GeneralSettings["species"])

                #println(GeneralSettings["species"][s])

                #Get settings
                settings = getsettings(GeneralSettings["lakes"][l], GeneralSettings["species"][s])
                push!(settings, "years" => parse.(Int64,GeneralSettings["years"])[1]) #add "years" from GeneralSettings
                push!(settings, "yearsoutput" => parse.(Int64,GeneralSettings["yearsoutput"])[1]) #add "years" from GeneralSettings
                push!(settings, "modelrun" => GeneralSettings["modelrun"][1]) #add "modelrun" from GeneralSettings

                #Add Spec & Lake number to output
                Spec_number_as_string=split(split(GeneralSettings["species"][s],".")[2],"_")[end]
                Spec_number=parse(Int, Spec_number_as_string)
                Lake_number_as_string=split(split(GeneralSettings["lakes"][l],".")[2],"_")[end]
                Lake_number=parse(Int, Lake_number_as_string)

                #Test if setting are logic; if not break
                if testSettings(settings)!=0
                      break
                end

                # Simulate environment
                dynamicData = Dict{Int16, DayData}()
                environment = simulateEnvironment(settings, dynamicData)

                # Get macrophytes in multiple depths
                result = simulateMultipleDepth_parallel(depths,settings,dynamicData) #Biomass, Number, indWeight, Height,
                # [depths][1=superInd][day*year,parameter]]

                #Write output
                for y in 0:(nyearsoutput-1)
                    #j=j+1 #geht nicht wegen parallelisierung
                    j=Int.(((lak[l]-1)*length(GeneralSettings["species"])*nyearsoutput) + ((s-1)*nyearsoutput) + (y+1))

                    Macroph[j,5]=Spec_number #species Number TODO for as many lines as years
                    Macroph[j,6]=Lake_number #lake Number TODO for as many lines as years
                    Macroph[j,7]=nyears[1]-y

                    junefirst=(nyears[1]-y)*365-(365-152) #N of output day
                    augustlast=(nyears[1]-y)*365-(365-243) #N of output day

                    # Export summer mean of total Biomass for each depth and year
                    for i = 1:4 #depths
                            Macroph[j,i] = mean(result[i][1][junefirst:augustlast,1]) #result[i][1][day,1:4]
                    end

                end
            end
        end
    return (Macroph) #Table For all lakes (&species) together
end

#CHARISMA_biomass_parallel_lastNyears()

"""
    CHARISMA_biomass_onedepth()

Function to run Charisma without saving output files for one depth

Arguments used from settings: none

Returns: Mean summer biomass for all lakes, species

"""
function CHARISMA_biomass_onedepth()
    #Set dir to home_dir of file

        # Get Settings for selection of lakes, species & depth
        cd(dirname(@__DIR__))
        GeneralSettings = parseconfigGeneral("./input/general.config.txt")
        depths = parse.(Float64, GeneralSettings["depths"])[1] #da nur 1 depth
        nyears = parse.(Int64, GeneralSettings["years"])
        nlakes = length(GeneralSettings["lakes"]) #VE
        nspecies = length(GeneralSettings["species"]) #VE


        # Define output structure
        Macroph = zeros(Float64, (nspecies*nlakes), 3) #1depth,
        j=0 # counter

        # Loop for model run for selected lakes, species and depths
        for l in 1:length(GeneralSettings["lakes"])

            #println(GeneralSettings["lakes"][l])

            for s in 1:length(GeneralSettings["species"])
                j=j+1 #counter
                #println(GeneralSettings["species"][s])

                #Get settings
                settings = getsettings(GeneralSettings["lakes"][l], GeneralSettings["species"][s])
                push!(settings, "years" => parse.(Int64,GeneralSettings["years"])[1]) #add "years" from GeneralSettings
                push!(settings, "yearsoutput" => parse.(Int64,GeneralSettings["yearsoutput"])[1]) #add "years" from GeneralSettings
                push!(settings, "modelrun" => GeneralSettings["modelrun"][1]) #add "modelrun" from GeneralSettings

                # Simulate environment
                dynamicData = Dict{Int16, DayData}()
                environment = simulateEnvironment(settings, dynamicData)

                # Get macrophytes in multiple depths
                result = simulate1Depth(depths,settings, dynamicData) #Biomass, Number, indWeight, Height,
                # [depths][1=superInd][day*year,parameter]]

                # Virtual Ecologist
                # (select total (seed+tuber) biomass [or all values] of fieldday of last year of simulation for all depths)
                junefirst=nyears[1]*365-(365-152) #N of output day
                augustlast=nyears[1]*365-(365-243) #N of output day

                # Export summer mean of total Biomass for each depth
                #for i = 1:4 #depths
                Macroph[j,1] = mean(result[1][junefirst:augustlast,1]) #result[i][1][day,1:4]
                #end

                Macroph[j,2]=s #species Number
                Macroph[j,3]=l #lake Number

            end
        end
    return (Macroph) #Table For all lakes (&species) together
end

#using Profile
#@time CHARISMA_biomass_onedepth()




"""
    CHARISMA_biomass_N_weight_hight_env()

Function to run Charisma without saving output files,

Arguments used from settings: none

Returns: Daily Biomass, Number of Individuals, indWeight, Height, for all lakes, species, and multiple depths in the last year of started simulation

"""



function CHARISMA_biomass_N_weight_hight_env()

        # Get Settings for selection of lakes, species & depth
        cd(dirname(@__DIR__))
        GeneralSettings = parseconfigGeneral("./input/general.config.txt")
        depths = parse.(Float64, GeneralSettings["depths"])
        nyears = parse.(Int64, GeneralSettings["years"])
        nlakes = length(GeneralSettings["lakes"]) #VE
        nspecies = length(GeneralSettings["species"]) #VE
        ndepths = length(depths)

        # Define output structure
        MacrophAll = []

        j=0 # counter

        # Loop for model run for selected lakes, species and depths
        for l in 1:length(GeneralSettings["lakes"])
            println(GeneralSettings["lakes"][l])

            for s in 1:length(GeneralSettings["species"])
                println(GeneralSettings["species"][s])
                j=j+1 #counter

                #Get settings
                settings = getsettings(GeneralSettings["lakes"][l], GeneralSettings["species"][s])

                #Test if setting are logic; if not break
                if testSettings(settings)!=0
                      break
                end

                push!(settings, "years" => parse.(Int64,GeneralSettings["years"])[1]) #add "years" from GeneralSettings
                push!(settings, "yearsoutput" => parse.(Int64,GeneralSettings["yearsoutput"])[1]) #add "years" from GeneralSettings
                push!(settings, "modelrun" => GeneralSettings["modelrun"][1]) #add "modelrun" from GeneralSettings

                # Simulate environment
                dynamicData = Dict{Int16, DayData}()
                environment = simulateEnvironment(settings, dynamicData)

                # Get macrophytes in multiple depths
                result = simulateMultipleDepth(depths,settings,dynamicData) #Biomass, Number, indWeight, Height,
                # [depths][1=superInd][day*year,parameter]]

                # Extract result of last year
                Res = []
                for d in 1:ndepths
                    push!(Res, result[d][1][(((nyears[1]-1)*365)+1):(nyears[1]*365),1:4])
                end

                # Combine outputs
                push!(MacrophAll,Res)
                push!(MacrophAll, environment) # add environment after loop for all species
            end

        end
    return (MacrophAll) #Table For all lakes (&species) together
end

# CHARISMA_biomass_N_weight_hight_env()

# Struktur (Output Plant1 in Lake1)
#          (Environment of lake1) [temp] [irradiance] [waterlevel] [lightAttenuation]
# ...
#          (Output PlantN in Lake1) [depths] [days, Biomass, Number of Individuals, indWeight, Height]
#          (Environment of lake1) [temp] [irradiance] [waterlevel] [lightAttenuation]

#          (Output Plant1in Lake2)
#          (Environment of lake2)
# ...
#          (Output PlantN in Lake2)
#          (Environment of lake2)
# ....


# Information inputs of vactors do not yet work
function CHARISMA_biomass_N_weight_hight_env2(lake,species,depth,years)

        # Get Settings for selection of lakes, species & depth
        cd(dirname(@__DIR__))
        GeneralSettings = parseconfigGeneral("./input/general.config.txt")
        GeneralSettings["depths"]= [depth]
        GeneralSettings["lakes"]= [lake]
        GeneralSettings["years"]= [years]
        GeneralSettings["species"]= [species]

        depths = parse.(Float64, GeneralSettings["depths"])
        nyears = parse.(Int64, GeneralSettings["years"])
        nlakes = length(GeneralSettings["lakes"]) #VE
        nspecies = length(GeneralSettings["species"]) #VE
        ndepths = length(depths)

        # Define output structure
        MacrophAll = []

        j=0 # counter

        # Loop for model run for selected lakes, species and depths
        for l in 1:length(GeneralSettings["lakes"])
            println(GeneralSettings["lakes"][l])

            for s in 1:length(GeneralSettings["species"])
                println(GeneralSettings["species"][s])
                j=j+1 #counter

                #Get settings
                settings = getsettings(GeneralSettings["lakes"][l], GeneralSettings["species"][s])

                #Test if setting are logic; if not break
                if testSettings(settings)!=0
                      break
                end

                push!(settings, "years" => parse.(Int64,GeneralSettings["years"])[1]) #add "years" from GeneralSettings
                push!(settings, "yearsoutput" => parse.(Int64,GeneralSettings["yearsoutput"])[1]) #add "years" from GeneralSettings
                push!(settings, "modelrun" => GeneralSettings["modelrun"][1]) #add "modelrun" from GeneralSettings

                # Simulate environment
                dynamicData = Dict{Int16, DayData}()
                environment = simulateEnvironment(settings, dynamicData)

                # Get macrophytes in multiple depths
                result = simulateMultipleDepth(depths,settings,dynamicData) #Biomass, Number, indWeight, Height,
                # [depths][1=superInd][day*year,parameter]]

                # Extract result of last year
                Res = []
                for d in 1:ndepths
                    push!(Res, result[d][1][(((nyears[1]-1)*365)+1):(nyears[1]*365),1:4])
                end

                # Combine outputs
                push!(MacrophAll,Res)
                push!(MacrophAll, environment)
            end
             # add environment after loop for all species
        end
    return (MacrophAll) #Table For all lakes (&species) together
end


#CHARISMA_biomass_N_weight_hight_env2("./input/lakes/lake_1.config.txt","./input/species/species_3.config.txt","-1.0","5")









"""
    CHARISMA_test_15lakes_4depths()

Function to testrun Charisma without saving output files for one depth

Arguments used from settings: none

Returns: Mean summer biomass for all lakes, species

"""
function CHARISMA_test_15lakes_4depths()

    Macroph = zeros(Float64, 15, 6)

    Macroph[:,1:4].= rand(0:500,15,4)
    Macroph[:,5].=1 #species Number
    Macroph[:,6].=collect(1:15) #lake Number

    return(Macroph)
end


function CHARISMA_parallel_test_15lakes_4depths()

    #cd(dirname(@__DIR__))
    #cd("model")

    # Include functions


    # Get Settings for selection of lakes, species & depth
    cd(dirname(@__DIR__))
    GeneralSettings = parseconfigGeneral("./input/general.config.txt")
    depths = parse.(Float64, GeneralSettings["depths"])
    nyears = parse.(Int64, GeneralSettings["years"])
    nlakes = length(GeneralSettings["lakes"]) #VE
    nspecies = length(GeneralSettings["species"]) #VE

    Macroph = zeros(Float64, (nspecies*nlakes), (length(GeneralSettings["depths"])+2))

    lak = zeros(nlakes)

    Threads.@threads for l in 1:length(GeneralSettings["depths"])
        Macroph[:,l].= rand(0:500,15)
    end
    Macroph[:,5].=1 #species Number
    Macroph[:,6].=collect(1:15) #lake Number

    return(Macroph)
end


"""
    CHARISMA_VE()

Function to run Charisma and the Virtual Ecologist Approach without saving output files

!OUT OF USE!

Arguments used from settings: none

Returns: Table with Kohler numbers for x depths (first x rows) for all lakes and species (each combination is a row)


function CHARISMA_VE()
    #Set dir to home_dir of file
        #cd(dirname(@__DIR__))
        #cd("model")

        # Include functions
        #include("defaults.jl")
        #include("input.jl")
        #include("functions.jl")
        #include("run_simulation.jl")
        #include("output.jl")

        # Get Settings for selection of lakes, species & depth
        cd(dirname(@__DIR__))
        GeneralSettings = parseconfigGeneral("./input/general.config.txt")
        depths = parse.(Float64, GeneralSettings["depths"])
        nyears = parse.(Int64, GeneralSettings["years"])
        nlakes = length(GeneralSettings["lakes"]) #VE
        nspecies = length(GeneralSettings["species"]) #VE

        # Set VE parameters - have to be output in the end
        pFindSpecies = rand(Uniform(0.5,1.0),nlakes)
        fieldday = 226 #

        # Define output structure
        mappedMacroph = zeros(Float64, (nspecies*nlakes), (length(GeneralSettings["depths"])+2))
        j=0 # counter

        # Loop for model run for selected lakes, species and depths
        for l in 1:length(GeneralSettings["lakes"])

            println(GeneralSettings["lakes"][l])

            for s in 1:length(GeneralSettings["species"])
                j=j+1 #counter
                println(GeneralSettings["species"][s])

                #Get settings
                settings = getsettings(GeneralSettings["lakes"][l], GeneralSettings["species"][s])
                push!(settings, "years" => parse.(Int64,GeneralSettings["years"])[1]) #add "years" from GeneralSettings
                push!(settings, "yearsoutput" => parse.(Int64,GeneralSettings["yearsoutput"])[1]) #add "years" from GeneralSettings
                push!(settings, "modelrun" => GeneralSettings["modelrun"][1]) #add "modelrun" from GeneralSettings

                # Get macrophytes in multiple depths
                result = simulateMultipleDepth(depths,settings) #Biomass, Number, indWeight, Height,
                # [depths][1=superInd][day*year,parameter]]

                # Virtual Ecologist
                # (select total (seed+tuber) biomass [or all values] of fieldday of last year of simulation for all depths)
                day=nyears[1]*365-(365-fieldday) #N of output day

                #Create Data structure & import data
                BiomNHeight = zeros(4,4) #biomass, N, indWeight, Height =rows; 4depths = columns

                for i = 1:4
                        BiomNHeight[i,:] = result[i][1][day,1:4]
                end

                # Check if Biomass < 0.01 and N<1 -> cannot be found
                for i = 1:4
                            if BiomNHeight[1,i]<0.01 && BiomNHeight[2,i]<1
                                    for j = 1:4
                                            BiomNHeight[j,i]=0
                                    end
                            end
                end

                # Find with probability
                #Create Data structure
                BiomNHeightMapped = zeros(4,4) #biomass, N, indWeight, Height =rows; 4depths = columns

                # Find with probability
                found =sample([1, 0], Weights([pFindSpecies[nlakes], 1-pFindSpecies[nlakes]]), 1)[1]
                for i=1:4 #TODO ndepths
                        if found == 1
                                BiomNHeightMapped[1,i]=BiomNHeight[1,i]
                        end
                        for j = 2:4
                                if BiomNHeightMapped[1,i] == 0
                                        BiomNHeightMapped[j,i] == 0
                                else BiomNHeightMapped[j,i] =BiomNHeight[j,i]
                                end
                        end
                end

                mappedMacroph[j,5]=s
                mappedMacroph[j,6]=l

                # Transform into Kohler value
                Kohler5_s = settings["Kohler5"]
                Kohler4 = Kohler5_s*64/125
                Kohler3 = Kohler5_s*27/125
                Kohler2 = Kohler5_s*8/125

                for d = 1:4
                        if BiomNHeightMapped[d,1]==0
                                mappedMacroph[j,d] = 0
                                #
                        elseif BiomNHeightMapped[d,1]>Kohler5_s
                                mappedMacroph[j,d] = 5
                        elseif BiomNHeightMapped[d,1]>Kohler4
                                mappedMacroph[j,d] = 4
                        elseif BiomNHeightMapped[d,1]>Kohler3
                                mappedMacroph[j,d] = 3
                        elseif BiomNHeightMapped[d,1]>Kohler2
                                mappedMacroph[j,d] = 2
                        elseif BiomNHeightMapped[d,1]>0
                                mappedMacroph[j,d] = 1
                        end
                end

            end
            # Combine result for all lakes
        end
    return (mappedMacroph) #For all lakes (&species) together
end

#CHARISMA_VE()
"""
