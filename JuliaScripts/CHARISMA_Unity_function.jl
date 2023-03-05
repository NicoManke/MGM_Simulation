# altered version of the CHARISMA_biomass_N_weight_hight_env() function which allows using two indices to determine upfront the lake and species that should be calculated
function CHARISMA_Unity(lakeNumber, speciesNumber)

    # Get Settings for selection of lakes, species & depth
    #cd(dirname(@__DIR__))
    GeneralSettings = parseconfigGeneral("./input/general.config.txt")
    depths = parse.(Float64, GeneralSettings["depths"])
    nyears = parse.(Int64, GeneralSettings["years"])
    nlakes = length(GeneralSettings["lakes"]) #VE
    nspecies = length(GeneralSettings["species"]) #VE
    ndepths = length(depths)

	indexLakes = lakeNumber
	indexSpecies = speciesNumber
	
	# Test if passed index is valid
	if (indexLakes > nlakes)
		println("Chosen lake index was out of bounds! Instead the index of the last available lake was used: ", nlakes)
		indexLakes = nlakes
	end
	
	# Test if passed index is valid
	if (indexSpecies > nspecies)
		println("Chosen species index was out of bounds! Instead the index of the last available species was used: ", nspecies)
		indexSpecies = nspecies
	end

    # Define output structure
    MacrophAll = []

    # Debug selected lake and species
    println(GeneralSettings["lakes"][indexLakes])
    println(GeneralSettings["species"][indexSpecies])

    # Get settings for selected lake and species 
    settings = getsettings(GeneralSettings["lakes"][indexLakes], GeneralSettings["species"][indexSpecies])

    # Test if setting are logic; if not break
    if !(testSettings(settings) != 0)

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
		push!(MacrophAll, Res)
		push!(MacrophAll, environment) # add environment after loop for all species
	end
	
    return (MacrophAll) #Table For all lakes (&species) together
end
