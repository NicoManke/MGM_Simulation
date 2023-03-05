# returns an array containing the individual amount of lakes and species
function count_lakes_and_species()

    #cd(dirname(@__DIR__)) # setting the working directory to the current directory 
	Lengths = Int64[]
	Config = parseconfigGeneral("./input/general.config.txt")
	#mein absoluter Pfad. Nur für Tests ;)
	#("D:WS21/GameLabII-WS21/Plugin/code/UnityCode/JuliaPlugin/input/general.config.txt")
	
	push!(Lengths, length(Config["lakes"]))
	push!(Lengths, length(Config["species"]))
	
	return Lengths
end
	