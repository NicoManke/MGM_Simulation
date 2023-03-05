#loads the dependencies required by the model
try 
	using HCubature
catch e
    import Pkg; Pkg.add("HCubature")
	using HCubature
end
	
try 
	using DelimitedFiles
catch e
    import Pkg; Pkg.add("DelimitedFiles")
	using DelimitedFiles
end

try 
	using Dates
catch e
    import Pkg; Pkg.add("Dates")
	using Dates
end

try 
	using Distributions
catch e
    import Pkg; Pkg.add("Distributions")
	using Distributions
end

try 
	using CSV
catch e
    import Pkg; Pkg.add("CSV")
	using CSV
end

#try 
#	using Random
#catch e
#	try
#		import Pkg; Pkg.add("Random")
#		using Random	
#	catch e
#		# it shouldn't be a problem, if this package is missing
#	end
#end

	#HCubature, #for Integration
	#DelimitedFiles, # for function writedlm, used to write output files
	#Dates, #to create output folder
	#Distributions, Random, #for killWithProbability
	#CSV#, #For virtual Ecologist
	#DataFrames, #For virtual Ecologist
	#StatsBase #For virtual Ecologist
