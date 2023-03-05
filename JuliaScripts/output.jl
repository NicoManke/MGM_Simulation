#Macrophyte Growth Model (MGM)
#
#Anne Lewerentz <anne.lewerentz@uni-wuerzburg.de>
#(c) 2021-2022, licensed under the terms of the MIT license
#
#Contains functions that are necessary to write output files

"""
    writeOutputMacrophytes(PlantResults)

Output functions for normal modelrun

Code source: Leidinger et al. (2021)

"""

function writeOutputMacrophytes(PlantResults,depths)
    #homdir = pwd()
    #dirname = settings["Lake"] * "_" *settings["Species"] * "_" *Dates.format(now(), "yyyy_m_d_HH_MM") # "folder" * * string(settings["latitude"])
    #cd(".\\output")
    #mkdir(dirname)
    #cd(dirname)

    for i in 1:length(depths)
        j=depths[i]
        writedlm("superInd$j.csv", PlantResults[i][1][:, :, 1], ',') #Biomass, Number, indWeight, Height, allocatedBiomass, SpreadBiomass
        writedlm("superIndSeed$j.csv", PlantResults[i][2][:, :, 1], ',') #Biomass, Number, indWeight, Height, allocatedBiomass, SpreadBiomass
        writedlm("superIndTuber$j.csv", PlantResults[i][3][:, :, 1], ',') #Biomass, Number, indWeight, Height, allocatedBiomass, SpreadBiomass
        writedlm("seeds$j.csv", PlantResults[i][4][:, :, 1], ',') #Biomass, Number, indWeight, Height, allocatedBiomass, SpreadBiomass
        writedlm("tubers$j.csv", PlantResults[i][5][:, :, 1], ',') #Biomass, Number, indWeight, Height, allocatedBiomass, SpreadBiomass
        writedlm("growthSeeds$j.csv", PlantResults[i][6][:, :, 1], ',') #Biomass, Number, indWeight, Height, allocatedBiomass, SpreadBiomass
        writedlm("growthTubers$j.csv", PlantResults[i][7][:, :, 1], ',') #Biomass, Number, indWeight, Height, allocatedBiomass, SpreadBiomass
    end
    #cd(homdir)
end



"""
    writeOutputMacrophytesLastXYears(PlantResults)

Output functions for normal modelrun - gives out last 5 years of simulation

Code source: Leidinger et al. (2021)

"""

function writeOutputMacrophytesLastXYears(settings, PlantResults, depths, Nyears)
    #homdir = pwd()
    #dirname = settings["Lake"] * "_" *settings["Species"] * "_" *Dates.format(now(), "yyyy_m_d_HH_MM") # "folder" * * string(settings["latitude"])
    #cd(".\\output")
    #mkdir(dirname)
    #cd(dirname)

    for i in 1:length(depths)
        j=depths[i]
        writedlm("superInd$j.csv", PlantResults[i][1][(Nyears*settings["yearlength"]-settings["yearsoutput"]*settings["yearlength"]):(Nyears*settings["yearlength"]), :, 1], ',') #Biomass, Number, indWeight, Height, allocatedBiomass, SpreadBiomass
        writedlm("superIndSeed$j.csv", PlantResults[i][2][(Nyears*settings["yearlength"]-settings["yearsoutput"]*settings["yearlength"]):(Nyears*settings["yearlength"]), :, 1], ',') #Biomass, Number, indWeight, Height, allocatedBiomass, SpreadBiomass
        writedlm("superIndTuber$j.csv", PlantResults[i][3][(Nyears*settings["yearlength"]-settings["yearsoutput"]*settings["yearlength"]):(Nyears*settings["yearlength"]), :, 1], ',') #Biomass, Number, indWeight, Height, allocatedBiomass, SpreadBiomass
        writedlm("seeds$j.csv", PlantResults[i][4][(Nyears*settings["yearlength"]-settings["yearsoutput"]*settings["yearlength"]):(Nyears*settings["yearlength"]), :, 1], ',') #Biomass, Number, indWeight, Height, allocatedBiomass, SpreadBiomass
        writedlm("tubers$j.csv", PlantResults[i][5][(Nyears*settings["yearlength"]-settings["yearsoutput"]*settings["yearlength"]):(Nyears*settings["yearlength"]), :, 1], ',') #Biomass, Number, indWeight, Height, allocatedBiomass, SpreadBiomass
        writedlm("growthSeeds$j.csv", PlantResults[i][6][(Nyears*settings["yearlength"]-settings["yearsoutput"]*settings["yearlength"]):(Nyears*settings["yearlength"]), :, 1], ',') #Biomass, Number, indWeight, Height, allocatedBiomass, SpreadBiomass
        writedlm("growthTubers$j.csv", PlantResults[i][7][(Nyears*settings["yearlength"]-settings["yearsoutput"]*settings["yearlength"]):(Nyears*settings["yearlength"]), :, 1], ',') #Biomass, Number, indWeight, Height, allocatedBiomass, SpreadBiomass
    end
    #cd(homdir)
end



"""
    writeOutputEnvironmentSettings(Env,Settings,)

Code source: Leidinger et al. (2021)

"""

function writeOutputEnvironmentSettings(
    Env,
    Settings,
)
    #homdir = pwd()
    #dirname = settings["Lake"] * "_" *settings["Species"] * "_" *Dates.format(now(), "yyyy_m_d_HH_MM") # "folder" * * string(settings["latitude"])
    #cd(".\\output")
    #mkdir(dirname)
    #cd(dirname)
    writedlm("Temp.csv", Env[1], ',')
    writedlm("Irradiance.csv", Env[2], ',')
    writedlm("Waterlevel.csv", Env[3], ',')
    writedlm("lightAttenuation.csv", Env[4], ',')

    writedlm("Settings.csv", Settings, ',')

    #cd(homdir)
end



"""
    writeGeneralSettings(Env,Settings,)

Code source: Leidinger et al. (2021)

"""

function writeGeneralSettings(
    GeneralSettings,
)
    #homdir = pwd()
    #dirname = settings["Lake"] * "_" *settings["Species"] * "_" *Dates.format(now(), "yyyy_m_d_HH_MM") # "folder" * * string(settings["latitude"])
    #cd(".\\output")
    #mkdir(dirname)
    #cd(dirname)
    writedlm("GeneralSettings.csv", GeneralSettings, ',')
    #cd(homdir)
end


"""
    writeOutput(settings::Dict{String, Any}, Env,PlantResults)

Creates the output directory and copies relevant files into it.

Code source: Leidinger et al. (2021)

"""

function writeOutput(settings::Dict{String, Any}, depth, Env, PlantResults, GeneralSettings,dest)
    homdir = pwd()
    #isdir(output) || mkdir(output)
    cd("./output")

    dirname = dest
    if isdir(dirname)
        #@warn "$(settings["modelrun"]) exists. Continuing anyway. Overwriting of files possible."
    else
        mkpath(dirname)
    end

    cd(dirname)
    writeGeneralSettings(GeneralSettings)
    dirname2 = settings["Lake"] * "_" *settings["Species"]
    if isdir(dirname2)
        #@warn "$(settings["modelrun"]) exists. Continuing anyway. Overwriting of files possible."
    else
        mkpath(dirname2)
    end
    cd(dirname2)
    #simlog("Setting up output directory $(settings["dest"])", settings)
    writeOutputEnvironmentSettings(Env, settings)
    #writeOutputMacrophytes(PlantResults, depths)
    writeOutputMacrophytesLastXYears(settings, PlantResults, depth, settings["years"])
    cd(homdir)

end
