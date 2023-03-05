#Macrophyte Growth Model (MGM)
#
#Anne Lewerentz <anne.lewerentz@uni-wuerzburg.de>
#(c) 2021-2022, licensed under the terms of the MIT license
#
#Contains all default settings


"""
Default global settings

Defines the list of configuration variables and returns their default values
in a Dict.
Code Source: Leidinger et al. (2021)
"""
function defaultSettingsGlobal()
    Dict(
        "yearlength" => 365, #Number of days per year [n]
        #"dest"  => string(Dates.format(now(), "yyyy_m_d_HH_MM")), #actual date
        )
end

"""
Default general settings

Defines the list of configuration variables and returns their default values
in a Dict.
Code Source: Leidinger et al. (2021)
"""
function defaultSettingsGeneral()
    Dict(
        "years" => 1, #Number of years to get simulated [n]
        "depths" => [-0.5,-1.0],
        "yearsoutput" => 2, #Number of years to get output for
        "species" => (".\\input\\species\\CharaAspera_1.config.txt"),
        "lakes" => (".\\input\\lakes\\TurbidWarmLakeNutrientrich.config.txt"),
        "modelrun" => "test", #name of outputfolder
        )
end

"""
Default settings for the environment

Defines the list of configuration variables and returns their default values
in a Dict.
Code Source: Leidinger et al. (2021)
"""
function defaultSettingsLake()
    # Return the default settings. All parameters must be registered here.
    Dict(
        #GENERAL
        "Lake" => "default", #Name of the lake

        #CARBONATE
        #"maxCarbonate" => 200.0, #Carbonate in water without plants [mg/l]

        #LIGHT
        "fracReflected" => 0.1, #Fraction of light reflectad at the water surface [fraction]
        "iDelay" => -10, #Days after 1st of January where I is minimal [d]
        "iDev" => 1.0, #Deviation factor to change total irradiation [fraction]
        "latitude" => 47.5, #Latitude of corresponding lake [°]
        "maxI" => 1500.0, #Maximal Irradiance [µE m^-2 s^-1]; 868 in CHARISMA
        "minI" => 150.0, #Minimal Irradiance [µE m^-2 s^-1]; 96 in CHARISMA
        "parFactor" => 0.5, #Fraction of total irradiation that is PAR [fraction]

        #NUTRIENT
        "maxNutrient" => 0.005, #Concentration of limiting nutrient in water without plants [mg l^-1]

        #TEMPERATURE
        "maxTemp" => 18.8, #Maximal mean daily temperature of a year in [°C]
        "minTemp" => 0.0, #Minimal mean daily temperature of a year in [°C]
        "tempDelay" => 23, #Days after 1st of January where Temp is minimal [d]
        "tempDev" => 1.0, #Factor for deviation of temperatur [fraction]

        #VERTUCAL LIGHT ATTENUATION / TURBIDITY
        "backgrKd" => 1.0, #Background light attenuation of water (Vertical light attenuation, turbidity) [m^-1]
        #"clearWaterFraction"
        #"clearWaterPeriod"
        #"clearWaterTiming"
        #"kd" => 2.0, #For alternative way to define cosine: Mean light attenuation coefficient (Kd) (cosine) [m^-1]
        "kdDelay" => -10.0, #Delay in cosine, the day number with the minimal light attenuation coefficient [d]
        "kdDev" => 1.0, #Deviation factor, a factor between 0 and 1 to change the whole light attenuation range [-]
        #"kdDiffusion" #for spatial process
        #"kdRange" #For alternative way to define cosine
        #"KdStochastic"
        "maxKd" => 8.0, #Maximum light attenuation coefficient [m^-1]
        "minKd" => 2.0, #Minimum light attenuation coefficient [m^-1]

        # WATER LEVEL
        "levelCorrection" => 0.0, #Correction for reference level [m]
        "maxW" => 0.0, #Maximal water level [m]
        "minW" => -0.0, #Minimal water level [m]
        #"WaterChange"
        #"WaterChangePeriod"
        #...
        "wDelay" => 280, #Delay of cosine of water level [m]
        #"wDev"
    )
end

"""
Default settings for the species

Defines the list of configuration variables and returns their default values
in a Dict.
Code Source: Leidinger et al. (2021)
"""
function defaultSettingsSpecies()
    # Return the default settings. All parameters must be registered here.
    Dict(
        #GENERAL
        "Species" => "default", #Name of the species

        #BIOMASS PARTIONING
        "seedsEndAge" => 60, #Age of the plants where seedFraction is reached [d]
        "seedsStartAge" => 30, #Age of the plants where seed formation starts [d]
        "tuberEndAge" => 60, #Age of the plants where tuberFraction is reached [d]
        "tuberStartAge" => 30, #Age of the plants where tuber formation starts [d]

        #CARBONATE
        #"hCarbonate" => 30, #Halfsaturation carbonate concentration of growth [mg/l]
        #"hCarboReduction" => 30, #Halfsaturation biomass of carbonate reduction by plants [g m^-1]
        #"pCarbonate" => 1, #Power of Hill function of carbonate dependent growth [-]

        #GROWTH
        "cTuber" => 0.1, #Fraction of tuber weight lost daily when sprout starts growing [fraction]
        "pMax" => 1.7, #Maximal gross photosynthesis=specific daily production of the plant top at 20°C in the absence of light limitation [g g^-1 h^-1] | 0.006 in CHARISMA
        "q10" => 2.0, #Q10 for maintenance respiration [-]
        "resp20" => 0.00193, #Respiration at 20°C [d^-1]

        #GROWTH FORM
        "heightMax" => 2.35,  #maximum length of macrophyte [m]
        "maxWeightLenRatio" => 0.001, #Weight of 1m young sprout [g m^-1]
        "rootShootRatio" => 0.1, #Proportion of plant allocated to the roots [fraction]
        #"spreadFrac" => 0.7, #Fraction of sprout weight spreaded under the water surface [fraction]

        #LIGHT
        "fracPeriphyton" => 0.2, #Fraction of light reduced by periphyton [fraction]
        "hPhotoDist" => 1.0, #Distance from plant top at which the photosy. is reduced factor 2 [m]
        "hPhotoLight" => 14.0, #Half saturation light intensity (PAR) for photosynthesis [µE m^-2 s^-1]
        "hPhotoTemp" => 14.0, #Half saturation temperature for photosynthesis [°C]
        "hTurbReduction" => 40.0, #Half saturation biomass of light attenuation reduction [g m^-1]
        "plantK" => 0.02, #Light attenuation of plant tissue [m^2 g^-1]
        "pPhotoTemp" => 3.0, #Exponent in temp effect (Hill function) for photosynthesis [-]
        "pTurbReduction" => 1.0, #Exponent in Hill function of light attenuation reduction [-]
        "sPhotoTemp" => 1.35,  #Scaling of temperature effect for photosynthesis [-]

        #MORTALITY
        "BackgroundMort" => 0.05, #Background mortality [mort d^-1]
        "cThinning" => 5950, #Constant of the thinning rule [-]
        "hWaveMort" => 0.1, #Half saturation depth for mortality [m]
        "maxAge" => 175, #Maximum age [d]
        #"maxDryDays" #Maximal dry days [d]
        "maxWaveMort" => 0.2, #Maximum loss of weight in shallow areas [fraction]
        "pWaveMort" => 4, #Power of Hill function for wave mortality [-]
        #"ThinAdjWeight" #Adjust the weight in the thinning rule [-]
        "thinning" => true, #Apply the thinning rule? [-]

        #NUTRIENT
        "hNutrient" => 0.0, #Halfsaturation nutrient concentration of growth [mg l^-1]
        #"hNutrReduction" => 200.0, #Half saturation biomass of nurtient reduction by plants [g m^-2]
        "pNutrient" => 1, #Power of Hill function of nutrient dependent growth [-]

        #REPRODUCTION
        "germinationDay" => 114,  #Day of germination of seeds [d]
        "reproDay" => 250, #Day of dispersal of seeds [d]
        "seedBiomass" => 0.00002, #Individual weight of seeds [g]
        "seedFraction" => 0.13, #Fraction of plant weight allocated to seeds [y^-1]
        "seedGermination" => 0.2, #Fraction of seeds that germinate [y^-1]
        #"SeedGrazingThres"
        #"SeedImport" => 0.1, #Import of seeds [N m^-2 y^-1]
        "seedInitialBiomass" => 2.0, #Initial biomass of seeds [g m^-2]
        "seedMortality" => 0.0, #Mortality of seeds [mort d^-1]
        #"SeedRadius"

        "tuberBiomass" => 0.00002,  #Individual weight of tubers [g]
        "tuberFraction" => 0.22, #Fraction of plant weight allocated to tubers [y^-1]
        "tuberGermination" => 1.0, #Fraction of tubers that germinate [y^-1]
        "tuberGerminationDay" => 114, #The day that tubers germinate [d]
        #"TuberGrazingThres"
        #"TuberImport"
        "tuberInitialBiomass" => 1.0, #Initial biomass of tubers [g m^-2]
        "tuberMortality" => 0.0, #Mortality of tubers [mort d^-1]

        #virtualEcologist parameter
        "Kohler5" => 100.0,
    )
end
