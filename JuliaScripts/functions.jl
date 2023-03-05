#Macrophyte Growth Model (MGM)
#
#Anne Lewerentz <anne.lewerentz@uni-wuerzburg.de>
#(c) 2021-2022, licensed under the terms of the MIT license
#
#Contains cifferent functions that are needed to run MGM

"""
    getDaylength(day; settings; dynamicData)

Takes a day and calculates for a distinct latitude the daylength.

Source:  R III Author: Robert J. Hijmans, r.hijmans@gmail.com # License GPL3 # Version 0.1  January 2009
Forsythe, William C., Edward J. Rykiel Jr., Randal S. Stahl, Hsin-i Wu and Robert M. Schoolfield, 1995.
A model comparison for daylength as a function of latitude and day of the year. Ecological Modeling 80:87-95.

Arguments used from settings: latitude

Return: daylength [h]
"""
function getDaylength(day, settings::Dict{String, Any}, dynamicData::Dict{Int16, DayData})
    if ismissing(dynamicData[day].daylength)
        settings["latitude"] > 90.0 ||
            settings["latitude"] < -90.0 &&
                return error("latitude must be between 90.0 & -90.0 Degree")
        p = asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (day - 186)))))
        a =
            (sin(0.8333 * pi / 180) + sin(settings["latitude"] * pi / 180) * sin(p)) /
            (cos(settings["latitude"] * pi / 180) * cos(p))
        if a < -1
            a = -1
        elseif a > 1
            a = 1
        end
        dynamicData[day].daylength = 24 - (24 / pi) * acos(a)
    end
    return (dynamicData[day].daylength) #[h]
end


"""
    getTemperature(day; settings; dynamicData)

Temperature gets modeled with a cosine function

Source: van Nes et al. (2003)

Arguments used from settings: yearlength,tempDev,maxTemp,minTemp,tempDelay

Result: Daily water temperature [°C]
"""
function getTemperature(day, settings::Dict{String, Any}, dynamicData::Dict{Int16, DayData})
    if ismissing(dynamicData[day].temperature)
        dynamicData[day].temperature =
            settings["tempDev"] * (
                settings["maxTemp"] -
                ((settings["maxTemp"] - settings["minTemp"]) / 2) *
                (1 + cos((2 * pi / settings["yearlength"]) * (day - settings["tempDelay"])))
            )
    end
    return (dynamicData[day].temperature)
end



"""
    getSurfaceIrradianceDay(day; settings; dynamicData)

Modeled with cosine function

Source: van Nes et al. (2003)

Arguments used from settings:  yearlength, maxI, minI, iDelay

Result: daily SurfaceIrradiance [μE m^-2 s^-1]
"""
function getSurfaceIrradianceDay(day, settings::Dict{String, Any}, dynamicData::Dict{Int16, DayData})
    if ismissing(dynamicData[day].surfaceIrradiance)
        dynamicData[day].surfaceIrradiance =
            settings["maxI"] - (
                ((settings["maxI"] - settings["minI"]) / 2) *
                (1 + cos((2 * pi / settings["yearlength"]) * (day - settings["iDelay"])))
            )
    end
    return (dynamicData[day].surfaceIrradiance)
end



"""
    getWaterlevel(day; settings; dynamicData)

Modeled with cosine function

Arguments used from settings: yearlength, maxW, minW, wDelay, levelCorrection

Result: Waterlevel above / below mean water level [m]
"""
function getWaterlevel(day, settings::Dict{String, Any}, dynamicData::Dict{Int16, DayData})
    if ismissing(dynamicData[day].waterlevel)
        dynamicData[day].waterlevel =
            - settings["levelCorrection"] + (
                settings["maxW"] -
                (settings["maxW"] - settings["minW"]) / 2 *
                (1 + cos((2 * pi / settings["yearlength"]) * (day - settings["wDelay"])))
            )
    end
    return (dynamicData[day].waterlevel) #[m]
end



"""
    reduceNutrientConcentration(Biomass; settings)

Reduction of Nutrient content if there is vegetation
## NOT USED

Source: van Nes et al. (2003)

Arguments used from settings: maxNutrient, hNutrReduction

Result: NutrientConcAdj [mg / l]
"""
function reduceNutrientConcentration(Biomass, settings::Dict{String, Any})
    NutrientConcAdj =
        settings["maxNutrient"] * settings["hNutrReduction"] /
        (settings["hNutrReduction"] + Biomass)
    return (NutrientConcAdj) #[mg / l]
end


"""
    getSurfaceIrradianceHour(day, hour; settings; dynamicData)

Daily total irradiation modeled as sine wave over the year

Source: van Nes et al. (2003)

Arguments used from settings: yearlength,latitude,maxI,minI,iDelay

Result: SurfaceIrradianceHour [μE m^-2 s^-1]
"""
function getSurfaceIrradianceHour(day, hour::Int64, settings::Dict{String, Any}, dynamicData::Dict{Int16, DayData}) #times in hour after sunset
    if ! (hour in keys(dynamicData[day].surfaceIrradianceHour))
        irradianceD = getSurfaceIrradianceDay(day, settings, dynamicData)
        daylength = getDaylength(day, settings, dynamicData)
        dynamicData[day].surfaceIrradianceHour[hour] =
            ((pi * irradianceD) / (2 * daylength)) * sin((pi * hour) / daylength)
    end
    return (dynamicData[day].surfaceIrradianceHour[hour]) #[μE m^-2 s^-1]
end

#getSurfaceIrradianceHour(120,5,settings,dynamicData)

"""
    getLightAttenuation(day; settings; dynamicData)

Modeled with a cosine function

Source: van Nes et al. (2003)

Arguments used from settings:kdDev, maxKd, minKd, yearlength, kdDelay

Returns: LightAttenuationCoefficient [m^-1]
"""
function getLightAttenuation(day, settings::Dict{String, Any}, dynamicData::Dict{Int16, DayData})
    if ismissing(dynamicData[day].lightAttenuation)
        dynamicData[day].lightAttenuation = (
            settings["kdDev"] * (
                settings["maxKd"] -
                (settings["maxKd"] - settings["minKd"]) / 2 *
                (1 + cos((2 * pi / settings["yearlength"]) * (day - settings["kdDelay"])))
            )
        )
    end
    return (dynamicData[day].lightAttenuation) # [m^-1]
end



"""
    getWaterDepth(day; settings; dynamicData)

Calcuates waterdepth dependent on Waterlevel and LevelOfGrid

Arguments used from settings: (yearlength,maxW,minW,wDelay,levelCorrection)

Returns: Waterdepth [m]
"""
function getWaterDepth(day, LevelOfGrid, settings::Dict{String, Any}, dynamicData::Dict{Int16, DayData})
    WaterDepth = getWaterlevel(day, settings, dynamicData) - LevelOfGrid
    return (WaterDepth) #[m]
end



"""
    getReducedLightAttenuation(day, Biomass; settings; dynamicData)

The Effect of vegetation on light attenuation : Reduction of turbidity due to plant Biomass ;
independent of von Growthform
########## NOT USED!!!

Source: van Nes et al. (2003)

Arguments used from settings: yearlength,kdDev,maxKd,minKd,kdDelay, backgrKd,hTurbReduction,pTurbReduction

Returns:  lightAttenuCoefAdjusted #[m^-1]
"""
function getReducedLightAttenuation(day, Biomass, settings::Dict{String, Any}, dynamicData::Dict{Int16, DayData})
    lightAttenuCoef = getLightAttenuation(day, settings, dynamicData)
    lightAttenuCoefAdjusted =
        settings["backgrKd"] +
        (lightAttenuCoef - settings["backgrKd"]) *
        (settings["hTurbReduction"]^settings["pTurbReduction"]) / (
            Biomass^settings["pTurbReduction"] +
            settings["hTurbReduction"]^settings["pTurbReduction"]
        )
    return (lightAttenuCoefAdjusted) #[m^-1]
end



"""
    getBiomassAboveZ(distWaterSurface, height1, height2, waterdepth, biomass1, biomass2)

Returns share of Biomass above distinct distance from water surface
height 1 and height 2 for two competing species/growth forms
biomass1 and Biomass2 also

Source: -

Arguments used from settings: none

Result: BiomassAboveZ [g/m^2]
"""
function getBiomassAboveZ(distWaterSurface, height1, height2, waterdepth, biomass1, biomass2)
    if height1>0
        BiomassAboveZ_1 = ((height1 - (waterdepth - distWaterSurface)) / height1) * biomass1
    else
        BiomassAboveZ_1 =0
    end
    if BiomassAboveZ_1 <0
        BiomassAboveZ_1=0
    end
    if height2>0
        BiomassAboveZ_2 = ((height2 - (waterdepth - distWaterSurface)) / height2) * biomass2
    else
        BiomassAboveZ_2=0
    end
    if BiomassAboveZ_2 <0
        BiomassAboveZ_2=0
    end
    BiomassAboveZ = BiomassAboveZ_1 + BiomassAboveZ_2
    return (BiomassAboveZ) #[g/m^2]
end


#getBiomassAboveZ(1.0,1.5,0.5,2.0,5.0,1.0)

"""
    getEffectiveIrradianceHour(day,hour,distWaterSurface,Biomass1, Biomass2,height1, height2; settings; dynamicData)

Description

Source: van Nes et al. (2003)

Arguments from settings: parFactor, fracReflected, iDev, plantK, fracPeriphyton, latitude, maxI, minI, iDelay,
yearlength,kdDev, maxKd, minKd, kdDelay, backgrKd, hTurbReduction, pTurbReduction, LevelOfGrid,
maxW, minW, wDelay, levelCorrection

Result: lightPlantHour=effectiveIrradiance #[µE/m^2*s]
"""
function getEffectiveIrradianceHour(
    day,
    hour::Int64,
    distWaterSurface,
    Biomass1,
    Biomass2,
    height1,
    height2,
    LevelOfGrid,
    settings::Dict{String, Any},
    dynamicData::Dict{Int16, DayData}
)
    irrSurfHr = getSurfaceIrradianceHour(day, hour, settings, dynamicData) #(µE m^-2*s^-1)
    irrSubSurfHr =
        irrSurfHr *
        (1 - settings["parFactor"]) * #PAR radiation
        (1 - settings["fracReflected"]) * # Reflection at water surface
        settings["iDev"] # Deviation factor
    #lightAttenuCoef = getReducedLightAttenuation(day, (Biomass1+Biomass2), settings)
    lightAttenuCoef = getLightAttenuation(day, settings, dynamicData) #ohne feedback auf kd durch Pflanzen
    waterdepth = getWaterDepth(day,LevelOfGrid, settings, dynamicData)
    if height1>waterdepth
        height1=waterdepth
    end
    higherbiomass = getBiomassAboveZ(distWaterSurface, height1, height2, waterdepth, Biomass1, Biomass2)
    lightWater =
        irrSubSurfHr *
        exp(1)^(-lightAttenuCoef * distWaterSurface - settings["plantK"] * higherbiomass) # LAMBERT BEER # ÂµE/m^2*s # MÃ¶glichkeit im Exponenten: (absorptivity*c_H2O_pure*dist_water_surface))
    lightPlantHour = lightWater - (lightWater * settings["fracPeriphyton"]) ## µE/m^2*s
    return lightPlantHour #[µE/m^2*s]
end

#getEffectiveIrradianceHour(180, 5, 1.0, 0.0, 0.0, 0.0, 0.0, -2.0,settings, dynamicData)




"""
    getRespiration(day, settings, dynamicData)

Temperature dependence of maintenance respiration is formulated using a Q10 of 2

Source: van Nes et al. (2003)

Arguments from settings: resp20, q10

Result: (Respiration) #[g g^-1 d^-1]
"""
function getRespiration(day, settings::Dict{String, Any}, dynamicData::Dict{Int16, DayData}) #DAILY VALUE
    Temper = getTemperature(day, settings, dynamicData)
    Respiration = settings["resp20"] * settings["q10"]^((Temper - 20.0) / 10)
    return (Respiration) #[g g^-1 d^-1]
end


"""
Test new functions for alternative TempFactor

using Plots
temp=25.0
mPhotoTemp=20.0
bPhotoTemp=8.0

function getTempFactor(temp,mPhotoTemp=25,bPhotoTemp=5)
    TempFactor=exp(-((temp-mPhotoTemp)^ 2)/(2*bPhotoTemp^2))
    return (TempFactor)
end

getTempFactor(temp,mPhotoTemp,bPhotoTemp)
plot(getTempFactor,0,30)


function getTempFactorBroad(temp,minPhotoTemp=15,bPhotoTemp=2,
    maxPhotoTemp=25)
    if temp<minPhotoTemp
        TempFactor=exp(-((temp-minPhotoTemp)^ 2)/(2*bPhotoTemp^2))
    elseif temp>maxPhotoTemp
        TempFactor=exp(-((temp-maxPhotoTemp)^ 2)/(2*bPhotoTemp^2))
    else
        TempFactor=1.0
    end
    return (TempFactor)
end
plot(getTempFactorBroad,0,30)
"""

"""
    getPhotosynthesis(day,hour,distWaterSurf,Biomass1, Biomass2, height1, height2,settings,dynamicData)

Calculation of PS every hour dependent on light, temperature, dist (plant aging), [Carbonate, Nutrients]

Source: van Nes et al. (2003)

Arguments from settings: yearlength, maxW, minW, wDelay, levelCorrection, hPhotoDist, parFactor,
fracReflected, iDev, plantK, fracPeriphyton, latitude, maxI, minI, iDelay, kdDev, maxKd, minKd,
kdDelay, backgrKd, hTurbReduction, pTurbReduction, hPhotoLight, tempDev, maxTemp, minTemp,
tempDelay, sPhotoTemp, pPhotoTemp, hPhotoTemp, #bicarbonateConc, #hCarbonate, #pCarbonate,
#nutrientConc, #pNutrient, #hNutrient, pMax

Result: psHour [g / g * h]
"""
#Photosynthesis (Biomass brutto growth) (g g^-1 h^-1)
function getPhotosynthesis(
    day,
    hour::Int64,
    distFromPlantTop,
    Biomass1,
    Biomass2,
    height1,
    height2,
    LevelOfGrid,
    settings::Dict{String,Any},
    dynamicData::Dict{Int16, DayData}
)

    waterdepth = getWaterDepth(day, LevelOfGrid, settings, dynamicData)
    distWaterSurf = waterdepth - height1 + distFromPlantTop
    if height1 > waterdepth
        height1 = waterdepth
    end
    distFactor = settings["hPhotoDist"] / (settings["hPhotoDist"] + distFromPlantTop) #m

    lightPlantHour = getEffectiveIrradianceHour(
        day,
        hour,
        distWaterSurf,
        Biomass1,
        Biomass2,
        height1,
        height2,
        LevelOfGrid,
        settings,
        dynamicData
    )
    lightFactor = lightPlantHour / (lightPlantHour + settings["hPhotoLight"]) #ÂµE m^-2 s^-1); The default half-saturation constants (C aspera 14 yE m-2s-1; P pectinatus 52) are based on growth experiments

    #lightFactor_new = exp(-((lightPlantHour-mPhotoLight)^ 2)/(2*bPhotoLight^2))

    temp = getTemperature(day, settings, dynamicData)
    tempFactor =
        (settings["sPhotoTemp"] * (temp^settings["pPhotoTemp"])) /
        ((temp^settings["pPhotoTemp"]) + (settings["hPhotoTemp"]^settings["pPhotoTemp"])) #Â°C

    #tempFactor_new = exp(-((temp-mPhotoTemp)^ 2)/(2*bPhotoTemp^2))

    #bicarbFactor = bicarbonateConc ^ pCarbonate / (bicarbonateConc ^ pCarbonate + hCarbonate ^ pCarbonate) # C.aspera hCarbonate=30 mg/l; P.pectinatus hCarbonate=60 mg/l

    #nutrientConc = reduceNutrientConcentration((Biomass1 + Biomass2), settings)
    nutrientConc=settings["maxNutrient"]
    nutrientFactor =
        (nutrientConc^settings["pNutrient"]) /
        (nutrientConc^settings["pNutrient"] + settings["hNutrient"]^settings["pNutrient"])

    psHour = settings["pMax"] * lightFactor * tempFactor * distFactor * nutrientFactor #* bicarbFactor # #(g g^-1 h^-1)

    return (psHour) ##[g / g * h]
end


"""
    getPhotosynthesisPLANTDay(day, height, Biomass; settings, dynamicData)

Calculation of daily PS

Source: van Nes et al. (2003)

Arguments used from settings: latitude,LevelOfGrid,yearlength,maxW, minW, wDelay, levelCorrection,
parFactor, fracReflected, iDev, plantK, fracPeriphyton, maxI, minI, iDelay, kdDev, maxKd, minKd,
kdDelay, backgrKd, hTurbReduction, pTurbReduction, hPhotoDist, hPhotoLight, tempDev,
maxTemp, minTemp, tempDelay, sPhotoTemp, pPhotoTemp, hPhotoTemp, pMax

Returns: PS daily [g / g * d]
"""
#using QuadGK
#using HCubature
function getPhotosynthesisPLANTDay(
    day,
    height1,
    height2,
    Biomass1,
    Biomass2,
    LevelOfGrid,
    settings::Dict{String,Any},
    dynamicData::Dict{Int16, DayData}
)

    daylength = getDaylength(day, settings, dynamicData)
    waterdepth = getWaterDepth((day), LevelOfGrid, settings, dynamicData)
    distPlantTopFromSurf = waterdepth - height1
    if height1 > waterdepth
        height1 = waterdepth
    end
    PS = 0
    if Biomass1 > 0.0
        for i = 1:floor(daylength) #Rundet ab # Loop über alle Stunden
            i = convert(Int64, i)
            for j in 0:0.1:1
                PS =
                    PS + getPhotosynthesis(
                            day,i,
                            j* height1,
                            Biomass1,Biomass2,
                            height1,height2,
                            LevelOfGrid,settings,dynamicData,
                        )*1/11 #because it is calculated in 11 steps, to calc the mean
            end
        end
    else
        PS = 0
    end
    return PS
end



#getPhotosynthesis(180,5,1.0,100.0,0.0,1.0,1.5,-2.0,settings)


"""
    growHeight(biomass; settings)

Height growth of plants

Source:

Arguments used from settings: maxWeightLenRatio

Returns: height [m]
"""
function growHeight(indBiomass::Float64, settings::Dict{String, Any})
    if indBiomass > 0
        height2 = indBiomass / settings["maxWeightLenRatio"]
    else
        height2 = 0
    end
    return height2
end


"""
    getDailyGrowth(seeds, biomass1, allocatedBiomass1, dailyPS, dailyRES, settings)

Calcualtion of daily growth

Source: van Nes et al. (2003)

Arguments used from settings: cTuber, rootShootRatio

Returns: daily biomass increase [g]
"""
function getDailyGrowth(
    seeds::Float64,
    biomass1::Float64,
    allocatedBiomass1::Float64,
    dailyPS::Float64,
    dailyRES::Float64,
    settings::Dict{String,Any},
)
    dailyGrowth =
        seeds * settings["cTuber"] + #Growth from seedBiomass
        (
            ((1 - settings["rootShootRatio"]) * biomass1 - allocatedBiomass1) * dailyPS - #GrossProduction : Growth from sprout
            biomass1 * dailyRES #Respiration
        )
    #if dailyGrowth > 0 # No negative growth allowed
    #    dailyGrowth = dailyGrowth
    #else
    #    dailyGrowth = 0
    #end
    return dailyGrowth
end


"""
    getNumberOfSeedsProducedByOnePlant(day, settings)

Description
Not used in that form in the code

Source: van Nes et al. (2003)

Arguments used from settings: seedFraction,seedBiomass

Returns: seedNumber [N]
"""
function getNumberOfSeedsProducedByOnePlant(Biomass, settings::Dict{String, Any})
    seedNumber = settings["seedFraction"] * Biomass / settings["seedBiomass"]
    #return round(seedNumber)
    return seedNumber
end

#getNumberOfSeedsProducedByOnePlant(0.4, settings)


"""
    getNumberOfSeeds(seedBiomass; settings)

Calculates number of Seeds by single seed biomass

Source: van Nes et al. (2003)

Arguments used from settings:

Returns: []
"""
function getNumberOfSeeds(seedBiomass, settings::Dict{String, Any})
    if settings["seedBiomass"]== 0
        seedNumber =0
    else
        seedNumber = seedBiomass / settings["seedBiomass"]
    end

    #return round(seedNumber)
    return (seedNumber)
end

"""
    getNumberOfTubers(tubersBiomass; settings)

Calculates number of Seeds by single seed biomass

Source: van Nes et al. (2003)

Arguments used from settings:

Returns: []
"""
function getNumberOfTubers(tubersBiomass, settings::Dict{String, Any})
    if settings["tuberBiomass"]==0
        tubersNumber=0
    else
        tubersNumber = tubersBiomass / settings["tuberBiomass"]
    end
    #return round(tubersNumber)
    return (tubersNumber)
end


"""
    getIndividualWeight(Biomass, Number)

Returns inidividual Weight of each plant represented by the Super-Individuum

Source: van Nes et al. (2003)

Arguments used from settings: none

Returns: indWeight [g]
"""
function getIndividualWeight(Biomass, Number)
    indWeight = Biomass / Number
    return indWeight
end



"""
    dieThinning(number, individualWeight)

Mortality due to competition at high plant denisties

Source: van Nes et al. (2003)

Arguments used from settings: none

Returns: numberAdjusted, individualWeightADJ []
"""
function dieThinning(number, individualWeight, settings::Dict{String, Any})
    numberAdjusted = (settings["cThinning"] / individualWeight)^(2 / 3)
    individualWeightADJ = (number / numberAdjusted) * individualWeight
    #return (round(numberAdjusted), individualWeightADJ)
    return (numberAdjusted, individualWeightADJ)
end

"""
    dieWaves(day,LevelOfGrid,settings, dynamicData)

Mortality due to wave damage; loss in number of plants untill reached water surface; Adult plants only lose weight

Source: van Nes et al. (2003)

Arguments used from settings: maxWaveMort,hWaveMort,pWaveMort

Returns: wave mortality [d^-1]
"""
function dieWaves(day, LevelOfGrid, settings::Dict{String, Any}, dynamicData::Dict{Int16, DayData})
    waterdepth = getWaterDepth(day, LevelOfGrid, settings, dynamicData)
    waveMortality =
        settings["maxWaveMort"] * (settings["hWaveMort"]^settings["pWaveMort"]) /
        ((settings["hWaveMort"]^settings["pWaveMort"]) + (waterdepth^settings["pWaveMort"]))
    return waveMortality
end




"""
    killWithProbability(Mort, N1)

Killing number of Plants by using a random number from Poisson distribution

Source: van Nes et al. (2003)

Arguments used from settings: none

Returns: Number of plants reduced
"""
function killWithProbability(Mort, N1)
    N2=0
    for i in 1:N1
        N2 = N2+rand(Binomial(1,1-Mort))[1]
    end
    return(N2)
end


"""
    killN(Mort, N1)

Killing number of Plants

Source: van Nes et al. (2003)

Arguments used from settings: Mort

Returns: Number of plants reduced
"""
function killN(Mort, N1)
    #N2=0
    N2= (1-Mort)*N1
    return(N2)
end
