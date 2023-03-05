#Macrophyte Growth Model (MGM)
#
#Anne Lewerentz <anne.lewerentz@uni-wuerzburg.de>
#(c) 2021-2022, licensed under the terms of the MIT license
#
#Contains functions that are necessary to run the simulation

"""
    simulate(LevelOfGrid; settings)

Simulation function for growth of macrophytes in one depth

Code following the description of Charisma by van Nes et al. (2003)

Arguments used from settings:
years, yearlength, germinationDay, heightMax, rootShootRatio, BackgroundMort, resp20,
q10, latitude, maxW, minW, wDelay, levelCorrection, parFactor, fracReflected,
iDev, plantK, fracPeriphyton, maxI, minI, iDelay, kdDev, maxKd, minKd, kdDelay,
hPhotoDist, hPhotoLight, tempDev, maxTemp, minTemp, tempDelay, sPhotoTemp,
pPhotoTemp, hPhotoTemp, pMax, maxAge, maxWeightLenRatio, seedInitialBiomass,
seedFraction, cTuber, seedGermination, seedBiomass, seedsStartAge, seedsEndAge,
reproDay, SeedMortality, spreadFrac, backgrKd, hTurbReduction, pTurbReduction,
thinning

Returns: [superInd]
"""
function simulate(LevelOfGrid, settings::Dict{String, Any}, dynamicData::Dict{Int16, DayData})
    #simlog("Starting simulation.", settings)

    #Initialisation
    seeds = zeros(Float64, settings["yearlength"], 3, settings["years"]) #SeedBiomass, SeedNumber, SeedsGerminatingBiomass
    tubers = zeros(Float64, settings["yearlength"], 3, settings["years"]) #tubersBiomass, tubersNumber
    superIndSeeds = zeros(Float64, settings["yearlength"], 6, settings["years"]) #Biomass, Number, indWeight, Height, allocatedSeedBiomass, allocatedTurionsBiomass
    superIndTubers = zeros(Float64, settings["yearlength"], 6, settings["years"]) #Biomass, Number, indWeight, Height, allocatedSeedBiomass, allocatedTurionsBiomass
    growthSeeds = zeros(Float64, settings["yearlength"], 3, settings["years"]) #dailyPS, dailyRES, dailyGrowth
    growthTubers = zeros(Float64, settings["yearlength"], 3, settings["years"]) #dailyPS, dailyRES, dailyGrowth
    #lightAttenuation = zeros(Float64, settings["yearlength"], 1, settings["years"]) #reducedlightAttenuation

    #Initialize variables from settings
    germinationDay = settings["germinationDay"]


    #LOOP OVER YEARS
    for y = 1:settings["years"] #Loop over years
        #println(y)

        # INITIALISATION
        if y == 1 #Initialize in first year SeedBiomass & SeedNumber
            seeds[1, 1, 1] = settings["seedInitialBiomass"] #initial SeedBiomass
            tubers[1, 1, 1] = settings["tuberInitialBiomass"] #initial TurionBiomass
            seeds[1, 2, 1] = getNumberOfSeeds(seeds[1, 1, 1],settings) #initial SeedNumber
            tubers[1, 2, 1] = getNumberOfTubers(tubers[1, 1, 1],settings) #initial tubersNumber
        else
            seeds[1, 1, y] = seeds[settings["yearlength"], 1, y-1] # get biomass of last day of last year
            seeds[1, 2, y] = getNumberOfSeeds(seeds[1, 1, y],settings)
            tubers[1, 1, y] = tubers[settings["yearlength"], 1, y-1] # get biomass of last day of last year
            tubers[1, 2, y] = getNumberOfTubers(tubers[1, 1, y],settings)
        end

        #BREAK CRITERIA to shorten runtime
        if seeds[1, 1, y]==0.0 && tubers[1, 1, y]==0.0
              #println("Break")
              break
        #else println(seeds[1, 1, y]+tubers[1, 1, y])

        #TODO if seed/tuber N stays the same - break

        end

        # LOOP OVER DAYS
        for d = 2:settings["yearlength"]
            WaterDepth = getWaterDepth(d, LevelOfGrid, settings, dynamicData)

            ########################################################################
            ##SEEDS: NO GROWTH UNTILL GERMINATION
            if superIndSeeds[d-1,1,y] == 0
                seeds[d, 1, y] = seeds[d-1, 1, y] - (seeds[d-1, 2, y] * settings["seedMortality"] *settings["seedBiomass"])
                seeds[d, 2, y] = getNumberOfSeeds(seeds[d, 1, y], settings) #SeedNumber

                if seeds[d, 2, y] < 1
                    seeds[d, 1, y] =0
                end

                #SEED GERMINATION
                if d == germinationDay
                    #SeedsGerminationBiomass determination
                    seeds[germinationDay, 3, y] = seeds[(germinationDay-1), 1, y] * settings["seedGermination"] #20% of the SeedsBiomass are transformed to SeedsGerminatingBiomass
                    seeds[germinationDay, 1, y] = #SeedBiomass update
                        seeds[germinationDay-1, 1, y] -
                        seeds[germinationDay, 3, y] -
                        seeds[germinationDay-1, 2, y] * settings["seedMortality"] *settings["seedBiomass"]#Remaining SeedsBiomass

                    seeds[germinationDay, 2, y] =#Update Number of left seeds
                        getNumberOfSeeds(seeds[germinationDay, 1, y],settings)

                    if seeds[germinationDay, 2, y] < 1
                        seeds[germinationDay, 1, y] =0
                    end

                    superIndSeeds[germinationDay, 2, y] = #Number of Germinated Individuals
                        getNumberOfSeeds(seeds[germinationDay, 3, y],settings)

                    superIndSeeds[germinationDay, 1, y] = #Calcualtion of Starting Plant Biomass
                        seeds[germinationDay, 3, y] * settings["cTuber"]

                    superIndSeeds[germinationDay, 3, y] = getIndividualWeight(
                        superIndSeeds[germinationDay, 1, y],
                        superIndSeeds[germinationDay, 2, y],
                    ) #Starting individualWeight

                    if superIndSeeds[germinationDay, 3, y] < 0
                        superIndSeeds[germinationDay, 3, y] =0
                    end

                    # Thinning, optional
                    if settings["thinning"] == true
                        thin =
                            dieThinning(superIndSeeds[germinationDay, 2, y],
                                superIndSeeds[germinationDay, 3, y], settings) #Adapts number of individuals [/m^2]& individual weight
                        #Rule to not get more individuals out of thinning
                        if (thin[1] < superIndSeeds[germinationDay, 2, y])
                            superIndSeeds[germinationDay, 2, y] = thin[1]
                            superIndSeeds[germinationDay, 3, y] = thin[2]
                        end
                    end

                    superIndSeeds[germinationDay, 4, y] =
                        growHeight(superIndSeeds[germinationDay, 3, y], settings)

                    if superIndSeeds[germinationDay, 2, y] == 0 #< 1
                        superIndSeeds[germinationDay, 2, y] = 0
                        superIndSeeds[germinationDay, 1, y] = 0
                        superIndSeeds[germinationDay, 3, y] = 0
                        superIndSeeds[germinationDay, 4, y] = 0
                    end
                end
            end


            ########################################################################
            #TUBERS: NO GROWTH UNTILL GERMINATION
            if superIndTubers[d-1,1,y] == 0
                tubers[d, 1, y] = tubers[d-1, 1, y] - (tubers[d-1, 2, y] * settings["tuberMortality"] * settings["tuberBiomass"]) #minus SeedMortality #SeedBiomass
                tubers[d, 2, y] = getNumberOfTubers(tubers[d, 1, y], settings) #SeedNumber

                if tubers[d, 2, y] < 1
                    tubers[d, 1, y] =0
                end
                #TUBERS GERMINATION
                if d == settings["tuberGerminationDay"]
                    tubers[settings["tuberGerminationDay"], 3, y] = #SeedsGerminationBiomass determination
                        tubers[(settings["tuberGerminationDay"]-1), 1, y] * settings["tuberGermination"] #20% of the SeedsBiomass are transformed to SeedsGerminatingBiomass
                    tubers[settings["tuberGerminationDay"], 1, y] = #SeedBiomass update
                        tubers[settings["tuberGerminationDay"]-1, 1, y] -
                        tubers[settings["tuberGerminationDay"], 3, y] -
                        tubers[settings["tuberGerminationDay"]-1, 2, y] * settings["tuberMortality"] * settings["tuberBiomass"]#Remaining tubersBiomass


                    tubers[settings["tuberGerminationDay"], 2, y] =#Update Number of left tubers
                        getNumberOfTubers(tubers[settings["tuberGerminationDay"], 1, y],settings)

                        if tubers[d, 2, y] < 1
                            tubers[d, 1, y] =0
                        end

                    superIndTubers[settings["tuberGerminationDay"], 2, y] = #Number of Germinated Individuals
                        getNumberOfTubers(tubers[settings["tuberGerminationDay"], 3, y],settings)

                    superIndTubers[settings["tuberGerminationDay"], 1, y] = #Calcualtion of Starting Plant Biomass
                        tubers[settings["tuberGerminationDay"], 3, y] * settings["cTuber"]

                    superIndTubers[settings["tuberGerminationDay"], 3, y] = getIndividualWeight(
                        superIndTubers[settings["tuberGerminationDay"], 1, y],
                        superIndTubers[settings["tuberGerminationDay"], 2, y],
                    ) #Starting individualWeight

                    if superIndTubers[settings["tuberGerminationDay"], 3, y] < 0
                        superIndTubers[settings["tuberGerminationDay"], 3, y] =0
                    end

                    # Thinning, optional
                    if settings["thinning"] == true
                        thin =
                            dieThinning(superIndTubers[settings["tuberGerminationDay"], 2, y],
                                superIndTubers[settings["tuberGerminationDay"], 3, y], settings) #Adapts number of individuals [/m^2]& individual weight
                        #Rule to not get more individuals out of thinning
                        if (thin[1] < superIndTubers[settings["tuberGerminationDay"], 2, y])
                            superIndTubers[settings["tuberGerminationDay"], 2, y] = thin[1]
                            superIndTubers[settings["tuberGerminationDay"], 3, y] = thin[2]
                        end
                    end

                    superIndTubers[settings["tuberGerminationDay"], 4, y] =
                        growHeight(superIndTubers[settings["tuberGerminationDay"], 3, y], settings)

                    if superIndTubers[settings["tuberGerminationDay"], 2, y] ==0 #< 1
                        superIndTubers[settings["tuberGerminationDay"], 2, y] = 0
                        superIndTubers[settings["tuberGerminationDay"], 1, y] = 0
                        superIndTubers[settings["tuberGerminationDay"], 3, y] = 0
                        superIndTubers[settings["tuberGerminationDay"], 4, y] = 0
                    end #no half individuals
                end
            end


            ########################################################################
            #GROWTH from SEEDS
            if superIndSeeds[d-1,1,y] > 0 && superIndTubers[d-1,1,y] == 0
                seeds[d, 1, y] = seeds[d-1, 1, y] - (seeds[d-1, 2, y] * settings["seedMortality"] * settings["seedBiomass"]) #minus SeedMortality #SeedBiomass
                seeds[d, 3, y] = (1 - settings["cTuber"]) * seeds[d-1, 3, y] #Reduction of allocatedBiomass untill it is used
                seeds[d, 2, y] = getNumberOfSeeds(seeds[d, 1, y],settings) #SeedNumber
                if seeds[d, 2, y] < 1
                    seeds[d, 1, y] =0
                end
                superIndSeeds[d, 2, y] = superIndSeeds[d-1, 2, y] #

                #GROWTH
                if superIndSeeds[d-1, 4, y] > WaterDepth
                    superIndSeeds[d-1, 4, y] = WaterDepth
                end

                growthSeeds[d, 2, y] = getRespiration(d, settings, dynamicData) #[g / g*d]

                growthSeeds[d, 1, y] = getPhotosynthesisPLANTDay( #[g / g*d]
                    d,
                    superIndSeeds[d-1, 4, y], #height yesterday
                    superIndTubers[d-1, 4, y],
                    ((1 - settings["rootShootRatio"]) * (superIndSeeds[d-1, 1, y]-superIndSeeds[d-1, 5, y]-superIndSeeds[d-1, 6, y])),
                    ((1 - settings["rootShootRatio"]) * (superIndTubers[d-1, 1, y]-superIndTubers[d-1, 5, y]-superIndTubers[d-1, 6, y])),
                    LevelOfGrid,
                    settings,
                    dynamicData,
                )

                growthSeeds[d, 3, y] = getDailyGrowth(
                    seeds[d, 3, y], #SeedGerminating Biomass
                    superIndSeeds[d-1, 1, y], #Yesterdays Biomass
                    (superIndSeeds[d-1, 5, y]), #Yesterdays allocatedBiomass from for seeds & turions
                    growthSeeds[d, 1, y], # Todays PS
                    growthSeeds[d, 2, y], # Todays Resp
                    settings,
                )

                ##Consequence of negative growth Seeds
                #if growthSeeds[d, 3, y] < 0 &&
                #   superIndSeeds[d-1, 3, y] > 0  #Check if indWeight>0
                #   Mort = (-growthSeeds[d, 3, y] / superIndSeeds[d-1, 3, y])
                #   if Mort >1
            #           Mort = 1
            #       end
            #        superIndSeeds[d, 2, y] = #Loss in number of Plants
            #            (1-Mort)*superIndSeeds[d,2,y]
            #        superIndSeeds[d, 1, y] = superIndSeeds[d-1, 1, y] #Biomass stays the same. Makes sense?
            #    else
            #        superIndSeeds[d, 1, y] = superIndSeeds[d-1, 1, y] + growthSeeds[d, 3, y] #Biomass
            #    end

                #Growth changed from original CHARIMA
                superIndSeeds[d, 1, y] = superIndSeeds[d-1, 1, y] + growthSeeds[d, 3, y] #If growth is negative: Biomass gets reduced

                #Mortality (N_Weight_Mortality)
                #Mort = dieWaves(d,LevelOfGrid,settings) + settings["BackgroundMort"] #+Herbivory
                #if superIndSeeds[d-1, 4, y] < settings["heightMax"] #Check if plant is not yet adult
                #    superIndSeeds[d, 2, y] = (1-Mort) * superIndSeeds[d, 2, y] # Loss in number of plants
                #else #Plant is adult -> lost of weight
                #    superIndSeeds[d, 1, y] = (1-Mort)*superIndSeeds[d, 1, y]
                #end

                #Mortality changed from original CHARISMA
                Mort = dieWaves(d,LevelOfGrid,settings,dynamicData) + settings["BackgroundMort"] #+Herbivory
                superIndSeeds[d, 1, y] = (1-Mort)*superIndSeeds[d, 1, y] # Loss in total Biomass
                superIndSeeds[d, 2, y] = (1-Mort) * superIndSeeds[d, 2, y] # Loss in number of plants

                # Calculate individual Weight
                superIndSeeds[d, 3, y] = getIndividualWeight(superIndSeeds[d, 1, y], superIndSeeds[d, 2, y]) #individualWeight = Biomass / Number

                if superIndSeeds[d, 3, y] < 0
                    superIndSeeds[d, 3, y] =0
                end

                #Thinning, optional
                if settings["thinning"] == true
                    thin = dieThinning(superIndSeeds[d, 2, y], superIndSeeds[d, 3, y], settings) #Adapts number of individuals [/m^2]& individual weight
                    if (thin[1] < superIndSeeds[d, 2, y]) #&& (Thinning[2] > 0)
                        superIndSeeds[d, 2, y] = thin[1] #N
                        superIndSeeds[d, 3, y] = thin[2] #indWeight
                    end
                end

                #Die-off if N==0
                if superIndSeeds[d, 2, y] ==0 #< 1
                    superIndSeeds[d, 2, y] = 0
                    superIndSeeds[d, 1, y] = 0
                    superIndSeeds[d, 3, y] = 0
                    superIndSeeds[d, 4, y] = 0
                end #no half individuals

                #Height calc
                superIndSeeds[d, 4, y] = growHeight(superIndSeeds[d, 3, y],settings)
                if superIndSeeds[d, 4, y] >= settings["heightMax"]
                    superIndSeeds[d, 4, y] = settings["heightMax"]
                end
                if superIndSeeds[d, 4, y] >= WaterDepth
                    superIndSeeds[d, 4, y] = WaterDepth
                end

                #ALLOCATION OF BIOMASS FOR SEED PRODUCTION
                if d > (germinationDay + settings["seedsStartAge"]) &&
                   d < (germinationDay + settings["seedsEndAge"]) #Age=Age of plant
                    superIndSeeds[d, 5, y] =
                        settings["seedFraction"] * ((d-germinationDay - settings["seedsStartAge"]) /
                        (settings["seedsEndAge"] - settings["seedsStartAge"])) * superIndSeeds[d, 1, y] # Copied from code!
                end
                if d >= (germinationDay + settings["seedsEndAge"]) &&
                   d <= settings["reproDay"]
                    superIndSeeds[d, 5, y] = settings["seedFraction"] * superIndSeeds[d, 1, y]  #allocatedBiomass - Fraction stays the same
                end
                #ALLOCATION OF BIOMASS FOR TUBER PRODUCTION
                if d > (settings["tuberGerminationDay"] + settings["tuberStartAge"]) &&
                   d < (settings["tuberGerminationDay"] + settings["tuberEndAge"]) #Age=Age of plant
                    superIndSeeds[d, 6, y] =
                        settings["tuberFraction"] * ((d-settings["tuberGerminationDay"] - settings["tuberStartAge"]) /
                        (settings["tuberEndAge"] - settings["tuberStartAge"])) * superIndSeeds[d, 1, y] # Copied from code!
                end
                if d >= (settings["tuberGerminationDay"] + settings["tuberEndAge"]) &&
                   d <= settings["reproDay"]
                    superIndSeeds[d, 6, y] = settings["tuberFraction"] * superIndSeeds[d, 1, y]  #allocatedBiomass - Fraction stays the same
                end

                if d == settings["reproDay"]
                    seeds[d, 1, y] =
                        seeds[d-1, 1, y] + superIndSeeds[d, 5, y] -
                        seeds[d-1, 1, y] * settings["seedMortality"]
                    tubers[d, 1, y] =
                            tubers[d-1, 1, y] + superIndSeeds[d, 6, y] -
                            tubers[d-1, 1, y] * settings["tuberMortality"]
                    superIndSeeds[d, 1, y] = superIndSeeds[d,1,y]-superIndSeeds[d, 5, y]-superIndSeeds[d, 6, y] #Loss of total Biomass as seeds are distributed
                    superIndSeeds[d, 5, y] = 0 #Allocated Biomass is lost
                    superIndSeeds[d, 6, y] = 0 #Allocated Biomass is lost
                end

                if d == (germinationDay + settings["maxAge"])
                    superIndSeeds[d, 1, y] = 0
                    superIndSeeds[d, 2, y] = 0
                    superIndSeeds[d, 3, y] = 0
                    superIndSeeds[d, 4, y] = 0
                    superIndSeeds[d, 5, y] = 0
                end
            end


            ########################################################################
            #GROWTH just from TUBERS
            if superIndSeeds[d-1,1,y] == 0 && superIndTubers[d-1,1,y] > 0
                tubers[d, 1, y] = tubers[d-1, 1, y] - (tubers[d-1, 2, y] * settings["tuberMortality"] * settings["tuberBiomass"]) #minus SeedMortality #SeedBiomass
                tubers[d, 3, y] = (1 - settings["cTuber"]) * tubers[d-1, 3, y] #Reduction of allocatedBiomass untill it is used
                tubers[d, 2, y] = getNumberOfTubers(tubers[d, 1, y], settings) #SeedNumber
                if tubers[d, 2, y] < 1
                    tubers[d, 1, y] =0
                end
                superIndTubers[d, 2, y] = superIndTubers[d-1, 2, y] #TODO do I need this line here?

                #GROWTH
                if superIndTubers[d-1, 4, y] > WaterDepth
                    superIndTubers[d-1, 4, y] = WaterDepth
                end

                growthTubers[d, 2, y] = getRespiration(d, settings, dynamicData) #[g / g*d]

                growthTubers[d, 1, y] = getPhotosynthesisPLANTDay( #[g / g*d]
                    d,
                    superIndTubers[d-1, 4, y],
                    superIndSeeds[d-1, 4, y],
                    ((1 - settings["rootShootRatio"]) * (superIndTubers[d-1, 1, y]-superIndTubers[d-1, 5, y]-superIndTubers[d-1, 6, y])),
                    ((1 - settings["rootShootRatio"]) * (superIndSeeds[d-1, 1, y]-superIndSeeds[d-1, 5, y]-superIndSeeds[d-1, 6, y])),
                    LevelOfGrid,
                    settings,
                    dynamicData,
                )

                growthTubers[d, 3, y] = getDailyGrowth(
                    tubers[d, 3, y], #SeedGerminating Biomass
                    superIndTubers[d-1, 1, y], #Yesterdays Biomass
                    (superIndTubers[d-1, 5, y]), #Yesterdays allocatedBiomass
                    growthTubers[d, 1, y], # Todays PS
                    growthTubers[d, 2, y], # Todays Resp
                    settings,
                )


                ##Consequence of negative growth TUBERS
                #if growthTubers[d, 3, y] < 0 &&
                #   superIndTubers[d-1, 3, y] > 0  #Check if indWeight>0
                #   Mort = (-growthTubers[d, 3, y] / superIndTubers[d-1, 3, y])
                #   if Mort >1
                #       Mort = 1
                #   end
                #   superIndTubers[d, 2, y] = (1-Mort)*superIndTubers[d, 2, y]
                #   superIndTubers[d, 1, y] = superIndTubers[d-1, 1, y] #Biomass stays the same. Makes sense?
                #else
                #    superIndTubers[d, 1, y] = superIndTubers[d-1, 1, y] + growthTubers[d, 3, y] #Biomass
                #end
                superIndTubers[d, 1, y] = superIndTubers[d-1, 1, y] + growthTubers[d, 3, y] #Biomass



                #Mortality (N_Weight_Mortality)
                Mort = dieWaves(d,LevelOfGrid,settings,dynamicData) + settings["BackgroundMort"] #+Herbivory
                #if superIndTubers[d-1, 4, y] < settings["heightMax"] #Check if plant is not yet adult
                #    superIndTubers[d, 2, y] = (1-Mort) * superIndTubers[d, 2, y] # Loss in number of plants
                #else #Plant is adult -> lost of weight
                #    superIndTubers[d, 1, y] = (1-Mort) * superIndTubers[d, 1, y] #Update total Biomass
                #end
                superIndTubers[d, 2, y] = (1-Mort) * superIndTubers[d, 2, y] # Loss in number of plants
                superIndTubers[d, 1, y] = (1-Mort) * superIndTubers[d, 1, y] #Update total Biomass

                superIndTubers[d, 3, y] = getIndividualWeight(superIndTubers[d, 1, y], superIndTubers[d, 2, y]) #individualWeight = Biomass / Number
                if superIndTubers[d, 3, y] < 0
                    superIndTubers[d, 3, y] =0
                end
                #Thinning, optional
                if settings["thinning"] == true
                    thin = dieThinning(superIndTubers[d, 2, y], superIndTubers[d, 3, y], settings) #Adapts number of individuals [/m^2]& individual weight
                    if (thin[1] < superIndTubers[d, 2, y]) #&& (Thinning[2] > 0)
                        superIndTubers[d, 2, y] = thin[1] #N
                        superIndTubers[d, 3, y] = thin[2] #indWeight
                    end
                end

                #Die-off if N<1
                if superIndTubers[d, 2, y] ==0 #< 1
                    superIndTubers[d, 2, y] = 0
                    superIndTubers[d, 1, y] = 0
                    superIndTubers[d, 3, y] = 0
                    superIndTubers[d, 4, y] = 0
                end #no half individuals

                #Height calc
                superIndTubers[d, 4, y] = growHeight(superIndTubers[d, 3, y],settings)
                if superIndTubers[d, 4, y] >= settings["heightMax"]
                    superIndTubers[d, 4, y] = settings["heightMax"]
                end
                if superIndTubers[d, 4, y] >= WaterDepth
                    superIndTubers[d, 4, y] = WaterDepth
                end

                #ALLOCATION OF BIOMASS FOR TUBER PRODUCTION
                if d > (settings["tuberGerminationDay"] + settings["tuberStartAge"]) &&
                   d < (settings["tuberGerminationDay"] + settings["tuberEndAge"]) #Age=Age of plant
                    superIndTubers[d, 6, y] =
                        settings["tuberFraction"] * ((d-settings["tuberGerminationDay"] - settings["tuberStartAge"]) /
                        (settings["tuberEndAge"] - settings["tuberStartAge"])) * superIndTubers[d, 1, y] # Copied from code!
                end
                if d >= (settings["tuberGerminationDay"] + settings["tuberEndAge"]) &&
                   d <= settings["reproDay"]
                    superIndTubers[d, 6, y] = settings["tuberFraction"] * superIndTubers[d, 1, y]  #allocatedBiomass - Fraction stays the same
                end
                #ALLOCATION OF BIOMASS FOR SEED PRODUCTION
                if d > (germinationDay + settings["seedsStartAge"]) &&
                   d < (germinationDay + settings["seedsEndAge"]) #Age=Age of plant
                    superIndTubers[d, 5, y] =
                        settings["seedFraction"] * ((d-germinationDay - settings["seedsStartAge"]) /
                        (settings["seedsEndAge"] - settings["seedsStartAge"])) * superIndTubers[d, 1, y] # Copied from code!

                end
                if d >= (germinationDay + settings["seedsEndAge"]) &&
                   d <= settings["reproDay"]
                    superIndTubers[d, 5, y] = settings["seedFraction"] * superIndTubers[d, 1, y]  #allocatedBiomass - Fraction stays the same

                end

                if d == settings["reproDay"]
                    tubers[d, 1, y] =
                        tubers[d-1, 1, y] + superIndTubers[d, 6, y] -
                        tubers[d-1, 1, y] * settings["tuberMortality"]
                    seeds[d,1,y]= seeds[d-1, 1, y] + superIndTubers[d, 5, y] -
                        seeds[d-1, 1, y] * settings["seedMortality"]
                    superIndTubers[d, 1, y] = superIndTubers[d,1,y]-superIndTubers[d, 5, y]-superIndTubers[d, 6, y] #Loss of total Biomass as Tubers are distributed
                    superIndTubers[d, 5, y] = 0 #Allocated Biomass is lost
                    superIndTubers[d, 6, y] = 0 #Allocated Biomass is lost
                end

                if d == (settings["tuberGerminationDay"] + settings["maxAge"])
                    superIndTubers[d, 1, y] = 0
                    superIndTubers[d, 2, y] = 0
                    superIndTubers[d, 3, y] = 0
                    superIndTubers[d, 4, y] = 0
                    superIndTubers[d, 5, y] = 0
                end
            end





            ########################################################################
            #GROWTH BOTH
            if superIndSeeds[d-1,1,y] > 0 && superIndTubers[d-1,1,y] > 0

                #Seedbank
                seeds[d, 1, y] = seeds[d-1, 1, y] - (seeds[d-1, 2, y] * settings["seedMortality"] * settings["seedBiomass"]) #minus SeedMortality #SeedBiomass
                seeds[d, 3, y] = (1 - settings["cTuber"]) * seeds[d-1, 3, y] #Reduction of allocatedBiomass untill it is used
                seeds[d, 2, y] = getNumberOfSeeds(seeds[d, 1, y],settings) #SeedNumber
                if seeds[d, 2, y] < 1
                    seeds[d, 1, y] =0
                end
                #Tuberbank? They should all "germinate" (TODO check)
                tubers[d, 1, y] = tubers[d-1, 1, y] - (tubers[d-1, 2, y] * settings["tuberMortality"] * settings["tuberBiomass"]) #minus SeedMortality #SeedBiomass
                tubers[d, 3, y] = (1 - settings["cTuber"]) * tubers[d-1, 3, y] #Reduction of allocatedBiomass untill it is used
                tubers[d, 2, y] = getNumberOfTubers(tubers[d, 1, y],settings) #SeedNumber
                if tubers[d, 2, y] < 1
                    tubers[d, 1, y] =0
                end
                #N individuals (TODO do I need this line here?)
                superIndSeeds[d, 2, y] = superIndSeeds[d-1, 2, y]
                superIndTubers[d, 2, y] = superIndTubers[d-1, 2, y]

                # Hight growth limitation in case of falling water level
                if superIndSeeds[d-1, 4, y] > WaterDepth
                    superIndSeeds[d-1, 4, y] = WaterDepth
                end
                if superIndTubers[d-1, 4, y] > WaterDepth
                    superIndTubers[d-1, 4, y] = WaterDepth
                end

                # GROWTH calculation for 2 superInds: Order: first Seeds,
                growthSeeds[d, 2, y] = getRespiration(d, settings, dynamicData) #[g / g*d]
                growthTubers[d, 2, y] = getRespiration(d, settings, dynamicData) #[g / g*d]

                growthSeeds[d, 1, y] = getPhotosynthesisPLANTDay( #[g / g*d]
                    d,
                    superIndSeeds[d-1, 4, y],
                    superIndTubers[d-1, 4, y],
                    ((1 - settings["rootShootRatio"]) * (superIndSeeds[d-1, 1, y]-superIndSeeds[d-1, 5, y]-superIndSeeds[d-1, 6, y])),
                    ((1 - settings["rootShootRatio"]) * (superIndTubers[d-1, 1, y]-superIndTubers[d-1, 5, y]-superIndTubers[d-1, 6, y])),
                    LevelOfGrid,
                    settings,
                    dynamicData,
                )
                growthTubers[d, 1, y] = getPhotosynthesisPLANTDay( #[g / g*d]
                    d,
                    superIndTubers[d-1, 4, y],
                    superIndSeeds[d-1, 4, y],
                    ((1 - settings["rootShootRatio"]) * (superIndTubers[d-1, 1, y]-superIndTubers[d-1, 5, y]-superIndTubers[d-1, 6, y])),
                    ((1 - settings["rootShootRatio"]) * (superIndSeeds[d-1, 1, y]-superIndSeeds[d-1, 5, y]-superIndSeeds[d-1, 6, y])),
                    LevelOfGrid,
                    settings,
                    dynamicData,
                )

                growthSeeds[d, 3, y] = getDailyGrowth(
                    seeds[d, 3, y], #SeedGerminating Biomass
                    superIndSeeds[d-1, 1, y], #Yesterdays Biomass
                    (superIndSeeds[d-1, 5, y]), #Yesterdays allocatedBiomass from for seeds
                    growthSeeds[d, 1, y], # Todays PS
                    growthSeeds[d, 2, y], # Todays Resp
                    settings,
                )
                growthTubers[d, 3, y] = getDailyGrowth(
                    tubers[d, 3, y], #TubersGerminating Biomass
                    superIndTubers[d-1, 1, y], #Yesterdays Biomass
                    (superIndTubers[d-1, 5, y]), #Yesterdays allocatedBiomass from for seeds
                    growthTubers[d, 1, y], # Todays PS
                    growthTubers[d, 2, y], # Todays Resp
                    settings,
                )

                ##Consequence of negative growth Seeds
                #if growthSeeds[d, 3, y] < 0 &&
                  # superIndSeeds[d-1, 3, y] > 0  #Check if indWeight>0
                  # Mort = (-growthSeeds[d, 3, y] / superIndSeeds[d-1, 3, y])
                 #  if Mort >1
                #       Mort = 1
                 #  end
                #    superIndSeeds[d, 2, y] = #Loss in number of Plants
                #        (1-Mort)*superIndSeeds[d,2,y]
                #    superIndSeeds[d, 1, y] = superIndSeeds[d-1, 1, y] #Biomass stays the same. Makes sense?
                #else
                #    superIndSeeds[d, 1, y] = superIndSeeds[d-1, 1, y] + growthSeeds[d, 3, y] #Biomass
                #end
                superIndSeeds[d, 1, y] = superIndSeeds[d-1, 1, y] + growthSeeds[d, 3, y]
                ##Consequence of negative growth TUBERS
                #if growthTubers[d, 3, y] < 0 &&
                #   superIndTubers[d-1, 3, y] > 0  #Check if indWeight>0
                #   Mort = (-growthTubers[d, 3, y] / superIndTubers[d-1, 3, y])
                #   if Mort >1
                #       Mort = 1
                #   end
                #   superIndTubers[d, 2, y] = (1-Mort)*superIndTubers[d, 2, y]
                #   superIndTubers[d, 1, y] = superIndTubers[d-1, 1, y] #Biomass stays the same. Makes sense?
                #else
                #    superIndTubers[d, 1, y] = superIndTubers[d-1, 1, y] + growthTubers[d, 3, y] #Biomass
                #end
                superIndTubers[d, 1, y] = superIndTubers[d-1, 1, y] + growthTubers[d, 3, y]

                #MORTALITY
                Mort = dieWaves(d,LevelOfGrid,settings,dynamicData) + settings["BackgroundMort"] #+Herbivory
                #if superIndSeeds[d-1, 4, y] < settings["heightMax"] #Check if plant is not yet adult
                #    superIndSeeds[d, 2, y] = (1-Mort)*superIndSeeds[d, 2, y] #Lost in N
                #else #Plant is adult -> lost of weight
                #    superIndSeeds[d, 1, y] = (1-Mort)*superIndSeeds[d, 1, y]
                #end
                superIndSeeds[d, 2, y] = (1-Mort)*superIndSeeds[d, 2, y] #Lost in N
                superIndSeeds[d, 1, y] = (1-Mort)*superIndSeeds[d, 1, y]
                #if superIndTubers[d-1, 4, y] < settings["heightMax"] #Check if plant is not yet adult
                #    superIndTubers[d, 2, y] = (1-Mort)*superIndTubers[d, 2, y]  #Lost in N
                #else #Plant is adult -> lost of weight
                #    superIndTubers[d, 1, y] = (1-Mort)*superIndTubers[d, 1, y]
                #end
                superIndTubers[d, 2, y] = (1-Mort)*superIndTubers[d, 2, y]  #Lost in N
                superIndTubers[d, 1, y] = (1-Mort)*superIndTubers[d, 1, y]
                # CALCULATION OF Ind Weight
                superIndSeeds[d, 3, y] = getIndividualWeight(superIndSeeds[d, 1, y], superIndSeeds[d, 2, y]) #individualWeight = Biomass / Number
                superIndTubers[d, 3, y] = getIndividualWeight(superIndTubers[d, 1, y], superIndTubers[d, 2, y]) #individualWeight = Biomass / Number

                if superIndSeeds[d, 3, y] < 0
                    superIndSeeds[d, 3, y] =0
                end
                if superIndTubers[d, 3, y] < 0
                    superIndTubers[d, 3, y] =0
                end

                #Thinning, optional #TODO SINNVOLL SO?
                if settings["thinning"] == true
                    thin = dieThinning(superIndSeeds[d, 2, y], superIndSeeds[d, 3, y], settings) #Adapts number of individuals [/m^2]& individual weight
                    if (thin[1] < superIndSeeds[d, 2, y]) #&& (Thinning[2] > 0)
                        superIndSeeds[d, 2, y] = thin[1] #N
                        superIndSeeds[d, 3, y] = thin[2] #indWeight
                    end
                    thin = dieThinning(superIndTubers[d, 2, y], superIndTubers[d, 3, y], settings) #Adapts number of individuals [/m^2]& individual weight
                    if (thin[1] < superIndTubers[d, 2, y]) #&& (Thinning[2] > 0)
                        superIndTubers[d, 2, y] = thin[1] #N
                        superIndTubers[d, 3, y] = thin[2] #indWeight
                    end
                end

                #Die-off if N<1
                if superIndSeeds[d, 2, y] ==0 #< 1
                    superIndSeeds[d, 2, y] = 0
                    superIndSeeds[d, 1, y] = 0
                    superIndSeeds[d, 3, y] = 0
                    superIndSeeds[d, 4, y] = 0
                end #no half individuals
                if superIndTubers[d, 2, y] ==0 #< 1
                    superIndTubers[d, 2, y] = 0
                    superIndTubers[d, 1, y] = 0
                    superIndTubers[d, 3, y] = 0
                    superIndTubers[d, 4, y] = 0
                end #no half individuals

                #HEIGHT
                superIndSeeds[d, 4, y] = growHeight(superIndSeeds[d, 3, y],settings)
                if superIndSeeds[d, 4, y] >= settings["heightMax"]
                    superIndSeeds[d, 4, y] = settings["heightMax"]
                end
                if superIndSeeds[d, 4, y] >= WaterDepth
                    superIndSeeds[d, 4, y] = WaterDepth
                end

                superIndTubers[d, 4, y] = growHeight(superIndTubers[d, 3, y],settings)
                if superIndTubers[d, 4, y] >= settings["heightMax"]
                    superIndTubers[d, 4, y] = settings["heightMax"]
                end
                if superIndTubers[d, 4, y] >= WaterDepth
                    superIndTubers[d, 4, y] = WaterDepth
                end

                #ALLOCATION OF BIOMASS FOR SEED PRODUCTION from both INDIVIDUUMS
                if d > (germinationDay + settings["seedsStartAge"]) &&
                   d < (germinationDay + settings["seedsEndAge"]) #Age=Age of plant
                    superIndSeeds[d, 5, y] =
                        settings["seedFraction"] * ((d-germinationDay - settings["seedsStartAge"]) /
                        (settings["seedsEndAge"] - settings["seedsStartAge"])) * superIndSeeds[d, 1, y] # Copied from code!

                    superIndTubers[d, 5, y] =
                        settings["seedFraction"] * ((d-germinationDay - settings["seedsStartAge"]) /
                        (settings["seedsEndAge"] - settings["seedsStartAge"])) * superIndTubers[d, 1, y] # Copied from code!

                end
                if d >= (germinationDay + settings["seedsEndAge"]) &&
                   d <= settings["reproDay"]
                    superIndSeeds[d, 5, y] = settings["seedFraction"] * superIndSeeds[d, 1, y]  #allocatedBiomass - Fraction stays the same
                    superIndTubers[d, 5, y] = settings["seedFraction"] * superIndTubers[d, 1, y]  #allocatedBiomass - Fraction stays the same

                end

                #ALLOCATION OF BIOMASS FOR TUBERS PRODUCTION from both Individuums
                if d > (settings["tuberGerminationDay"] + settings["tuberStartAge"]) &&
                   d < (settings["tuberGerminationDay"] + settings["tuberEndAge"]) #Age=Age of plant
                    superIndSeeds[d, 6, y] =
                        settings["tuberFraction"] * ((d-settings["tuberGerminationDay"] - settings["tuberStartAge"]) /
                        (settings["tuberEndAge"] - settings["tuberStartAge"])) * superIndSeeds[d, 1, y] # Copied from code!
                    superIndTubers[d, 6, y] =
                        settings["tuberFraction"] * ((d-settings["tuberGerminationDay"] - settings["tuberStartAge"]) /
                        (settings["tuberEndAge"] - settings["tuberStartAge"])) * superIndTubers[d, 1, y] # Copied from code!
                end
                if d >= (settings["tuberGerminationDay"] + settings["seedsEndAge"]) &&
                   d <= settings["reproDay"]
                    superIndSeeds[d, 6, y] = settings["tuberFraction"] * superIndSeeds[d, 1, y]  #allocatedBiomass - Fraction stays the same
                    superIndTubers[d, 6, y] = settings["tuberFraction"] * superIndTubers[d, 1, y]  #allocatedBiomass - Fraction stays the same
                end

                # REPRODUCTION DAY
                if d == settings["reproDay"]
                    seeds[d, 1, y] =
                        seeds[d-1, 1, y] + superIndSeeds[d, 5, y] + superIndTubers[d, 5, y] -
                        seeds[d-1, 1, y] * settings["seedMortality"]
                    superIndSeeds[d, 5, y] = 0 #Allocated Biomass is lost
                    superIndTubers[d, 5, y] = 0 #Allocated Biomass is lost

                    tubers[d, 1, y] =
                        tubers[d-1, 1, y] + superIndTubers[d, 6, y] + superIndSeeds[d, 6, y]-
                        tubers[d-1, 1, y] * settings["tuberMortality"]
                    superIndSeeds[d, 6, y] = 0 #Allocated Biomass is lost
                    superIndTubers[d, 6, y] = 0 #Allocated Biomass is lost

                    superIndSeeds[d, 1, y] = superIndSeeds[d,1,y] - superIndSeeds[d, 5, y]  - superIndSeeds[d, 6, y] #Loss of total Biomass as seeds are distributed
                    superIndTubers[d, 1, y] = superIndTubers[d,1,y] - superIndTubers[d, 5, y]- superIndTubers[d, 6, y] #Loss of total Biomass as seeds are distributed

                end

                #Seasonal die-off
                if d == (settings["tuberGerminationDay"] + settings["maxAge"])
                    superIndTubers[d, 1, y] = 0
                    superIndTubers[d, 2, y] = 0
                    superIndTubers[d, 3, y] = 0
                    superIndTubers[d, 4, y] = 0
                    superIndTubers[d, 5, y] = 0
                end

                if d == (germinationDay + settings["maxAge"])
                    superIndSeeds[d, 1, y] = 0
                    superIndSeeds[d, 2, y] = 0
                    superIndSeeds[d, 3, y] = 0
                    superIndSeeds[d, 4, y] = 0
                    superIndSeeds[d, 5, y] = 0
                end

            end #End growth both superInds

        end #Loop over days
    end #Loop over years
    superInd = superIndSeeds + superIndTubers
    return (superInd,superIndSeeds, superIndTubers, seeds, tubers, growthSeeds, growthTubers) #SINNVOLL? ,seeds, growth, lightAttenuation
    #return(superInd)
end

"""
    simulate1Depth(settings)

Simulates 1 depth and returns results in a sturctured manner
"""
function simulate1Depth(depth, settings::Dict{String,Any}, dynamicData::Dict{Int16, DayData})
    #println(depth)
    Res = simulate(depth, settings, dynamicData)
    ResA = Res[1][:, :, 1] #superInd[day,parameter,year]
    ResB = Res[2][:, :, 1] #superIndSeeds[day,parameter,year]
    ResC = Res[3][:, :, 1] #superIndTubers[day,parameter,year]
    ResD = Res[4][:, :, 1] #seeds[day,parameter,year]
    ResE = Res[5][:, :, 1] #tubers[day,parameter,year]
    ResF = Res[6][:, :, 1] #growthSeeds[day,parameter,year]
    ResG = Res[7][:, :, 1] #growthTubers[day,parameter,year]
    for y = 2:settings["years"]
        ResA = vcat(ResA, Res[1][:, :, y]) #superInd
        ResB = vcat(ResB, Res[2][:, :, y]) #superIndSeeds
        ResC = vcat(ResC, Res[3][:, :, y]) #superIndTubers
        ResD = vcat(ResD, Res[4][:, :, y]) #seeds
        ResE = vcat(ResE, Res[5][:, :, y]) #tubers
        ResF = vcat(ResF, Res[6][:, :, y]) #growthSeeds
        ResG = vcat(ResG, Res[7][:, :, y]) #growthTubers
    end
    return ResA, ResB, ResC, ResD, ResE, ResF, ResG #superInd,superIndSeeds, superIndTubers, seeds, tubers, growthSeeds, growthTubers
end #[day*year,parameter], [day*year,parameter]





"""
    simulateMultipleDepth(settings)

Simulates multiple depth and returns Res[depth][dataset][day,parameter]
"""
function simulateMultipleDepth(depths,settings::Dict{String,Any}, dynamicData::Dict{Int16, DayData})
    Res = []
    for d in depths
        push!(Res, simulate1Depth(d,settings,dynamicData))
    end
    return Res
end

"""
    simulateMultipleDepth_parallel(settings)

Simulates multiple depth and returns Res[depth][dataset][day,parameter]
"""
function simulateMultipleDepth_parallel(depths,settings::Dict{String,Any}, dynamicData::Dict{Int16, DayData})
    Res = []
    de = zeros(length(depths))
    Threads.@threads for d in 1:length(depths)
        de[d]=depths[d]
        push!(Res, simulate1Depth(de[d],settings,dynamicData))
        #println(de[d])
    end
    return Res
end



"""
    simulateEnvironment(settings)

Function to calculate the environment with identical input variables as function simulate

Code following the description of Charisma by van Nes et al. (2003)

Arguments used from settings: yearlength, ...

Returns: temp, irradiance, waterlevel, lightAttenuation []
"""

function simulateEnvironment(settings::Dict{String, Any}, dynamicData::Dict{Int16, DayData})
    temp = Float64[]
    irradiance = Float64[]
    waterlevel = Float64[]
    lightAttenuation = Float64[]
    #for y = 1:settings["years"]
        for d = 1:settings["yearlength"]
            dynamicData[d] = DayData()
            push!(temp, getTemperature(d,settings,dynamicData))
            push!(irradiance, getSurfaceIrradianceDay(d,settings,dynamicData))
            push!(waterlevel, getWaterlevel(d,settings,dynamicData))
            push!(lightAttenuation, getLightAttenuation(d,settings, dynamicData))
        end
    #end
    return (temp, irradiance, waterlevel, lightAttenuation)
end
