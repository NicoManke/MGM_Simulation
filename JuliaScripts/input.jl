#Macrophyte Growth Model (MGM)
#
#Anne Lewerentz <anne.lewerentz@uni-wuerzburg.de>
#(c) 2021-2022, licensed under the terms of the MIT license
#
#Contains functions that are necessary to read in input files

"""
    getsettings()

Combines all configuration options to produce a single settings dict.
Order of precedence: config file - default values

Code source: Leidinger et al. (2021)
"""

function getsettings(configfileLake::String = "",configfileSpecies::String = "",configfileGeneral::String = "",)
    defaultsGlobal = defaultSettingsGlobal()
    defaultsLake = defaultSettingsLake()
    defaultsSpecies = defaultSettingsSpecies()

    if !isempty(configfileLake) && isfile(configfileLake)
        configsLake = parseconfigLake(configfileLake)
    else
        configsLake = Dict{String,Any}()
    end
    if !isempty(configfileSpecies) && isfile(configfileSpecies)
        configsSpecies = parseconfigSpecies(configfileSpecies)
    else
        configsSpecies = Dict{String,Any}()
    end

    settings = merge(defaultsGlobal, defaultsLake, configsLake, defaultsSpecies, configsSpecies)
    return settings
end



"""
    testSettings
Function to test for logic of settings to break the look for faster optimization
Further improvement: Print output if not fulfilled
"""
function testSettings(settings::Dict{String, Any})
    # Check setting for logic input
    x=0
    # Environment
    if settings["minI"] > settings["maxI"]  # what would be wrong
        x=+1
        println("ERROR: minI > maxI")
    end
    if settings["minTemp"] > settings["maxTemp"]
        x=+1
        println("ERROR: minTemp > maxTemp")
    end
    if settings["minKd"] > settings["maxKd"]
        x=+1
        println("ERROR: minKd > maxKd")
    end
    if settings["minW"] > settings["maxW"]
        x=+1
        println("ERROR: minW > maxW")
    end

    # Species
    if settings["seedsEndAge"] <= settings["seedsStartAge"]
        x=+1
        println("ERROR: seedsEndAge < seedsStartAge")
    end
    if settings["tuberEndAge"] <= settings["tuberStartAge"]
        x=+1
        println("ERROR: tuberEndAge < tuberStartAge")
    end
    if settings["maxAge"] <= settings["seedsEndAge"]  #add for tubers?
        x=+1
        println("ERROR: maxAge < seedsStartAge")
    end
    if settings["reproDay"] <= settings["germinationDay"] + settings["seedsEndAge"]
        x=+1
        println("ERROR: reproDay < germinationDay + seedsEndAge")
    end
    if settings["reproDay"] >= settings["germinationDay"] + settings["maxAge"]
        x=+1
        println("ERROR: reproDay > germinationDay + maxAge")
    end
    return x
end


#testSettings(settings)


"""
    basicparser(filename)

Do elementary parsing on a config or map file.

Reads in the file, strips whole-line and inline comments
and separates lines by whitespace.
Returns a 2d array representing the tokens in each line.

Code source: Leidinger et al. (2021)

"""
function basicparser(filename::String)
    # Read in the file
    lines = String[]
    open(filename) do file
        lines = readlines(file)
    end
    # Remove comments and tokenize
    lines = map(x -> strip(x), lines)
    filter!(x -> !isempty(x), lines)
    filter!(x -> (x[1] != '#'), lines)
    lines = map(s -> strip(split(s, '#')[1]), lines)
    lines = map(split, lines)
    map(l -> map(s -> convert(String, s), l), lines)
end


"""
    advancedparser(filename)

Do elementary parsing on a config or map file.

Code source: Leidinger et al. (2021)

"""
function advancedparser(filename::String)
    # Read in the file
    lines = String[]
    open(filename) do file
        lines = readlines(file)
    end
    # Remove comments and tokenize
    lines = map(x -> strip(x), lines)
    filter!(x -> !isempty(x), lines)
    filter!(x -> (x[1] != '#'), lines)
    lines = map(s -> strip(split(s, '#')[1]), lines)
    lines = map(split, lines)
    map(l -> map(s -> convert(String, s), l), lines)
end



"""
    parseconfigLake(filename)

Parse a configuration file and return a settings dict.

The config syntax is very simple: each line consists of a parameter
name and a value (unquoted), e.g. `nniches 2`. `#` is the comment character.

Code source: Leidinger et al. (2021)

"""

function parseconfigLake(configfilename::String)
    config = basicparser(configfilename)
    settings = Dict{String, Any}()
    defaults = defaultSettingsLake()
    for c in config
        if length(c) != 2
            #simlog("Bad config file syntax: $c", settings, 'w', "")
        elseif c[1] in keys(defaults)
            value = c[2]
            if !(typeof(defaults[c[1]]) <: AbstractString)
                try
                    value = parse(typeof(defaults[c[1]]), c[2]) # or Meta.parse with the old functionality
                catch
                    #simlog("$(c[1]) not of type $(typeof(defaults[c[1]])).",
                    #       settings, 'w', "")
                end
            end
            settings[c[1]] = value
        else
            # XXX maybe parse anyway
            #simlog(c[1]*" is not a recognized parameter!", settings, 'w', "")
        end
    end
    settings
end



"""
    parseconfigSpecies(filename)

Parse a configuration file and return a settings dict.

The config syntax is very simple: each line consists of a parameter
name and a value (unquoted), e.g. `nniches 2`. `#` is the comment character.

Code source: Leidinger et al. (2021)

"""

function parseconfigSpecies(configfilename::String)
    config = basicparser(configfilename)
    settings = Dict{String, Any}()
    defaults = defaultSettingsSpecies()
    for c in config
        if length(c) != 2
            #simlog("Bad config file syntax: $c", settings, 'w', "")
        elseif c[1] in keys(defaults)
            value = c[2]
            if !(typeof(defaults[c[1]]) <: AbstractString)
                try
                    value = parse(typeof(defaults[c[1]]), c[2]) # or Meta.parse with the old functionality
                catch
                    #simlog("$(c[1]) not of type $(typeof(defaults[c[1]])).",
                    #       settings, 'w', "")
                end
            end
            settings[c[1]] = value
        else
            # XXX maybe parse anyway
            #simlog(c[1]*" is not a recognized parameter!", settings, 'w', "")
        end
    end
    settings
end


"""
    parseconfigGeneral(filename)

Parse a configuration file and return a settings dict.

The config syntax is very simple: each line consists of a parameter
name and a value (unquoted), e.g. `nniches 2`. `#` is the comment character.

Code source: Leidinger et al. (2021)
"""
function parseconfigGeneral(configfilename::String)
    config = advancedparser(configfilename)
    settings = Dict{String, Any}()
    defaults = defaultSettingsGeneral()
    for c in config
        if c[1] in keys(defaults)
            value=String[]
            for i in 2:length(c)
                push!(value,c[i])
            end
            settings[c[1]] = value
        end
    end
    settings
end
