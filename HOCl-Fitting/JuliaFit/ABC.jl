function defineInternalCoordinates(zMatrixCoordinates::Vector{Float64})::Vector{Float64}
    convertToRadians::Float64 = pi/180.0
    # Trivial case - no dihedrals to convert just degrees to radians
    internalCoordinates::Vector{Float64} = zeros(length(zMatrixCoordinates))
    internalCoordinates[1:2] = zMatrixCoordinates[1:2]
    internalCoordinates[3] = zMatrixCoordinates[3]*convertToRadians
    return internalCoordinates
end

function convertAnglesToRadians(structuralParameters::Vector{Float64})::Vector{Float64}
    convertToRadians::Float64 = pi/180.0
    structuralParameters[3] = structuralParameters[3]*convertToRadians
    return structuralParameters
end

function convertAnglesToDegrees(structuralParameters::Vector{Float64})::Vector{Float64}
    convertToRadians::Float64 = pi/180.0
    structuralParameters[3] = structuralParameters[3]/convertToRadians
    return structuralParameters
end

function defineXiCoordinates(internalCoordinates::Vector{Float64}, structuralParameters)
    T::DataType = eltype(structuralParameters)
    convertToRadians::Float64 = pi/180
    numberOfModes::Int64 = length(internalCoordinates)
    numberOfAtoms::Int64 = (numberOfModes + 6)/3
    numberOfStretches::Int64 = numberOfAtoms - 1
    xi::Vector{T} = zeros(numberOfModes)
    equilibriumParameters = structuralParameters[1:numberOfModes]
    morseParameters = structuralParameters[numberOfModes+1:end]
    
    # Define Morse-like coordinates for stretches
    stretchDisplacement = internalCoordinates[1:numberOfStretches] - equilibriumParameters[1:numberOfStretches]
    xi[1:numberOfStretches] = 1 .- exp.(-morseParameters.*stretchDisplacement)
    # xi[1:numberOfStretches] = stretchDisplacement

    # xi[numberOfStretches+1:end] = internalCoordinates[numberOfStretches+1:end] - equilibriumParameters[numberOfStretches+1:end]
    # # Trigonometric-type function for the bending 
    xi[numberOfStretches+1:end] = cos.(internalCoordinates[numberOfStretches+1:end]) .- cos.(equilibriumParameters[numberOfStretches+1:end])
    # xi[numberOfStretches+1:end] = cos.(pi .- internalCoordinates[numberOfStretches+1:end].*convertToRadians) .- cos.(pi .- equilibriumParameters[numberOfStretches+1:end].*convertToRadians)
    # xi[numberOfStretches+1:end] = sin.(pi .- internalCoordinates[numberOfStretches+1:end].*convertToRadians) .- cos.(pi .- equilibriumParameters[numberOfStretches+1:end].*convertToRadians)
    return xi
end

# A function for pre-scaling for Chebyshev polynomials, use integration ranges defined in TROVE
function defineXiCoordinatesScaled(internalCoordinates::Vector{Float64}, structuralParameters)
    T::DataType = eltype(structuralParameters)
    convertToRadians::Float64 = pi/180
    numberOfModes::Int64 = length(internalCoordinates)
    numberOfAtoms::Int64 = (numberOfModes + 6)/3
    numberOfStretches::Int64 = numberOfAtoms - 1
    xi::Vector{T} = zeros(numberOfModes)
    equilibriumParameters = structuralParameters[1:numberOfModes]
    morseParameters = structuralParameters[numberOfModes+1:end]
    
    # Define Morse-like coordinates for stretches
    stretchDisplacement = internalCoordinates[1:numberOfStretches] - equilibriumParameters[1:numberOfStretches]
    xi[1:numberOfStretches] = 1 .- exp.(-morseParameters.*stretchDisplacement)
    # xi[1:numberOfStretches] = stretchDisplacement

    # xi[numberOfStretches+1:end] = internalCoordinates[numberOfStretches+1:end] - equilibriumParameters[numberOfStretches+1:end]
    # # Trigonometric-type function for the bending 
    xi[numberOfStretches+1:end] = cos.(internalCoordinates[numberOfStretches+1:end]) .- cos.(equilibriumParameters[numberOfStretches+1:end])
    # xi[numberOfStretches+1:end] = cos.(pi .- internalCoordinates[numberOfStretches+1:end].*convertToRadians) .- cos.(pi .- equilibriumParameters[numberOfStretches+1:end].*convertToRadians)
    # xi[numberOfStretches+1:end] = sin.(pi .- internalCoordinates[numberOfStretches+1:end].*convertToRadians) .- cos.(pi .- equilibriumParameters[numberOfStretches+1:end].*convertToRadians)

    # Rescale within -1 and 1 for Chebyshev basis
    upperBounds::Vector{Float64} = zeros(3)
    lowerBounds::Vector{Float64} = zeros(3)

    # -0.34,2.5
    # -0.40,1.8
    upperBounds[1] = 1 - exp(morseParameters[1]*0.34)
    upperBounds[2] = 1 - exp(morseParameters[2]*0.40)
    lowerBounds[1] = 1 - exp(-morseParameters[1]*2.5)
    lowerBounds[2] = 1 - exp(-morseParameters[2]*1.8)
    bendingRangeForRescale::Vector{Float64} = collect(LinRange(40, 180, 200))
    bendingRangeForRescale = cos.(bendingRangeForRescale.*convertToRadians) .- cos(equilibriumParameters[end])
    upperBounds[3] = maximum(bendingRangeForRescale)
    lowerBounds[3] = minimum(bendingRangeForRescale)

    xi = (2*xi .- (upperBounds .+ lowerBounds))./(upperBounds .- lowerBounds)

    return xi
end