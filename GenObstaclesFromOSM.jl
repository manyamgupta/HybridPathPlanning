using LightXML
using Geodesy
using PyPlot
using FileIO

pygui(true)
xdoc = parse_file("../BenchmarkMaps/newyork0.osm")

xroot = root(xdoc)

bldgWayIds = String[]

for e in xroot["way"]   
    for cd in child_elements(e)
        if has_attribute(cd,"k")
            if attributes_dict(cd)["k"] == "building" && attributes_dict(cd)["v"] == "yes"
                # println(attributes_dict(cd)["v"])
                # println(attributes_dict(e)["id"])
                push!(bldgWayIds, attributes_dict(e)["id"])
            end
        end
    end    
end

bldgsNodeIdsVec = Vector{Vector{String}}()
allBldgNodeIds = Vector{String}(undef,0)
for e in xroot["way"]   
    if attributes_dict(e)["id"] in bldgWayIds
        @show e
        bldgNodeIds = Vector{String}(undef,0)
        for cd in child_elements(e)
            if has_attribute(cd,"ref")
                # @show attributes_dict(cd)["ref"]
                push!(bldgNodeIds, attributes_dict(cd)["ref"])
                push!(allBldgNodeIds, attributes_dict(cd)["ref"])
            end
        end
        push!(bldgsNodeIdsVec, bldgNodeIds)
    end
end

obsCoordsVec = Matrix{Float64}[]
nodeToCoordsDict = Dict("String"=>[0.;0.])

for k = 1:length(xroot["node"])
    global nodeToCoordsDict    
    n = xroot["node"][k]
    nodeId = attributes_dict(n)["id"] # in allBldgNodeIds
        # @show n
    lat = parse(Float64, attributes_dict(n)["lat"])
    lon = parse(Float64, attributes_dict(n)["lon"])
    node_lla = LLA(lat, lon,0)
    # @show node_lla
    node_ecef = ECEF(node_lla, wgs84)
    # @show node_ecef
    nodeCoords = Vector(node_ecef)
    # @show nodeCoords
    # obsCoords = [obsCoords; nodeCoords']
    nodeToCoordsDict = merge(nodeToCoordsDict, Dict(nodeId=>nodeCoords))

end

for bldgNodeIds in bldgsNodeIdsVec
    obsCoords = Matrix{Float64}(undef,0,2)
    for nodeId in bldgNodeIds
        obsCoords = [obsCoords; nodeToCoordsDict[nodeId][1:2]']
    end
    push!(obsCoordsVec, obsCoords)
end

xmin = minimum([minimum(obsCoordsVec[i][:,1]) for i=1:length(obsCoordsVec)])
ymin = minimum([minimum(obsCoordsVec[i][:,2]) for i=1:length(obsCoordsVec)])

obsList = Matrix{Float64}[]
for i=1:length(obsCoordsVec)
    obs = [obsCoordsVec[i][:,1].-xmin obsCoordsVec[i][:,2].-ymin]*20
    push!(obsList, obs)
end

noObs = length(obsList)
for j=1:noObs
    obs = obsList[j]
    obs = [obs;obs[1,:]']
    PyPlot.plot(obs[:,1], obs[:,2], color="red", linewidth=1.5)

end

FileIO.save("newyork0.jld","obsList",obsList)