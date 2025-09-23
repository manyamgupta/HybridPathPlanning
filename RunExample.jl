# Urban AirMobility path planning 
# Modeled as hybrid path planning that includes uav traveling using gas or electric
# edges in quiet zones are required to use electric mode only
# decsions involves to use gas or electric or partly gas/electric on non-qiet edges
# travel cost is the fuel cost, electric segments does not have a cost
# hybrid path planning graph is generated with discretized battery charge levels
# script to run the instances
using LinearAlgebra
using PyCall
using PyPlot
using FileIO, JLD, JLD2
using Graphs
using LightGraphs
include("utils.jl")

include("HybPP_LB_v3.jl")
include("HybPP_UB_v3.jl")

pygui(true)


mapFileName = string("boston2.jld")
map = FileIO.load(mapFileName)

obstaclesList = map["obsList"]
sceneFileName = string("boston2_scenarios.jld2")
scene = FileIO.load(sceneFileName)
startPtsVec = scene["startPtsVec"]
goalPtsVec = scene["goalPtsVec"]

k=2
startPos = startPtsVec[k,:]
goalPos = goalPtsVec[k,:]

delLength = 400.
# number of partitions the charge levels are divided into for the lower bound
numchargePartns = 20
n_discs_charge = numchargePartns+1
numSubPartns = 5
discharge_rate = .04 #battery discharge rate per unit distance in electric mode
recharge_rate = .02 #battery charge rate per unit distance in gas mode
min_charge_at_dest = 50. #minimum charge required at destination
qmin = 20
qmax = 100

battery_params = (discharge_rate=discharge_rate, recharge_rate=recharge_rate, min_charge_at_dest=min_charge_at_dest, qmin=qmin, qmax=qmax)
deltaQ = (qmax-qmin)/numchargePartns
chargePartitions = [((i-1)*deltaQ+qmin, i*deltaQ+qmin) for i=1:numchargePartns]


resultsHeaders = ["Inst" "LB" "UB" "UB_nfqz" "PcGap" "pctCostSvng" "TimeOfflineLB" "TimeOnlineLB" "TimeTotalLB" "TimeOfflineUB" "TimeOnlineUB"]
resultsList = zeros(50,11)

plotFlag=true


stats1 = @timed HppBaseGrLB(obstaclesList, delLength, battery_params, chargePartitions)
baseGrParamsLB = stats1.value
compTimeOffl_LB = stats1.time

stats2 = @timed HppLB(baseGrParamsLB, startPos, goalPos, battery_params, numSubPartns)
resultsLB = stats2.value
compTimeOnl_LB = stats2.time

pathCostLB = resultsLB.pathCostLB

# if plotFlag
#     @show pathCostLB    
#     spLB = resultsLB.sp
#     nodeCoordsListLB = resultsLB.grParams.nodeCoordsList    
#     PyPlot.figure()
#     PlotObstacles(obstaclesList)
#     PlotGraphNodesLB(baseGrParamsLB.vertCoordsList)
#     PlotPath(spLB, nodeCoordsListLB, battery_params)
#     PyPlot.axis("equal")

#     PyPlot.figure()
#     PlotChargeProfileLB(spLB, nodeCoordsListLB, battery_params)
#     PyPlot.xlabel("distance traveled")
#     PyPlot.ylabel("State of Charge")
# end
    # PrintPathNodes(spLB, nodeCoordsListLB, battery_params)

flyOverQuietZones = true
stats3 = @timed HppBaseGrUB(obstaclesList, delLength, battery_params, n_discs_charge, [])
baseGrParamsUB = stats3.value
compTimeOffl_UB = stats3.time

stats4 = @timed HppUB(baseGrParamsUB, startPos, goalPos, battery_params, n_discs_charge)
resultsUB = stats4.value
compTimeOnl_UB = stats4.time
pathCostUB = resultsUB.pathCostUB

println("Lower bound of the optimum: ", pathCostLB)
println("Solution cost (upper bound): ", pathCostUB)

if plotFlag    
    spUB = resultsUB.sp
    nodeCoordsListUB = resultsUB.grParams.nodeCoordsList
    PyPlot.figure()
    PlotObstacles(obstaclesList)
    PlotGraphNodesUB(baseGrParamsUB.vertCoordsList)
    PlotPath(spUB, nodeCoordsListUB, battery_params)
    PyPlot.axis("equal")

    PyPlot.figure()
    PlotChargeProfileUB(spUB, nodeCoordsListUB, battery_params)
    PyPlot.xlabel("distance traveled")
    PyPlot.ylabel("State of Charge")
end
   
