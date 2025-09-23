# Author: S Gupta Manyam
using LinearAlgebra
import LazySets.LineSegment
import LazySets.is_intersection_empty
using SparseArrays
# using MetaGraphs

mutable struct Obstacle
    vertexList::Matrix{Float64} #vertices of the polygonal obstacle
    delLen::Float64 # descretization parameter 
    discretePointsList::Matrix{Float64} #each edge is dicretized into a list of points    
    id::Int64
    sideIdsList::Matrix{Int64}

    function Obstacle(vList::Matrix{Float64}, delL::Float64, obsId::Int64)
        # dpList, sIdsList = DiscretizeEdges(vList, delL)
        # new(vList, delL, dpList, sIdsList, obsId)
        dpList, sIdsList = DiscretizeEdges(vList, delL)
        new(vList, delL, dpList, obsId, sIdsList)    
    end

    function DiscretizeEdges(vList::Matrix{Float64}, delL::Float64)

        numEdgs = size(vList, 1)
        vList = [vList; vList[1,:]']
        dps = Matrix{Float64}(undef,0,2)
        sIds = Matrix{Int64}(undef,0,2)
        
        for i=1:size(vList,1)-1
            edgLen = LinearAlgebra.norm(vList[i,:]-vList[i+1,:])
            numDiscPoints = Int(floor(edgLen/delL))+1
            numDiscPoints = numDiscPoints<2 ? 2 : numDiscPoints
            
            lamVec = LinRange(0,1,numDiscPoints)
            lamVec = lamVec[1:end-1]
            dps_edge = (-lamVec.+1).*vList[i,:]' + lamVec.*vList[i+1,:]'
            
            dps = [dps; dps_edge]
            sIds_curEdg = Int64.([(i*ones(numDiscPoints-1)) zeros(numDiscPoints-1) ])
            sIds_curEdg[1,2] = i-1 == 0 ? numEdgs : i-1
            sIds = [sIds; sIds_curEdg]

        end
        return dps, sIds
    end
end

function IntersectionLineSegments(p1,p2,p3,p4)
    # https://math.stackexchange.com/questions/3176543/intersection-point-of-2-lines-defined-by-2-points-each
    # function returns the point of intersection of line segments p1p2 and p3p4
    n = p2-p1
    m = p3-p4
    p = p3-p1

    D = n[1]*m[2]-n[2]*m[1]
    if abs(D) < 1e-6
        return NaN
    end
    Qx = m[2]*p[1]-m[1]*p[2]
    t = Qx/D
    
    intPt = p1+t*(p2-p1)

    return intPt
end

function InInt(lb, ub, t )
    
    lb = mod(lb, 2*pi)
    ub = mod(ub, 2*pi)
    t = mod(t, 2*pi) 
    if (lb >= ub)
        if(t >= lb || t <= ub)
            return true
        else
            return false
        end
    else
        if(t >= lb && t <= ub)
            return true
        else
            return false
        end
    end
end


function DistanceBetweenNodes(n1,n2,coordsList)
    return LinearAlgebra.norm(coordsList[n1,:]-coordsList[n2,:])
end


function CheckChrgLvlsFeas(b1,b2, batteryParams, edgLen)
    # assumption: first lamda*edgLen distance is on gas, rest is electric
    feasCheck = false
    alpha = batteryParams.discharge_rate
    beta = batteryParams.recharge_rate
    q_min = batteryParams.qmin
    q_max = batteryParams.qmax

    # if b1 == q_max && b2 == q_max # all gas mode
    #     feasCheck = true
    #     lamda = 1.
    #     lambdaVec = [lamda]
    if b2 >= min(b1+beta*edgLen, q_max) # not feasible coz final battery level not attainable
        feasCheck = false
        lamda = 1.
        lambdaVec = [lamda]
    elseif b1-alpha*edgLen > q_min && b2 <= b1-alpha*edgLen #all electric mode
        feasCheck = true
        lamda = 0.
        lambdaVec = [lamda]
    else
        lambda1 = (q_max-b1)/(beta*edgLen)
        lambda2 = (q_max-b2)/(alpha*edgLen)

        if lambda1>=0 && lambda1<=1 && lambda2>=0 && lambda2<=1 && lambda1+lambda2<1
            # lamda = 1.0-lambda2
            # switching twice is cheaper in this case
            lamda = lambda1 + (alpha/(alpha+beta))*(1-lambda1-lambda2)
            lambdaVec = [lambda1, beta*(1-lambda1-lambda2)/(alpha+beta), alpha*(1-lambda1-lambda2)/(alpha+beta) ]
        else
            lamda = (b2-b1+edgLen*alpha)/(edgLen*(alpha+beta))
            lambdaVec = [lamda]
        end
        if lamda>=0 && lamda <= 1
            feasCheck = true
        else
            feasCheck = false
        end
    end
    
    return feasCheck, lambdaVec
end

function CheckLineIntersectsSides(p1,p2,sidesList)

    retVal = false
    p1 = convert(Array{Float64,1}, p1)
    p2 = convert(Array{Float64,1}, p2)
    line_nodes = LineSegment(p1,p2)
    for k=1:length(sidesList)
        side = sidesList[k]
        side_p = convert(Array{Float64,1}, side[1])
        side_q = convert(Array{Float64,1}, side[2])

        line_side = LineSegment(side_p, side_q)
        int_res = is_intersection_empty(line_side, line_nodes, true)
        if int_res[1] == false 
            if round.(int_res[2],digits=4) != round.(p1, digits=4) && round.(int_res[2],digits=4) != round.(p2,digits=4)
                # println("int_res2: $(int_res[2])")
                # println("p1: $p1")
                # println("p2: $p2")
                return true
            end
        end
    end
    return retVal
end

function CheckLineIntersectsSides2(p1,p2,sidesList)
    retVal = false
    p1 = (x=p1[1], y=p1[2])
    p2 = (x=p2[1], y=p2[2])
    
    for k=1:length(sidesList)
        side = sidesList[k]
        p3 = (x=side[1][1], y=side[1][2])
        p4 = (x=side[2][1], y=side[2][2])
        
        retVal = CheckIntersect(p1,p2,p3,p4)
        if retVal
            return true
        end
        # if intrsctFlag == true
        #     intrsPt = IntersectionLineSegs(p1,p2,p3,p4) 
        #     if round.(intrsPt,digits=4) != round.(p1, digits=4) && round.(intrsPt,digits=4) != round.(p2,digits=4)
        #         return true
        #     end
        # end
    end
    return retVal
end

function CheckLineIntersectsSidesStrtGl(p1,p2,sidesList)
    retVal = false
    p1 = (x=p1[1], y=p1[2])
    p2 = (x=p2[1], y=p2[2])
    
    for k=1:length(sidesList)
        side = sidesList[k]
        p3 = (x=side[1][1], y=side[1][2])
        p4 = (x=side[2][1], y=side[2][2])
        
        retVal = CheckIntersectStGl(p1,p2,p3,p4)
        if retVal
            return true
        end
    end
    return retVal
end

function CCW(A,B,C)
    return (C.y-A.y) * (B.x-A.x) > (B.y-A.y) * (C.x-A.x)
end
function CheckIntersect(A,B,C,D)
    if CheckLinear(C,D, A) || CheckLinear(C,D,B)
        return true
    end
    return CCW(A,C,D) != CCW(B,C,D) && CCW(A,B,C) != CCW(A,B,D)
end

function CheckIntersectStGl(A,B,C,D)
    if CheckLinear(C,D, A) || CheckLinear(C,D,B)
        return false
    end
    return CCW(A,C,D) != CCW(B,C,D) && CCW(A,B,C) != CCW(A,B,D)
end

function CheckLinear(A,B,C)
    # returns true if C lies on linesegment AB
    if abs((C.y-A.y) * (B.x-A.x) - (B.y-A.y) * (C.x-A.x)) < 1e-6
        return true
    end
    return false
end

function PlotGraph(coordsList, edgList)
    PyPlot.scatter(coordsList[:,1], coordsList[:,2], marker="x")
    for i = 1:size(coordsList,1)
        PyPlot.text(coordsList[i,1], coordsList[i,2],"$i" )
    end

    for edg in edgList        
        u = edg.src
        v = edg.dst
        PyPlot.plot( [coordsList[u,1],coordsList[v,1] ], [coordsList[u,2],coordsList[v,2] ], linewidth=1,color="blue" )        
    end
end

function EnumEdges(coordsList, edgList, edgCostMatrix)
    for edg in edgList
        u = edg.src
        v = edg.dst
        println("$(coordsList[u,:]) -> $(coordsList[v,:]): cost: $(edgCostMatrix[u,v])")
    end

end

function PlotObstacles(obsList, attr=(color="red", linewidth=1.5))

    noObs = length(obsList)
    for j=1:noObs
        # @show j
        obs = obsList[j]
        obs = [obs;obs[1,:]']
        PyPlot.plot(obs[:,1], obs[:,2], color=attr.color, linewidth=attr.linewidth)
        # readline()
    end
   
end



function PlotPath(sp, nodeCoordsList, batteryParams)

    for i=1:length(sp)
        
        u = sp[i].src
        v = sp[i].dst
        node1 = nodeCoordsList[u,:]
        node2 = nodeCoordsList[v,:]

        if length(node1) == 4
            node1 = node1[[1,2,4]]
            node2 = node2[[1,2,3]]
            # feas, switchPoint = CheckBatteryFeasibilityLB(node1, node2, batteryParams)
        end
        feas, switchPoints = CheckBatteryFeasibilityUB_chrgprf(node1, node2, batteryParams)
        prevSwitchPt = node1[1:2]
        for k=1:size(switchPoints,1)
            switchPoint = switchPoints[k,:]

            xs_gas = [prevSwitchPt[1], switchPoint[1] ]        
            ys_gas = [prevSwitchPt[2], switchPoint[2] ]

            pathcolor = isodd(k) ? "magenta" : "green"
            PyPlot.plot(xs_gas, ys_gas, color=pathcolor, linewidth=2)
            prevSwitchPt = switchPoint
            
        end
        xs_elec = [switchPoints[end,1], node2[1] ]
        ys_elec = [switchPoints[end,2], node2[2] ]
        PyPlot.plot(xs_elec, ys_elec, color="green", linewidth=2)
        
    end
end

function PrintPathNodes(sp, nodeCoordsList, batteryParams)
    for i=1:length(sp)
        
        u = sp[i].src
        v = sp[i].dst
        node_src = nodeCoordsList[u,:]
        node_dst = nodeCoordsList[v,:]
        @show node_src
        @show node_dst
        if length(node_src) == 3
            feas, switchPoint = CheckBatteryFeasibilityUB(node_src, node_dst, batteryParams)
        elseif length(node_src) == 4
            feas, switchPoint = CheckBatteryFeasibilityLB(node_src, node_dst, batteryParams)
        end
        @show switchPoint
        edgCost = norm(switchPoint-node_src[1:2])
        edgLength = norm(node_dst[1:2]-node_src[1:2])
        @show edgCost
        @show edgLength

        
    end
end
function ComputePathCost(sp, costMat::SparseMatrixCSC)

    pathCost = 0.
    for edg in sp
        edgCost = costMat[edg.src, edg.dst]
        edgCost = edgCost == 0.00001 ? 0 : edgCost
        pathCost += edgCost
    end

    return pathCost
end

# function ComputePathCost(sp, graph::MetaDiGraph)

#     pathCost = 0.
#     for edg in sp
#         edgCst = get_prop(graph, edg.src, edg.dst, :weight)
#         pathCost += edgCst
#     end        

#     return pathCost
# end

function PlotCircle(cntr, rho, fmt)
    alvec = LinRange(0,2*pi, 1000)

    xv = cntr[1].+rho.*cos.(alvec)
    yv = cntr[2].+rho.*sin.(alvec)
    PyPlot.plot(xv, yv, color=fmt.color, linewidth=fmt.linewidth)
end

function DiscPosVerts(obsList, delL, start=(NaN, NaN), goal=(NaN, NaN), noFlyInds = [])

    # obsDiscList = Array{Obstacle, 1}(undef, 0)
    sidesList = Array{Tuple{Vector{Float64}, Vector{Float64}}}(undef, 0)
    
    coordsList = Matrix{Float64}(undef,0,2)
    coordsList = [start[1] start[2]; goal[1] goal[2]] #list of all the vertices coordinates the index is the nodeId    
    obsIdList = [-1;0] #list mapping the nodeId to obstacle id, the index is the nodeId, obsId is -1 for start and 0 for the goal
    sideIdsList = [-1 -1;0 0] #list mapping the nodeId to side id, the index is the nodeId, sideId is -1 for start and 0 for the goal, sideId is the index of the side in the sidesList
    lastSideId = 0    
    delLInf = 1e8
    for k = eachindex(obsList)
        obsPolyVerts = obsList[k]
        obsPolyVerts_app = [obsPolyVerts; obsPolyVerts[1,:]']
        # for j = 1:size(obsPolyVerts,1)
        #     push!(sidesList, (obsPolyVerts_app[j,:], obsPolyVerts_app[j+1,:]) )
        # end
        slObs = [(obsPolyVerts_app[j,:], obsPolyVerts_app[j+1,:]) for j in 1:size(obsPolyVerts,1)]
        append!(sidesList, slObs )
        if k in noFlyInds
            obsDisc = Obstacle(obsPolyVerts,delLInf,k) # with very high delL, the discretization only gives vertices of each edge
        else
            obsDisc = Obstacle(obsPolyVerts,delL,k)
        end
        # push!(obsDiscList, obsDisc)
        
        coordsList = [coordsList; obsDisc.discretePointsList ]                                 
        obsIdList = [obsIdList; Int64.(k*ones(size(obsDisc.discretePointsList,1)))]

        lastSideIdList = lastSideId*ones(size(obsDisc.sideIdsList))
        for k in 1:length(lastSideIdList)
            if obsDisc.sideIdsList[k]==0
                lastSideIdList[k]=0
            end
        end
        sideIdsListCur = Int64.(lastSideIdList+ obsDisc.sideIdsList)
        sideIdsList = [sideIdsList; sideIdsListCur]
        lastSideId = sideIdsList[end,1]
    end
    goalObsId = size(obsList,1)+1 #this is the id used for replicas of goal node for each battery level
    obsIdList = [obsIdList; Int64.(goalObsId)]
    
    coordsList = [coordsList; [goal[1] goal[2]] ]
    sideIdsList = [sideIdsList; 0 0]
    return coordsList, obsIdList, sidesList, sideIdsList
end

function CheckPtLiesInObsList(pt, obsList)
    # pt is a Luxor Point
    # obsList is list of obstacles, vector of vector of points
    flag = false
    for k = 1:length(obsList)
        flag = isinside(pt, obsList[k])
        if flag
            return true
        end
    end
    return flag
end