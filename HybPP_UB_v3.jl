# This version can address the no-fly zones
# The function HppBaseGrUB takes noFlyZoneIndices as input, this indicates which among the obstacles are no-flyzones

function HppBaseGrUB(obsList, delL, batteryParams, numDiscs, noFlyInds=Vector{Int64}(undef,0))
    
    vertCoordsList, obsIdList, sidesList, sideIdsList = DiscPosVerts(obsList, delL, (NaN, NaN), (NaN, NaN), noFlyInds)
    
    grPosSmpl, smplGrEdgType = PosSmplgGraphUB(vertCoordsList, obsIdList, sidesList, sideIdsList, noFlyInds)
    
    baseGrParams = SmplgGraphUB(grPosSmpl, smplGrEdgType, vertCoordsList, obsIdList, batteryParams, numDiscs)

    baseGrParams = merge(baseGrParams, (sidesList=sidesList, sideIdsList=sideIdsList, grPosSmpl=grPosSmpl, smplGrEdgType=smplGrEdgType, obsIdList=obsIdList))
    
    return baseGrParams
    
end

function HppUB(baseGrParams, start, goal, batteryParams, numDiscs)

    nodeCoordsList = baseGrParams.nodeCoordsList
    edgWtsMat = baseGrParams.edgWtsMat
    gr_uam = baseGrParams.gr_uam
    chrgSamplesVec = baseGrParams.chrgSamplesVec
    sidesList = baseGrParams.sidesList
    nodesOfVertList = baseGrParams.nodesOfVertList
    nodeIdtoVertIdMap = baseGrParams.nodeIdtoVertIdMap

    sp = LightGraphs.SimpleGraphs.SimpleEdge{Int64}[]
    pathCostUB=0.

    grParams = Tuple{}
    baseGrParams = AddstartGoalToGrPosSmplUB(baseGrParams, (start, goal))
    grParams = AddStartGoalToGraphUB(gr_uam, edgWtsMat, nodeCoordsList, baseGrParams.vertCoordsList, (start, goal), chrgSamplesVec, batteryParams, sidesList)
        
    sp = LightGraphs.a_star(grParams.gr_uam, 1, 2, sparse(grParams.edgWtsMat))
        
    pathCostUB = ComputePathCost(sp, grParams.edgWtsMat)
    # @show pathCostUB

    return (pathCostUB=pathCostUB, sp=sp, chrgSamplesVec=chrgSamplesVec, grParams=grParams, nodesOfVertList=nodesOfVertList, nodeIdtoVertIdMap=nodeIdtoVertIdMap, baseGrParams=baseGrParams)
end

function RemoveEdges(grParams, nodesToBePatnd)
    gr_uam = grParams.gr_uam
    edgWtsMat = grParams.edgWtsMat
    for nd in nodesToBePatnd
        in_nbrs = copy(LightGraphs.inneighbors(gr_uam, nd))
        out_nbrs = copy(LightGraphs.outneighbors(gr_uam, nd))
        [LightGraphs.rem_edge!(gr_uam, nbr_nd, nd) for nbr_nd in in_nbrs]
        [edgWtsMat[nbr_nd, nd]=0. for nbr_nd in in_nbrs]

        [LightGraphs.rem_edge!(gr_uam, nd, nbr_nd) for nbr_nd in out_nbrs]
        [edgWtsMat[nd, nbr_nd]=0. for nbr_nd in out_nbrs]
    end   
    grParams = merge(grParams, (gr_uam=gr_uam, edgWtsMat=edgWtsMat))
    return grParams
end

function AddNewNodesToList(grParams, nodesToBePatnd, numSubPartns, nodeIdtoVertIdMap )
    lastNodeId = length(nodeIdtoVertIdMap)
    newNodesRange = lastNodeId+1:lastNodeId+numSubPartns*length(nodesToBePatnd)
    nodeCoordsListNew = copy(grParams.nodeCoordsList)
    for node_p in nodesToBePatnd
        curNodeCoords = nodeCoordsListNew[node_p,:]
        delQ = (curNodeCoords[4]-curNodeCoords[3])/numSubPartns

        newNodeCoords = [ ones(numSubPartns)*curNodeCoords[1] ones(numSubPartns)*curNodeCoords[2] curNodeCoords[3].+collect(0:numSubPartns-1)*delQ curNodeCoords[3].+collect(1:numSubPartns)*delQ]

        nodeCoordsListNew = [nodeCoordsListNew; newNodeCoords]
        curVertId = nodeIdtoVertIdMap[node_p]
        append!(nodeIdtoVertIdMap, Int64.(ones(numSubPartns)*curVertId) )
    end
    
    grParams = merge(grParams, (nodeCoordsList=nodeCoordsListNew,))
    return grParams, newNodesRange, nodeIdtoVertIdMap
end
function AddstartGoalToGrPosSmplUB(baseGrParams, startgoalCoords)
    grPosSmpl = baseGrParams.grPosSmpl
    smplGrEdgType = baseGrParams.smplGrEdgType
    vertCoordsList = baseGrParams.vertCoordsList
    sidesList = baseGrParams.sidesList
    sideIdsList = baseGrParams.sideIdsList
    vertCoordsList[1,1:2] = startgoalCoords[1]
    vertCoordsList[2,1:2] = startgoalCoords[2]
    vertCoordsList[end,1:2] = startgoalCoords[2]

    for k=3:size(vertCoordsList,1)
        # sidesListTrimmed = copy(sidesList)
        # sideIdP2 = sideIdsList[k]
        # deleteat!(sidesListTrimmed,sideIdP2)
        if ~CheckLineIntersectsSidesStrtGl(startgoalCoords[1], vertCoordsList[k,1:2], sidesList)
            LightGraphs.add_edge!(grPosSmpl,1,k)  
            # set_prop!(grPosSmpl, 1, k, :edgType, "Unrestricted")
            smplGrEdgType[1,k] = 3.
        end
    end
    for k=1:size(vertCoordsList,1)-1
        if k!= 2
            if ~CheckLineIntersectsSidesStrtGl(vertCoordsList[k,1:2], startgoalCoords[2], sidesList)
                LightGraphs.add_edge!(grPosSmpl,k, size(vertCoordsList,1))  
                # set_prop!(grPosSmpl, k, size(vertCoordsList,1), :edgType, "Unrestricted")
                smplGrEdgType[k,size(vertCoordsList,1)] = 3.
            end
        end
    end
    baseGrParams = merge(baseGrParams, (grPosSmpl=grPosSmpl,))
    return baseGrParams
end
function AddStartGoalToGraphUB(gr_uam, edgWtsMat, nodeCoordsList, vertCoordsList, startgoalCoords, chrgSamplesVec, batteryParams, sidesList)

    alpha = batteryParams.discharge_rate
    beta = batteryParams.recharge_rate
    numNodes = size(nodeCoordsList,1)
    numPartnsGoal = length(chrgSamplesVec)
    
    # vertCoordsList[1,1:2] = startgoalCoords[1]
    # vertCoordsList[2,1:2] = startgoalCoords[2]
    # vertCoordsList[end,1:2] = startgoalCoords[2]
    nodeCoordsList[1,1:2] = startgoalCoords[1]
    nodeCoordsList[2,1:2] = startgoalCoords[2]
    [nodeCoordsList[end-j+1,1:2] = startgoalCoords[2] for j=1:numPartnsGoal]
    # @show startgoalCoords[2]
    # @show numPartnsGoal
    # @show nodeCoordsList
    for k = 3:numNodes
        if ~CheckLineIntersectsSidesStrtGl(nodeCoordsList[1,1:2], nodeCoordsList[k,1:2], sidesList)
            batFeas, switchPt=  CheckBatteryFeasibilityUB(nodeCoordsList[1,:], nodeCoordsList[k,:], batteryParams)
            if batFeas   
                gasEdgLen = nodeCoordsList[1,1:2] == switchPt ? 0.00001 : norm(switchPt[1:2]-nodeCoordsList[1,1:2] )
                LightGraphs.add_edge!(gr_uam,1,k)                            
                edgWtsMat[1, k] = gasEdgLen                                         
            end
        end
    end

    for j = 3:numNodes-numPartnsGoal
        for k=numNodes-numPartnsGoal+1:numNodes
            if ~CheckLineIntersectsSidesStrtGl(nodeCoordsList[j,1:2], nodeCoordsList[k,1:2], sidesList)
                batFeas, switchPt=  CheckBatteryFeasibilityUB(nodeCoordsList[j,:], nodeCoordsList[k,:], batteryParams)
                if batFeas    
                    gasEdgLen = nodeCoordsList[j,1:2] == switchPt ? 0.00001 : norm(switchPt[1:2]-nodeCoordsList[j,1:2] )                                   
                    LightGraphs.add_edge!(gr_uam,j,k)                            
                    edgWtsMat[j, k] = gasEdgLen                                         
                end
            end
        end
    end

    return(gr_uam=gr_uam, edgWtsMat=edgWtsMat, nodeCoordsList=nodeCoordsList, vertCoordsList=vertCoordsList)
end


function PosSmplgGraphUB(vertCoordsList, obsIdList, sidesList, sideIdsList, noFlyInds=Vector{Int64}(undef,0))
# sideIdsList is the sideId of the vertices in the vertCoordsList

    grSmpl = LightGraphs.DiGraph(size(vertCoordsList,1))
    
    goalObsId = obsIdList[end]
    smplGrEdgType = sparse(zeros(size(vertCoordsList,1),size(vertCoordsList,1) ) )
    # edgType, LastLeg:1, QuietZone:2, Unrestricted:3
    smplGrEdgType[LightGraphs.nv(grSmpl),2] = 1.

    LightGraphs.add_edge!(grSmpl,LightGraphs.nv(grSmpl),2) #edge from last node (goal replica) to goal node
    

    for i=3:LightGraphs.nv(grSmpl)-1
        for j=3:LightGraphs.nv(grSmpl)-1
            p1 = vertCoordsList[i,1:2]
            p2 = vertCoordsList[j,1:2]        
            
            if i!=j && obsIdList[i] == obsIdList[j] && p1 != p2 # quiet zone edges 
                # these are the edges connecting vertices of same obstacle

                
                sideIdsP1 = sideIdsList[i,:]
                sideIdsP2 = sideIdsList[j,:]
                commonSideId = intersect(sideIdsP1, sideIdsP2)
                if length(commonSideId)>0
                    LightGraphs.add_edge!(grSmpl,i,j)
                    smplGrEdgType[i,j] = 2.
                elseif ~(obsIdList[i] in noFlyInds)
                    LightGraphs.add_edge!(grSmpl,i,j)                
                    smplGrEdgType[i,j] = 2.
                end
            elseif i!=j && obsIdList[i] != 0 && obsIdList[j]>0
                sidesListTrimmed = copy(sidesList)
                sideIdsP1P2 = sort([sideIdsList[i,:]; sideIdsList[j,:]])
                sideIdsP1P2 = sideIdsP1P2[findall(sideIdsP1P2.!=0)]
                deleteat!(sidesListTrimmed,sideIdsP1P2)
                if ~CheckLineIntersectsSides2(p1,p2,sidesListTrimmed)
                    LightGraphs.add_edge!(grSmpl,i,j)                      
                    smplGrEdgType[i,j] = 3.

                end
            end
        end
    end
    return grSmpl, smplGrEdgType
end

function SmplgGraphUB(grPosSmpl, smplGrEdgType, vertCoordsList, obsIdList, batteryParams, numDiscs, flyOverQZ=true)
    # There are placeholders for nodes corresponding to start, goal, and goal replicas 
    # The edges are not added between start, goal and other nodes
    # Zero cost edges between goal replicas and goal are added

    nodeCoordsList = Matrix{Float64}(undef,0,3)
    nodeObsIdList = Vector{Int64}(undef,0)
    nodeIdtoVertIdMap = Vector{Int64}(undef,0)
    nodesOfVertList = Vector{Int64}[] # each element V_i is a vector of nodes corresponding to the vertex i
    lastNodeId = 0
    chrgSamplesVec = LinRange(batteryParams.qmin, batteryParams.qmax,numDiscs)
    for k = 1:size(vertCoordsList,1)
        if k<=2
            nodeCoordsList = [nodeCoordsList; [vertCoordsList[k,1] vertCoordsList[k,2] 100 ] ]  
                
            push!(nodeObsIdList, obsIdList[k])
            push!(nodeIdtoVertIdMap, k)
            push!(nodesOfVertList, [lastNodeId+1] ) 
        else
            for i=1:length(chrgSamplesVec)
                nodeCoordsList = [nodeCoordsList; [vertCoordsList[k,1] vertCoordsList[k,2] chrgSamplesVec[i] ] ]  
                
                push!(nodeObsIdList, obsIdList[k])
                push!(nodeIdtoVertIdMap, k)
            end   
            push!(nodesOfVertList, collect(lastNodeId+1:lastNodeId+length(chrgSamplesVec) ) )    
        end
        
        lastNodeId = nodesOfVertList[end][end]   
        
    end
    numNodes = size(nodeCoordsList,1)
    
    gr_uam = LightGraphs.SimpleDiGraph(numNodes) 

    edgWtsMat = zeros(numNodes,numNodes)
    for edg in LightGraphs.edges(grPosSmpl)
        u = edg.src
        v = edg.dst
        nodes_vert_u = nodesOfVertList[u]
        nodes_vert_v = nodesOfVertList[v]        
        # edgTp = get_prop(grPosSmpl, u, v, :edgType)
        edgTp = smplGrEdgType[u,v]
        for node_u in nodes_vert_u
            for node_v in nodes_vert_v
                gr_uam, edgWtsMat = AddEdgeToGrUB(gr_uam, edgWtsMat, node_u, node_v, edgTp, nodeCoordsList, batteryParams, flyOverQZ)
            end
        end
    end    
    edgWtsMat = sparse(edgWtsMat)

    return (gr_uam=gr_uam, edgWtsMat=edgWtsMat, vertCoordsList=vertCoordsList, nodeCoordsList=nodeCoordsList, chrgSamplesVec=chrgSamplesVec, nodesOfVertList=nodesOfVertList, nodeIdtoVertIdMap=nodeIdtoVertIdMap)
end


function AddEdgeToGrUB(gr_uam, edgWtsMat, node_u, node_v, edgTp, nodeCoordsList, batteryParams, flyOverQZ)

    if edgTp == 1.
        if nodeCoordsList[node_u,3] >= batteryParams.min_charge_at_dest
            LightGraphs.add_edge!(gr_uam,node_u,node_v)                           
            edgWtsMat[node_u, node_v] = 0.00001                  
        end
    else
        batFeas, switchPt=  CheckBatteryFeasibilityUB(nodeCoordsList[node_u,:], nodeCoordsList[node_v,:], batteryParams)
        if batFeas    
            gasEdgLen = norm(switchPt-nodeCoordsList[node_u,1:2] )
            if edgTp == 3.
                LightGraphs.add_edge!(gr_uam,node_u,node_v)  
                gasEdgLen = gasEdgLen==0 ? 0.00001 : gasEdgLen
                edgWtsMat[node_u, node_v] = gasEdgLen                  
            elseif edgTp == 2. && gasEdgLen<1e-6 #quiet zone edges
                LightGraphs.add_edge!(gr_uam,node_u,node_v)    
                if flyOverQZ                                        
                    edgWtsMat[node_u, node_v] = 0.00001
                else                                    
                    edgWtsMat[node_u, node_v] = 1e8
                end

            end
        end
    end
    return gr_uam, edgWtsMat
end

function PlotGraphNodesUB(vertCoordsList)
    

    PyPlot.scatter(vertCoordsList[1:2,1], vertCoordsList[1:2,2], marker="x",color="blue")
    PyPlot.scatter(vertCoordsList[3:end,1], vertCoordsList[3:end,2], marker=".",color="blue")

    PyPlot.text(vertCoordsList[1,1], vertCoordsList[1,2]-400, "s" )
    PyPlot.text(vertCoordsList[2,1], vertCoordsList[2,2]-400, "t" )
  

end


function CheckBatteryFeasibilityUB(node1, node2, batteryParams)
    edgLen = norm(node1[1:2]-node2[1:2])

    if edgLen < 1e-6        
        return (false, node1[1:2])
    end
    b1 = node1[3]
    b2 = node2[3]

    feasCheck, lamdaVec = CheckChrgLvlsFeas(b1,b2,batteryParams, edgLen)
    if length(lamdaVec) == 1
        lamda = lamdaVec[1]        
    else
        lamda = lamdaVec[1]+lamdaVec[3]
    end
    switchPoint = (1-lamda)*node1[1:2] + lamda[1]*node2[1:2]
    return (feasCheck, switchPoint)
end

function CheckBatteryFeasibilityUB_chrgprf(node1, node2, batteryParams)
    edgLen = norm(node1[1:2]-node2[1:2])

    if edgLen < 1e-6        
        return (false, Matrix(node1[1:2]') ) 
    end
    b1 = node1[3]
    b2 = node2[3]

    feasCheck, lambdaVec = CheckChrgLvlsFeas(b1,b2,batteryParams, edgLen)
    if length(lambdaVec) == 1
        switchPoints = (1-lambdaVec[1])*node1[1:2] + lambdaVec[1]*node2[1:2]
        switchPoints = Matrix(switchPoints')
    else
        lambdaVec = cumsum(lambdaVec)
        switchPoints = (1.0 .-lambdaVec).*node1[1:2]' + lambdaVec.*node2[1:2]'
    end
    return (feasCheck, switchPoints)
end

function PlotChargeProfileUB(sp, nodeCoordsList, batteryParams)
    distTrvld = 0.
    curChrg = 100.
    distVsChrg = [distTrvld curChrg]
    edgDists = [0]
    alpha = batteryParams.discharge_rate
    beta = batteryParams.recharge_rate
    for edg in sp
    
        node_src = nodeCoordsList[edg.src,:]
        node_dst = nodeCoordsList[edg.dst,:]
        # edgDists = [edgDists; edgDists[end]+norm(node_src[1:2]-node_dst[1:2])]
        edgLen = norm(node_src[1:2]-node_dst[1:2])
        edgDists = [edgDists; edgLen]
        q_src = nodeCoordsList[edg.src,3]
        q_dst = nodeCoordsList[edg.dst,3]
        feas, switchPoints = CheckBatteryFeasibilityUB_chrgprf(node_src, node_dst, batteryParams)
        prevSwitchPt = node_src[1:2]
        curChrg = q_src
        for k = 1:size(switchPoints,1)
            
            # @show switchPoints          
            switchPoint = switchPoints[k,1:2]  
            distVsChrg = [distVsChrg; distTrvld curChrg]
            distTrvld += norm(prevSwitchPt-switchPoint)
            if isodd(k)
                curChrg = curChrg+norm(prevSwitchPt-switchPoint)*beta
                curChrg = curChrg>100. ? 100. : curChrg
            else
                curChrg = curChrg-norm(prevSwitchPt-switchPoint)*alpha
            end
            distVsChrg = [distVsChrg; distTrvld curChrg]    


            prevSwitchPt = switchPoint
        end
        
        distTrvld += norm(switchPoints[end,:]-node_dst[1:2])
        curChrg = curChrg-norm(switchPoints[end,:]-node_dst[1:2])*alpha
        distVsChrg = [distVsChrg; distTrvld q_dst]
    end
    distVsChrg = distVsChrg[1:end-1,:]
    PyPlot.plot(distVsChrg[:,1], distVsChrg[:,2], linewidth=2, color="blue")
    edgDists = cumsum(edgDists)
    PyPlot.vlines(edgDists,0,100, color="blue", linewidth=1, linestyle="--")
end