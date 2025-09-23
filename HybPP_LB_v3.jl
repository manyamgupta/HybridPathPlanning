#v3 includes no-flyzones

function HppBaseGrLB(obsList, delL, batteryParams, chargePartns, noFlyInds=Vector{Int64}(undef,0))
    
    vertCoordsList, obsIdList, sidesList, sideIdsList = DiscPosVerts(obsList, delL, (NaN, NaN), (NaN, NaN), noFlyInds)
    numVerts = size(vertCoordsList,1)
    chrgPrtnsVec = [chargePartns for i=1:numVerts-2]
    insert!(chrgPrtnsVec,1,[(100.,100.)])
    insert!(chrgPrtnsVec,1,[(100.,100.)])
    
    grPosSmpl, smplGrEdgType = PosSmplgGraphLB(vertCoordsList, obsIdList, sidesList, sideIdsList, noFlyInds)
    
    baseGrParams = SmplgGraphLB(grPosSmpl, smplGrEdgType, vertCoordsList, obsIdList, batteryParams, chrgPrtnsVec)

    baseGrParams = merge(baseGrParams, (sidesList=sidesList, grPosSmpl=grPosSmpl, smplGrEdgType=smplGrEdgType))
    
    return baseGrParams
    
end

function HppLB(baseGrParams, start, goal, batteryParams, numSubPartns)

    nodeCoordsList = baseGrParams.nodeCoordsList
    edgWtsMat = baseGrParams.edgWtsMat
    gr_uam = baseGrParams.gr_uam
    chrgPrtnsVec = baseGrParams.chrgPrtnsVec
    sidesList = baseGrParams.sidesList
    nodesOfVertList = baseGrParams.nodesOfVertList
    nodeIdtoVertIdMap = baseGrParams.nodeIdtoVertIdMap

    sp = LightGraphs.SimpleGraphs.SimpleEdge{Int64}[]
    pathCostLB=0.
    nodesToBePatnd = Int64[]
    nodesRemoved = Vector{Int64}(undef,0)
    grParams = Tuple{}
    # tic = time()
    baseGrParams = AddstartGoalToGrPosSmplLB(baseGrParams, (start, goal), sidesList)
    # compTime_AddSGPosGr = time()-tic
    
    # tic = time()
    grParams = AddStartGoalToGraphLB(gr_uam, baseGrParams, edgWtsMat, nodeCoordsList, (start, goal), chrgPrtnsVec, batteryParams, sidesList)
    # compTime_AddSGGrLB = time()-tic

    # @show compTime_AddSGPosGr
    # @show compTime_AddSGGrLB

    for k=1:1
        if k>=2
            # chrgPrtnsVec, nodesOfVertList = IterPartnChrgLvls(nodesToBePatnd, chrgPrtnsVec,numSubPartns, nodesOfVertList, nodeIdtoVertIdMap )
            grParams, newNodesRange, nodeIdtoVertIdMap = AddNewNodesToListLB(grParams, nodesToBePatnd, numSubPartns, nodeIdtoVertIdMap )
            
            # tic = time()
            grParams = AddNewEdgesLB(grParams, baseGrParams, newNodesRange, batteryParams, nodesRemoved, nodeIdtoVertIdMap)
            # compTime_AddEdgs = time()-tic
            # @show compTime_AddEdgs
            # tic = time()
            grParams = RemoveEdgesLB(grParams, nodesToBePatnd)
            # compTime_RemEdgs = time()-tic
            # @show compTime_RemEdgs
            # return (pathCostLB=0, sp=[], chrgPrtnsVec=chrgPrtnsVec, grParams=grParams,newNodesRange=newNodesRange,nodeIdtoVertIdMap=nodeIdtoVertIdMap)
        end
        
        sp = LightGraphs.a_star(grParams.gr_uam, 1, 2, sparse(grParams.edgWtsMat))
        
        
        pathCostLB = ComputePathCost(sp, grParams.edgWtsMat)
        
        # @show pathCostLB
        nodesToBePatnd = IntermediateNodes(sp)
        append!(nodesRemoved, nodesToBePatnd)
    end
    return (pathCostLB=pathCostLB, sp=sp, chrgPrtnsVec=chrgPrtnsVec, grParams=grParams)
end

function RemoveEdgesLB(grParams, nodesToBePatnd)
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

function AddNewEdgesLB(grParams, baseGrParams, newNodesRange, batteryParams, nodesRemoved, nodeIdtoVertIdMap)

    # vertCoordsList = grParams.vertCoordsList

    # @show nodeCoordsListNew
    # @show vertCoordsList
    # @show nodesOfVertList
    # @show nodeIdtoVertIdMap
    grPosSmpl = baseGrParams.grPosSmpl
    smplGrEdgType = baseGrParams.smplGrEdgType
    gr_uam = grParams.gr_uam
    nodeCoordsList = grParams.nodeCoordsList
    edgWtsMat = grParams.edgWtsMat
    nv_old = LightGraphs.nv(gr_uam)    
    numNewNodes = length(newNodesRange)
    nv_new = nv_old+numNewNodes
    LightGraphs.add_vertices!(gr_uam,numNewNodes)
    edgWtsMat = [[edgWtsMat zeros(nv_old, numNewNodes)];zeros(numNewNodes, nv_old+numNewNodes)]
    tic = time()
    count_CallsAddEdg = 0
    # for node_u in 1:nv_old
    #     if ~(node_u in nodesRemoved)
    #         for node_v in newNodesRange
    #             if LightGraphs.has_edge(grPosSmpl, nodeIdtoVertIdMap[node_u], nodeIdtoVertIdMap[node_v])
                    
    #                 edgTp = smplGrEdgType[nodeIdtoVertIdMap[node_u], nodeIdtoVertIdMap[node_v]]
    #                 AddEdgeToGrLB(gr_uam, edgWtsMat, node_u, node_v, edgTp, nodeCoordsList, batteryParams)
    #                 count_CallsAddEdg += 1
    #             end
    #             if LightGraphs.has_edge(grPosSmpl, nodeIdtoVertIdMap[node_v], nodeIdtoVertIdMap[node_u])
                    
    #                 edgTp = smplGrEdgType[nodeIdtoVertIdMap[node_v], nodeIdtoVertIdMap[node_u]]
    #                 AddEdgeToGrLB(gr_uam, edgWtsMat, node_v, node_u, edgTp, nodeCoordsList, batteryParams)
    #                 count_CallsAddEdg += 1
    #             end
    #         end
    #     end
    # end
    # time_addnewedges1 = time()-tic

    # tic = time()
    # for node_u in newNodesRange        
    #     for node_v in newNodesRange
    #         if LightGraphs.has_edge(grPosSmpl, nodeIdtoVertIdMap[node_u], nodeIdtoVertIdMap[node_v])
                
    #             edgTp = smplGrEdgType[nodeIdtoVertIdMap[node_u], nodeIdtoVertIdMap[node_v]]
    #             AddEdgeToGrLB(gr_uam, edgWtsMat, node_u, node_v, edgTp, nodeCoordsList, batteryParams)
    #             count_CallsAddEdg += 1
    #         end
    #     end        
    # end
    # time_addnewedges2 = time()-tic

    tic = time()
    newNodes = setdiff(1:nv_new, nodesRemoved)
    for node_u in newNodes
        for node_v in newNodes
            if ~LightGraphs.has_edge(gr_uam, node_u, node_v)
                if LightGraphs.has_edge(grPosSmpl, nodeIdtoVertIdMap[node_u], nodeIdtoVertIdMap[node_v])
                    
                    edgTp = smplGrEdgType[nodeIdtoVertIdMap[node_u], nodeIdtoVertIdMap[node_v]]
                    AddEdgeToGrLB(gr_uam, edgWtsMat, node_u, node_v, edgTp, nodeCoordsList, batteryParams)
                    count_CallsAddEdg += 1
                end
            end

        end
    end
    time_AddNewEdgesLB = time()-tic
    # @show time_AddNewEdgesLB
    
    grParams = merge(grParams, (gr_uam=gr_uam, edgWtsMat=edgWtsMat))
    
    # @show count_CallsAddEdg
    return grParams
end

function AddNewNodesToListLB(grParams, nodesToBePatnd, numSubPartns, nodeIdtoVertIdMap )
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
function AddstartGoalToGrPosSmplLB(baseGrParams, startgoalCoords, sidesList)
    grPosSmpl = baseGrParams.grPosSmpl
    smplGrEdgType = baseGrParams.smplGrEdgType
    vertCoordsList = baseGrParams.vertCoordsList
    vertCoordsList[1,1:2] = startgoalCoords[1]
    vertCoordsList[2,1:2] = startgoalCoords[2]
    vertCoordsList[end,1:2] = startgoalCoords[2]

    for k=3:size(vertCoordsList,1)
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
    baseGrParams = merge(baseGrParams, (grPosSmpl=grPosSmpl,smplGrEdgType=smplGrEdgType))
    return baseGrParams
end
function AddStartGoalToGraphLB(gr_uam, baseGrParams, edgWtsMat, nodeCoordsList,  startgoalCoords, chrgPrtnsVec, batteryParams, sidesList)

    alpha = batteryParams.discharge_rate
    beta = batteryParams.recharge_rate
    numNodes = size(nodeCoordsList,1)
    numPartnsGoal = length(chrgPrtnsVec[end])
    vertCoordsList = baseGrParams.vertCoordsList
    nodeIdtoVertIdMap = baseGrParams.nodeIdtoVertIdMap
    
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
        # if LightGraphs.has_edge(baseGrParams.grPosSmpl, nodeIdtoVertIdMap[1], nodeIdtoVertIdMap[k] )
            batFeas, switchPt=  CheckBatteryFeasibilityLB(nodeCoordsList[1,:], nodeCoordsList[k,:], batteryParams)
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
            # if LightGraphs.has_edge(baseGrParams.grPosSmpl, nodeIdtoVertIdMap[j], nodeIdtoVertIdMap[k] )
                batFeas, switchPt=  CheckBatteryFeasibilityLB(nodeCoordsList[j,:], nodeCoordsList[k,:], batteryParams)
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


function PosSmplgGraphLB(vertCoordsList, obsIdList, sidesList, sideIdsList, noFlyInds)

    grSmpl = LightGraphs.DiGraph(size(vertCoordsList,1))
    # mgrSmpl = MetaDiGraph(grSmpl,0.)
    mgrSmpl = grSmpl
    smplGrEdgType = sparse(zeros(size(vertCoordsList,1),size(vertCoordsList,1) ) )
    # edgType, LastLeg:1, QuietZone:2, Unrestricted:3
    goalObsId = obsIdList[end]
    LightGraphs.add_edge!(mgrSmpl,LightGraphs.nv(grSmpl),2)
    
    # set_prop!(mgrSmpl,Graphs.nv(grSmpl),2, :edgType, "LastLeg")
    smplGrEdgType[LightGraphs.nv(grSmpl),2] = 1.
    for i=3:LightGraphs.nv(grSmpl)-1
        for j=3:LightGraphs.nv(grSmpl)-1
            p1 = vertCoordsList[i,1:2]
            p2 = vertCoordsList[j,1:2]        
            # if obsIdList[i] == goalObsId && obsIdList[j] == 0 # edge from dummy goal to goal
            #     Graphs.add_edge!(mgrSmpl,i,j)
            #     set_prop!(mgrSmpl,i,j, :edgType, "LastLeg")
            # elseif obsIdList[i] == goalObsId # no edges from dummy goal to other nodes
            #     continue                
            if i!=j && obsIdList[i] == obsIdList[j] && p1 != p2
                # quiet zone edges 
                # these are the edges connecting vertices of same obstacle
                # LightGraphs.add_edge!(mgrSmpl,i,j)                
                # smplGrEdgType[i,j] = 2.
                sideIdsP1 = sideIdsList[i,:]
                sideIdsP2 = sideIdsList[j,:]
                commonSideId = intersect(sideIdsP1, sideIdsP2)
                if length(commonSideId)>0
                    LightGraphs.add_edge!(mgrSmpl,i,j)
                    smplGrEdgType[i,j] = 2.
                elseif ~(obsIdList[i] in noFlyInds)
                    LightGraphs.add_edge!(mgrSmpl,i,j)                
                    smplGrEdgType[i,j] = 2.
                end

                                
   
            elseif i!=j && obsIdList[i] != 0 && obsIdList[j]>0
                sidesListTrimmed = copy(sidesList)                
                sideIdsP1P2 = sort([sideIdsList[i,:]; sideIdsList[j,:]])
                sideIdsP1P2 = sideIdsP1P2[findall(sideIdsP1P2.!=0)]
                deleteat!(sidesListTrimmed,sideIdsP1P2)

                if ~CheckLineIntersectsSides2(p1,p2,sidesListTrimmed)
                    LightGraphs.add_edge!(mgrSmpl,i,j)  
                    # set_prop!(mgrSmpl,i,j, :edgType, "Unrestricted")
                    smplGrEdgType[i,j] = 3.

                end
            end
        end
    end
    return mgrSmpl, smplGrEdgType
end

function SmplgGraphLB(grPosSmpl, smplGrEdgType, vertCoordsList, obsIdList, batteryParams, chrgPrtnsVec)
    # There are placeholders for nodes corresponding to start, goal, and goal replicas 
    # The edges are not added between start, goal and other nodes
    # Zero cost edges between goal replicas and goal are added

    nodeCoordsList = Matrix{Float64}(undef,0,4)
    nodeObsIdList = Vector{Int64}(undef,0)
    nodeIdtoVertIdMap = Vector{Int64}(undef,0)
    nodesOfVertList = Vector{Int64}[] # each element V_i is a vector of nodes corresponding to the vertex i
    alpha = batteryParams.discharge_rate
    beta = batteryParams.recharge_rate
    lastNodeId = 0
    for k = 1:size(vertCoordsList,1)
      
        for i=1:length(chrgPrtnsVec[k])
            nodeCoordsList = [nodeCoordsList; [vertCoordsList[k,1] vertCoordsList[k,2] chrgPrtnsVec[k][i][1] chrgPrtnsVec[k][i][2]] ]  
            
            push!(nodeObsIdList, obsIdList[k])
            push!(nodeIdtoVertIdMap, k)
        end   
        push!(nodesOfVertList, collect(lastNodeId+1:lastNodeId+length(chrgPrtnsVec[k]) ) )
        lastNodeId = nodesOfVertList[end][end]     
    end
    numNodes = size(nodeCoordsList,1)
    goalObsId = obsIdList[end]
    # @show numNodes
    # @show nodesOfVertList

    # gr_temp = DiGraph(numNodes) 
    # gr_uam = MetaDiGraph(gr_temp, 0.0)
    gr_uam = LightGraphs.SimpleDiGraph(numNodes) 

    edgWtsMat = zeros(numNodes,numNodes)
    tic = time()
    count_callAddEdgs_FirstGr = 0
    for edg in LightGraphs.edges(grPosSmpl)
        u = edg.src
        v = edg.dst
        nodes_vert_u = nodesOfVertList[u]
        nodes_vert_v = nodesOfVertList[v]
        # edgTp = get_prop(grPosSmpl, u, v, :edgType)
        edgTp = smplGrEdgType[u,v]
        for node_u in nodes_vert_u
            for node_v in nodes_vert_v
                AddEdgeToGrLB(gr_uam, edgWtsMat, node_u, node_v, edgTp, nodeCoordsList, batteryParams)
                count_callAddEdgs_FirstGr += 1
            end
        end
    end    
    edgWtsMat = sparse(edgWtsMat)
    compTime_AddEdges_FirstGr = time()-tic
    # @show compTime_AddEdges_FirstGr
    # @show count_callAddEdgs_FirstGr
    
    return (gr_uam=gr_uam, edgWtsMat=edgWtsMat, vertCoordsList=vertCoordsList, nodeCoordsList=nodeCoordsList, chrgPrtnsVec=chrgPrtnsVec, nodesOfVertList=nodesOfVertList, nodeIdtoVertIdMap=nodeIdtoVertIdMap)
end


function AddEdgeToGrLB(gr_uam, edgWtsMat, node_u, node_v, edgTp, nodeCoordsList, batteryParams)

    if edgTp == 1.
        if nodeCoordsList[node_u,4] >= batteryParams.min_charge_at_dest
            LightGraphs.add_edge!(gr_uam,node_u,node_v)                           
            edgWtsMat[node_u, node_v] = 0.00001                  
        end
    else
        batFeas, switchPt=  CheckBatteryFeasibilityLB(nodeCoordsList[node_u,:], nodeCoordsList[node_v,:], batteryParams)
        if batFeas    
            gasEdgLen = norm(switchPt-nodeCoordsList[node_u,1:2] )
            if edgTp == 3.
                LightGraphs.add_edge!(gr_uam,node_u,node_v)  
                gasEdgLen = gasEdgLen==0 ? 0.00001 : gasEdgLen
                edgWtsMat[node_u, node_v] = gasEdgLen                  
            elseif edgTp == 2. && gasEdgLen<1e-6
                LightGraphs.add_edge!(gr_uam,node_u,node_v)                            

                edgWtsMat[node_u, node_v] = 0.00001                                                                           
            end
        end
    end
    return
end


function CheckVertsOnSameSide(p1, p2, side1, side2, )
    if SlopeOfSide((p1,p2)) == SlopeOfSide(side1) || SlopeOfSide((p1,p2)) == SlopeOfSide(side2) || SlopeOfSide((p2,p1)) == SlopeOfSide(side1) || SlopeOfSide((p2,p1)) == SlopeOfSide(side2)
        return true
    else
        return false
    end
end

function SlopeOfSide(side)

    return round((side[2][2]-side[1][2])/(side[2][1]-side[1][1]), digits=2)

end


function IterPartnChrgLvls(nodesToBePatnd, chrgPrtnsVec, numSubPartns, nodesOfVertList, nodeIdtoVertIdMap)

    # @show chrgPrtnsVec
    # @show nodesOfVertList
    # @show nodesToBePatnd
    # @show nodeIdtoVertIdMap
       
    chrgPrtnsVecNew = copy(chrgPrtnsVec)
    numVertsCumSum = cumsum(length.(chrgPrtnsVec))
    lastNodeId = nodesOfVertList[end][end] 
    newNodesRange = lastNodeId+1:lastNodeId+numSubPartns*length(nodesToBePatnd)
    for nodePrtn in nodesToBePatnd
        
        rowInd, elemIndInRow = FindNodePartnIndex(nodePrtn, numVertsCumSum)
        rowToBeSubPartnd = copy(chrgPrtnsVec[rowInd])
        nodesRowSubPrtd = copy(nodesOfVertList[rowInd])
        newNodesRow = [nodesRowSubPrtd;collect(lastNodeId+1:lastNodeId+numSubPartns)]
        lastNodeId += numSubPartns
        deleteat!(chrgPrtnsVecNew,rowInd)
        deleteat!(nodesOfVertList, rowInd)
       
        chrgIntervalSubPartnd = rowToBeSubPartnd[elemIndInRow]
        delQ = (chrgIntervalSubPartnd[2]-chrgIntervalSubPartnd[1])/numSubPartns
        newIntervals = [(chrgIntervalSubPartnd[1] + (i-1)*delQ, chrgIntervalSubPartnd[1] + i*delQ) for i=1:numSubPartns]

        rowToBeSubPartnd = [rowToBeSubPartnd; newIntervals]
        insert!(chrgPrtnsVecNew, rowInd, rowToBeSubPartnd)
        insert!(nodesOfVertList, rowInd, newNodesRow)
        [push!(nodeIdtoVertIdMap, rowInd) for i=1:numSubPartns]

    end
    # @show chrgPrtnsVecNew
    # @show nodesOfVertList
    # @show nodeIdtoVertIdMap
    return chrgPrtnsVecNew, nodesOfVertList, nodeIdtoVertIdMap, newNodesRange

end

function FindNodePartnIndex(nodePrtn, numVertsCumSum)

    rowInd = findfirst(x->x>=nodePrtn, numVertsCumSum)
    elemInd = nodePrtn-numVertsCumSum[rowInd-1]
    return rowInd, elemInd
end

function IntermediateNodes(sp)
    
    nodesList = Vector{Int64}(undef,0)
    for edg in sp[1:end-1]
        append!(nodesList, edg.dst)
    end
    return nodesList
end



function ConstructUBFromLB(sp, nodeCoordsList, alpha, beta)

    # pathNodes = [edg.src for edg in sp]
    # push!(pathNodes, sp[end].dst )
    # nodeCoordsList = nodeCoordsList[pathNodes, :]
    
    pathCost = 0.
    chrg_v = nodeCoordsList[sp[end-1].dst,4]
    chrgVec = [chrg_v]
    for i = 1:length(sp)-1        
        u = sp[end-i].src
        v = sp[end-i].dst
        nodeCoords_u = nodeCoordsList[u,1:2]
        nodeCoords_v = nodeCoordsList[v,1:2]        
        feasCheck = true
        chrg_u_cands = [0.5*sum(nodeCoordsList[u,3:4]), nodeCoordsList[u,4], nodeCoordsList[u,3]]
        for k=1:3
            chrg_u = chrg_u_cands[k]
            feasCheck, switchPoint = CheckBatteryFeasibility(push!(nodeCoords_u[1:2], chrg_u), push!(nodeCoords_v[1:2], chrg_v), alpha, beta)        
            if feasCheck
                pathCost += norm(nodeCoordsList[u,1:2]-switchPoint)
                chrg_v = chrg_u
                push!(chrgVec, chrg_u)
                break
            end
        end
        if !feasCheck
            println("Charge levels are not feasible")
            return NaN, nodeCoordsListUB
        end            
    end

    return pathCost, reverse(chrgVec)
end

function PlotGraphNodesLB(vertCoordsList, chrgPrtnsVec=[])
    numPartns = length.(chrgPrtnsVec)

    PyPlot.scatter(vertCoordsList[1:2,1], vertCoordsList[1:2,2], marker="x",color="blue")
    PyPlot.scatter(vertCoordsList[3:end,1], vertCoordsList[3:end,2], marker=".",color="blue")

    PyPlot.text(vertCoordsList[1,1], vertCoordsList[1,2]-400, "s" )
    PyPlot.text(vertCoordsList[2,1], vertCoordsList[2,2]-400, "t" )

    # for i = 1:size(vertCoordsList,1)
    #     PyPlot.text(vertCoordsList[i,1]+200, vertCoordsList[i,2]+200,"$(numPartns[i])" )
    # end

end

function PlotChargeProfileLB(sp, nodeCoordsList, batteryParams)
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
        q_src = nodeCoordsList[edg.src,4]
        q_dst = nodeCoordsList[edg.dst,3]
        feas, switchPoints = CheckBatteryFeasibilityLB_chrgprf(node_src, node_dst, batteryParams)
        prevSwitchPt = node_src[1:2]
        curChrg = q_src

        for k = 1:size(switchPoints,1)
            
                 
            switchPoint = switchPoints[k,1:2]  
            # @show switchPoint
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
        # distVsChrg = [distVsChrg; distTrvld q_src]
        # distTrvld += norm(node_src[1:2]-switchPoint)
        # curChrg = q_src+norm(node_src[1:2]-switchPoint)*beta
        # curChrg = curChrg>100. ? 100. : curChrg        
        # distVsChrg = [distVsChrg; distTrvld curChrg]
    
        distTrvld += norm(switchPoints[end,:]-node_dst[1:2])
        curChrg = curChrg-norm(switchPoints[end,:]-node_dst[1:2])*alpha
        distVsChrg = [distVsChrg; distTrvld q_dst]
        
    end
    distVsChrg = distVsChrg[1:end-1,:]
    PyPlot.plot(distVsChrg[:,1], distVsChrg[:,2], linewidth=2, color="brown")
    edgDists = cumsum(edgDists)
    PyPlot.vlines(edgDists,0,100, linewidth=1, linestyle="--", color="brown")
    
end

function CheckBatteryFeasibilityLB_chrgprf(node1, node2, batteryParams)
    edgLen = norm(node1[1:2]-node2[1:2])

    if edgLen < 1e-6        
        return (false, Matrix(node1[1:2]') ) 
    end
    b1 = node1[4]
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

function CheckBatteryFeasibilityLB(node1, node2, batteryParams)
    edgLen = norm(node1[1:2]-node2[1:2])

    if edgLen < 1e-6        
        return (false, node1[1:2])
    end
    b1 = node1[4]
    b2 = node2[3]

    feasCheck, lamdaVec = CheckChrgLvlsFeas(b1,b2,batteryParams, edgLen)
    if length(lamdaVec) == 1
        lamda = lamdaVec[1]        
    else
        lamda = lamdaVec[1]+lamdaVec[3]
    end
    switchPoint = (1-lamda)*node1[1:2] + lamda*node2[1:2]
    return (feasCheck, switchPoint)
end