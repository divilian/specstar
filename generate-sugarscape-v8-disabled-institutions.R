library(R6)
library(ineq)
library(assertthat)
library(foreach)
library(parallel)

library(doMC)
registerDoMC(detectCores())
##registerDoMC(1)


rm(list = ls()) ## clear the environment
set.seed(111) #100 set the random-number-generator's seed to allow replication

generate_sugarscape <- function(side, carcap){
    m <- matrix(numeric(side * side), nrow=side)
    blrow <- side/3
    blcol <- side/3
    trrow <- (2*side)/3
    trcol <- (2*side)/3
    
    for (rownum in (1:side)){
        for(colnum in (1:side)){
            dist.b <- sqrt((rownum-blrow)^2 + (colnum-blcol)^2)
            dist.t <- sqrt((rownum-trrow)^2 + (colnum-trcol)^2)          
            m[rownum, colnum] <- ifelse(dist.b < dist.t,
                                        carcap/dist.b^0.3,
                                        carcap/dist.t^0.3) 
        }
    }
    m[blrow, blcol] <- carcap
    m[trrow, trcol] <- carcap
    return(m)
}


generate_agents <- function(popn, initial_sug_level, metabol_range,
                            vision_range, side){
    allcoords <- expand.grid(1:side, 1:side)
    names(allcoords) <- c("x", "y")
    loccoords <- sample(1:nrow(allcoords), popn)
    
    lapply((1:popn),
           function(agnum){
               
               Agent$new(agnum, initial_sug_level,
                         ceiling(runif(1, metabol_range[1],
                                       metabol_range[2])),
                         ceiling(runif(1, vision_range[1],
                                       vision_range[2])), 
                         as.vector(t(allcoords[loccoords[agnum],]))[1],
                         as.vector(t(allcoords[loccoords[agnum],]))[2]
                         )})
}

add_new_agents <- function(lagents, sscape, birthrate, inboundrate,
                           initial_sug_level, metabol_range, vision_range){

    no.new.agents <- round(length(lagents) * birthrate) +
        round(length(lagents) * inboundrate)

    ds <- expand.grid(1:nrow(sscape), 1:nrow(sscape))
    l.allspots <- lapply(1:nrow(ds),
                         function(rownum){
                             c(ds[rownum, 1], ds[rownum,2])
                         })
    l.takenspots <- get_taken_spots(lagents)

    l.takenspots <- l.takenspots[lapply(l.takenspots,
                                        length) > 0]
    
    l.freespots <- setdiff(l.allspots, l.takenspots)
    agents.to.add <- min(no.new.agents, length(l.freespots))

    startagentid <- max(unlist(lapply(lagents, function(agent){agent$agentid}))) + 1

    
    l.prov.spots <- l.freespots[sample(1:length(l.freespots), size=agents.to.add)]
    l.newagents <- lapply((startagentid : (startagentid + agents.to.add - 1)), 
                          function(agnum){
                          Agent$new(agnum, initial_sug_level,
                                    ceiling(runif(1, metabol_range[1],
                                                  metabol_range[2])),
                                    ceiling(runif(1, vision_range[1],
                                                  vision_range[2])),
                                    l.prov.spots[[agnum - startagentid + 1]][1],
                                    l.prov.spots[[agnum - startagentid + 1]] [2])
                          }
                          )
    return(l.newagents)                           
}

get_taken_spots <- function(lst_agents){
    lapply(lst_agents, function(agent){c(agent$location.x, agent$location.y)})
}

perform_agent_motions <- function(lagents, sugarscape, v.institution.objs,
                                  threshold, timeperiod){
    lagents <- sample(lagents, length(lagents))
    lapply(lagents,
           function(agent){
               agent$move_load_consume_deposit(sugarscape, get_taken_spots(lagents),
                                       threshold, timeperiod, v.institution.objs)
               return(agent)
           }
           )
}

collate_suglevels <- function(lobjects){
    sugvals <- unlist(lapply(lobjects,
                             function(object){
                                 object$sugarlevel
                             }))
    sugvals[sugvals >= 0]
}


collate_results <- function(side, carcap, regenrt, adensity,
                            initial_sug_level, metabol_range, vision_range,
                            birthrate, inboundrate, outboundrate, threshold,
                            l.mats.suglevels, l.inst.mats.suglevels){

    tryCatch({
        vct.ginis <- sapply(l.mats.suglevels, function(mat){ineq(mat, type="Gini")})
        vct.inst.ginis <- sapply(l.inst.mats.suglevels, function(mat){ineq(mat, type="Gini")})
        c(side, carcap, regenrt, adensity, initial_sug_level, metabol_range[2], vision_range[2],
                     birthrate, inboundrate, outboundrate, threshold, vct.ginis, vct.inst.ginis)
    },
    error=function(cond){
        print(paste("error condition:", cond))
    },
    finally={}
    )
    
    
}

get_distances_df <- function(lagents){
    num.rows <- choose(length(lagents), 2)
    df <- data.frame(firstagent = rep(NA, num.rows),
                     secondagent = rep(NA, num.rows),
                     distance = rep(-1, num.rows))
    rowid <- 1

    for(agid1 in 1:(length(lagents)-1)){
        for(agid2 in (agid1 + 1):(length(lagents))){
            agent1 <- lagents[[agid1]]
            agent2 <- lagents[[agid2]]

            dist <- abs(agent1$location.x - agent2$location.x) +
                abs(agent1$location.y - agent2$location.y)

            df$firstagent[rowid] <- agent1$agentid
            df$secondagent[rowid] <- agent2$agentid
            df$distance[rowid] <- dist

            rowid <- rowid + 1
        }
    }
    return(df)
}

perform_alive_updates <- function(lagents, v.inst.objs){
    tryCatch({
        for(inst.obj in v.inst.objs){
            if(inst.obj$sugarlevel <= 0){
                print(paste("Institution object", inst.obj$institution.id,
                            "has zero sugar and has been de-activated"))
                inst.obj$alive = FALSE ## make non-alive
                ## reset all member agent's memberships
                for(aid in inst.obj$v.member.ids){
                    for(agobj in Filter(function(agent.obj){
                        agent.obj$agentid == aid}, lagents)){
                        agent.obj$institution.id <- 0
                    }
                } 
            } ## sugarlevel <= 0
            else{
                ## print(paste("institute object", inst.obj$institution.id,
                ##             "'s sugarlevel:", inst.obj$sugarlevel, collapse=""))
            }
            ## next, reset the institution object's vector of member IDs
            inst.obj$v.member.ids <- NULL
        }## main for loop
    }, ## end of try expression
    error = function(){
    },
    finally = {
        return(v.inst.objs)}
    ) ## end tryCatch
}

perform_institution_actions <- function(lagents, sugarscape, threshold,
                                        timeperiod, v.institution.objs){
    distances.df <- get_distances_df(lagents)
    ## complete the formation-of-institutions actions
    distances.df <- distances.df[distances.df$distance == 1, 1:2] ## fetch rows with interagent dist == 1
    
    tryCatch({
        if(nrow(distances.df) > 0){
            new.id <- 0
            for(rownum in (1:nrow(distances.df))){
                firstagent <- Filter(function(agent){agent$agentid == distances.df[rownum,]$firstagent},
                                     lagents)[[1]]
                secondagent <- Filter(function(agent){agent$agentid == distances.df[rownum,]$secondagent},
                                      lagents)[[1]]
                
                if(firstagent$institution.id == 0 && secondagent$institution.id == 0){
                    if(firstagent$sugarlevel > threshold && secondagent$sugarlevel > threshold){
                        if(is.null(v.institution.objs)){ ## need to create the first institution
                            new.id <- 1
                        }
                        else{
                            v.ids <- sapply(v.institution.objs,
                                            function(instobj){instobj$institution.id})
                            new.id <- max(v.ids) + 1
                        }
                        
                        ## add their IDs to the institution obj
                        ## make deposits to the institution obj's ledger - happens via new, which calls deposit()
                        institution.obj <- Institution$new(new.id,
                                                           c(firstagent$sugarlevel - threshold,
                                                             secondagent$sugarlevel - threshold),
                                                           c(firstagent$agentid, secondagent$agentid),
                                                           timeperiod)
                        ## add the institution object to the vector of institution objects
                        if(is.null(v.institution.objs)){
                            v.institution.objs <- c(institution.obj) 
                        }
                        else {
                            v.institution.objs <- c(v.institution.objs, institution.obj) 
                        }
                        ## add its ID to the two agents
                        firstagent$institution.id <- new.id
                        secondagent$institution.id <- new.id
                        ## reduce the amounts from the two agents
                        firstagent$sugarlevel <- firstagent$sugarlevel - threshold
                        secondagent$sugarlevel <- secondagent$sugarlevel - threshold
                    }
                    else{
                        ## not enough sugar in at least one of the agents,
                        ## so can't create an institution object
                    }
                } ## main if
                else if(secondagent$institution.id == 0 && firstagent$institution.id > 0 &&
                        secondagent$sugarlevel > threshold){
                    ## assign agent 2's inst. obj. id to agent 1's institution.obj.id
                    secondagent$institution.id <- firstagent$institution.id
                    institution.obj <- Filter(function(instobj){instobj$institution.id ==
                                                                    firstagent$institution.id},
                                              v.institution.objs)[[1]]
                    if(is.null(institution.obj) || institution.obj$institution.id < 1){
                        print("Institution object not found in the case of second agent being added!")
                        print(paste("Institution id sought:", firstagent$institution.id, " IDs present in v.institution.objs:",
                                    paste(sapply(v.institution.objs, function(instobj){instobj$institution.id}))))
                    }
                    else{
                        ## have second agen make a deposit
                        institution.obj$deposit(secondagent$sugarlevel - threshold, timeperiod, secondagent$agentid)
                        ## reduce second agent's sugarlevel
                        secondagent$sugarlevel <- secondagent$sugarlevel - threshold 
                    }
                } ## main's first else if
                else if(firstagent$institution.id == 0 && secondagent$institution.id > 0 &&
                        firstagent$sugarlevel > threshold){
                    ## assign agent 2's inst. obj. id to agent 1's institution.obj.id
                    firstagent$institution.id <- secondagent$institution.id
                    institution.obj <- Filter(function(instobj){instobj$institution.id ==
                                                                    secondagent$institution.id},
                                              v.institution.objs)[[1]]
                    if(is.null(institution.obj) || institution.obj$institution.id < 1){
                        print("Institution object not found in the case of first agent being added!")
                        print(paste("Institution id sought:", secondagent$institution.id, " IDs present in v.institution.objs:",
                                    paste(sapply(v.institution.objs, function(instobj){instobj$institution.id}))))
                    }
                    else{
                        ## have agent 1 make a deposit
                        institution.obj$deposit(firstagent$sugarlevel - threshold, timeperiod, firstagent$agentid)
                        ## reduce agent 1's sugarlevel
                        firstagent$sugarlevel <- firstagent$sugarlevel - threshold
                    }
                } ## main's second else if

                else{
                    ## both agents are associated with institutions
                    ## so do nothing
                } ## main's final else
                
            } ## for processing rows in distances.df
            
        }## if nrow of distances.df > 0
        v.institution.objs <- perform_alive_updates(lagents, v.institution.objs)
    }, #main tryCatch block
    error=function(condition){
        print("An error occurred")
        print(condition)
    },
    finally={
        return(v.institution.objs) ## return the vector of institution objects
    } ##finally
    ) ## end of tryCatch 
} ## end of perform_institution_actions


animate <- function(lagents, sugarscape, carcap, regenrt, periods, adensity,
                    birthrate, inboundrate, outboundrate, initial_sug_level,
                    metabol_range, vision_range, threshold, v.institution.objs){

    l.matrices.suglevels <- list()
    ds.agent.counts <- data.frame((1:periods), (1:periods))
    names(ds.agent.counts) <- c("alive", "dead")

    l.inst.matrices.suglevels <- list()
                                       
    for(tstep in (1:periods)){
        lagents <- perform_agent_motions(lagents, sugarscape, v.institution.objs,
                                         threshold, tstep)
      
        ## v.institution.objs <- perform_institution_actions(lagents, sugarscape,
        ##                                                   threshold, tstep,
        ##                                                   v.institution.objs)
        l.newagents <- add_new_agents(lagents, sugarscape, birthrate,
                                      inboundrate, initial_sug_level,
                                      metabol_range, vision_range)
        
        l.matrices.suglevels[[tstep]] <- collate_suglevels(lagents)
        if(is.null(v.institution.objs)){
            print(paste("looks like v.institution.objs is empty! see below"))
            print(v.institution.objs)
            l.inst.matrices.suglevels[[tstep]] <- rep(0, length(l.matrices.suglevels))
        } 
        else{
            l.inst.matrices.suglevels[[tstep]] <- collate_suglevels(v.institution.objs)
        }
        
        sugarscape <- generate_sugarscape(nrow(sugarscape), carcap)
    }

    results <- collate_results(nrow(sugarscape), carcap, regenrt, adensity,
                               initial_sug_level, metabol_range, vision_range,
                               birthrate, inboundrate, outboundrate, threshold,
                               l.matrices.suglevels, l.inst.matrices.suglevels)
    print(paste("No. of elements after running the current sim.:", length(results)))
    results
}

runsim <- function(rep.num){
    timesteps <- 30 ## is kept constant across all combinations - may need to increase
    source("agent-class-v8.R") ## contains the definition of the Agent class
    source("institution-v8.R")
        
    combos.df <- read.csv("parameter-ranges-testing.csv",
                          header = TRUE, sep=",")
    print(paste("There are total of", nrow(combos.df), "combinations to process."))
    tryCatch({
        dfrows <- foreach(rownum = 1:nrow(combos.df), .combine='rbind', 
                          .packages=c('ineq', 'assertthat', 'R6'),
                          .inorder=FALSE) %dopar% {
                              ## browser()
                              print(paste("Going to process combination:", rownum))
                              side <- combos.df$Side[rownum] ## number of cells on each side of the sugarscape
                              capacity <- combos.df$Capacity[rownum] ## sugar capacity in each cell
                              regenerationrate <- combos.df$RegRate[rownum] ## amt. of sugar growthin each time step, at each cell
                              agentdensity <- combos.df$Adensity[rownum] ## how densely the sugarscape if populated by aents
                              
                              popsize <- ceiling(agentdensity * side^2) ## compute the agent population
                              
                              metabol_range <- c(1, combos.df$MtblRate[rownum]) ## range of possible agent metabolic rates
                              vision_range <- c(1, combos.df$VsnRng[rownum]) ## range of how agent's vision in NSWE directions
                              initial_sug_level <- combos.df$InitSgLvl[rownum] ## initial sugarlevel of each agentsimulation
                              birthrate <- combos.df$Birthrate[rownum] ## birthrate of adding new agents to sugarscape
                              inboundrate <- combos.df$InbndRt[rownum] ## inbound migration rate
                              outboundrate <- combos.df$OtbndRt[rownum] ## outbound migration rate
                              threshold <- combos.df$Threshold[rownum] ## threshold for forming/joining new institution
                              sugarscape <- generate_sugarscape(side, capacity)
                              lst_agents <- generate_agents(popsize, initial_sug_level, metabol_range,
                                                            vision_range, side)
                              
                              v.institution.objs <- NULL
                              
                              results <- animate(lst_agents, sugarscape, capacity, regenerationrate,
                                                 timesteps, agentdensity, birthrate, inboundrate,
                                                 outboundrate, initial_sug_level, metabol_range,
                                                 vision_range, threshold, v.institution.objs)
                              print(paste("Finished processing combination:", rownum))
                              print("Results:")
                              print(results)
                              results
                              ## readline("Enter enter")
                          } ## end of foreach-dopar

        print(dfrows)
        
        dfrows <- data.frame(dfrows)
        print(paste("No. of columns in dfrows is:", ncol(dfrows)))
        names(dfrows) <- c("side", "carryingcapacity","regenrate", "agentdensity",
                           "initialsuglevel", "metabolrangemax", "visionrangemax",
                           "birthrate", "inboundrate", "outboundrate", "threshold",
                           c(unlist(sapply(1 : timesteps,
                                           function(tstep){
                                               paste("period_", tstep, sep="",
                                                     collapse="")}))),
                           c(unlist(sapply(1 : timesteps,
                                           function(tstep){
                                               paste("inst_period_", tstep, sep="",
                                                     collapse="")})))
                           )
        fname <- paste("replication-", rep.num, "-results.csv",sep="", collapse="")
        write.table(dfrows, file=fname, sep=";", row.names=FALSE)}, ## try part

        error = function(cond){
            print("Encountered an error")
            print(cond)
        },
        
        finally={
        print("In finally!!")}
        ) ## end try catch
}

## replicate the simulation several times, for each combination of input parameters
num.repetitions <- 20 #10
for(repetition in 1:num.repetitions){
    runsim(repetition)
    print(paste("Completed repetition:", repetition))
}

## collate the dataset
finaldf <- foreach(rep.num = 1:num.repetitions, .combine='rbind', .inorder=TRUE) %dopar%{
    fname <- paste("replication-", rep.num, "-results.csv",sep="", collapse="")
    read.csv(fname, header=T, sep=";")
}
## store the results to disk
write.table(finaldf, "collated-data.csv", sep=";", row.names=FALSE)
