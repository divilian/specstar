library(R6)

Agent <- R6Class("Agent",
                 public = list(
                     sugarlevel = 0,
                     agentid = NA,
                     metabol = 0,
                     vision = 0,
                     location.x = 0,
                     location.y = 0,
                     alive = TRUE,
                     institution.id = 0,
                     
                     initialize = function(aid, suglev, mtblr, vis,
                                           loc.x, loc.y){
                         self$agentid <- aid
                         self$sugarlevel <- suglev 
                         self$metabol <- mtblr
                         self$vision <- vis
                         self$location.x <- loc.x
                         self$location.y <- loc.y
                     }, ##end initialize()

                     move_load_consume_deposit = function(sugarscape, l_takenspots,
                                                  threshold, timeperiod,
                                                  v.institution.objs){
                         ## first, move, load-up, and consume
                         tryCatch({
                             ## browser()
                             if(self$alive){
                                 nextlocation <- private$get_next_loc(sugarscape,
                                                                      l_takenspots)
                                 ## if a new viable location is found
                                 if(!(self$location.x == nextlocation[1]) && 
                                    !(self$location.y == nextlocation[2])){
                                     ## set new location
                                     self$location.x <- nextlocation[1]
                                     self$location.y <- nextlocation[2]
                                 } ## new location if
                                 ## collect and consume sugar
                                 self$sugarlevel <- self$sugarlevel +
                                     sugarscape[nextlocation[1], nextlocation[2]] -
                                     self$metabol

                                 ## reduce sugarscape's sugarlevel at that location
                                 sugarscape[nextlocation[1], nextlocation[2]] <- 0

                                 ## if need more sugar, try borrowing
                                 if(self$sugarlevel <= 0) {
                                     if(self$institution.id > 0){ ## @@
                                         institution.obj <- Filter(function(instobj){
                                             instobj$institution.id ==
                                                 self$institution.id},
                                             v.institution.objs)[[1]]
                                         
                                         if(is.null(institution.obj)){## | is.na(institution.obj)){
                                             print(paste("Failed to obtain institute object",
                                                         "with id:", self$institution.id))
                                         }
                                         else{
                                             ## borrow.amt == self$metabol
                                             if(institution.obj$withdrawal_possible(self$metabol)){
                                                 institution.obj$withdraw(self$metabol, timeperiod,
                                                                          self$agentid)
                                                 self$sugarlevel <- self$sugarlevel + self$metabol
                                                 self$alive <- TRUE
                                             }
                                         }
                                     } ## @@
                                     else{
                                         ## print(paste("Agent", self$agentid, "could not",
                                                    ## "borrow because it wasn't a member."))
                                         self$alive <- FALSE
                                     }
                                 } ## self$sugarlevel <=0
                                 else{ ## have excess sugar 
                                     deposit.amt <- self$sugarlevel -
                                         (self$metabol + threshold)
                                     if(deposit.amt > 0){
                                         ## make a deposit
                                         if(self$institution.id > 0){ ## @@
                                             institution.obj <- Filter(function(instobj){
                                                 instobj$institution.id ==
                                                     self$institution.id},
                                                 v.institution.objs)[[1]]
                                             
                                             if(is.null(institution.obj)){## | is.na(institution.obj)){
                                                 print(paste("Failed to obtain institute object",
                                                             "with id:", self$institution.id))
                                             }
                                             else{
                                                 institution.obj$deposit(deposit.amt, timeperiod,
                                                                         self$agentid)

                                                 self$sugarlevel <- self$sugarlevel - deposit.amt
                                             }
                                         }
                                     }
                                     else{
                                         ## print(paste("Agent", self$agentid,
                                         ##             "did not have enough to deposit"))
                                     }
                                 } ## excess sugar else
                             } ## alive
                         }, ##try expression
                         error= function(condition){
                             print("Error encountered in move_load_consume_deposit")
                             print(condition)
                         }, ## error
                         finally={

                             invisible(self)                         
                         } ##finally
                         ) ## tryCatch
                     } ## move_load_consume_deposit

                 ), ##end of public
                 
                 private = list(
                     get_next_loc = function(sugarscape, l_takenspots){
                         ds <- expand.grid(1:nrow(sugarscape), 1:nrow(sugarscape))
                         l_allspots <- lapply(1:nrow(ds),
                                              function(rownum){
                                                  c(ds[rownum, 1], ds[rownum,2])
                                              })
                         l_takenspots <- l_takenspots[lapply(l_takenspots,
                                                             length) > 0]
                         l_freespots <- setdiff(l_allspots, l_takenspots)
                         ds <- data.frame(matrix(unlist(l_freespots),
                                                 nrow=length(l_freespots),
                                                 byrow=TRUE))
                         
                         dist <- unlist(lapply(1:nrow(ds),
                                               function(rownum){
                                                   sqrt((ds[rownum,1] -
                                                         self$location.x)^2 +
                                                        (ds[rownum,2] -
                                                         self$location.y)^2)
                                               }
                                               ))
                         ds <- data.frame(ds, dist)
                         names(ds) <- c("xloc", "yloc", "dist")
                         
                         ds <- ds[ds$dist <= self$vision, 1:2]
                         if(nrow(ds) == 0){
                             return(c(self$location.x, self$location.y))
                         }
                         ds <- data.frame(ds,
                                          unlist(lapply(1:nrow(ds),
                                                        function(rowno){
                                                            sugarscape[ds[rowno,1],
                                                                       ds[rowno,2]]}
                                                        )))
                         names(ds) <- c("xloc","yloc","suglevel")
                         ds <- ds[ds$suglevel > self$sugarlevel &
                                   ds$suglevel == max(ds$suglevel),]
 
                         if(nrow(ds) == 0){
                              return(c(self$location.x, self$location.y))
                         } else{
                             return(c(ds[1, 1], ds[1,2]))
                         }
                         
                     } ## end of get_next_loc
                     
                 ) ## end of private
                 ) ## end definition of Agent
