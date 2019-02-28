library(R6)

Institution <- R6Class("Institution",
                  public = list(
                      institution.id = 0,
                      sugarlevel = 0,
                      v.member.ids = integer(),
                      alive = FALSE,
                      transaction.log = NA,

                      initialize = function(inst.id, vct.suglevels, vct.agent.ids,
                                            timeperiod){
                          self$institution.id <- inst.id
                          self$transaction.log <- data.frame(agentid = NA,
                                                             suglevel = 0,
                                                             action = NA,
                                                             timeperiod = NA)
                          self$deposit(vct.suglevels[1], timeperiod,
                                       vct.agent.ids[1], TRUE)
                          self$deposit(vct.suglevels[2], timeperiod,
                                       vct.agent.ids[2], FALSE)
                          
                          self$v.member.ids <- c(self$v.member.ids,
                                                 vct.agent.ids)
                          self$alive = TRUE
                      }, ##initialize

                      deposit = function(sugaramt, timeperiod, agentid,
                                         isfirstdeposit = FALSE){
                          if(!isfirstdeposit){
                              self$transaction.log <- rbind(self$transaction.log,
                                                       c(agentid, sugaramt,
                                                         "deposit", timeperiod))
                          }
                          else{
                              self$transaction.log$agentid[1] <- agentid
                              self$transaction.log$suglevel[1] <- sugaramt
                              self$transaction.log$action[1] <- "deposit"
                              self$transaction.log$timeperiod[1] <- timeperiod
                          }
                          
                          self$sugarlevel <- self$sugarlevel + sugaramt
                          invisible(self)
                      }, ## deposit

                  withdraw = function(sugaramt, timeperiod, agentid){
                      self$transaction.log <- rbind(self$transaction.log,
                                                    c(agentid, sugaramt,
                                                      "withdrawal", timeperiod))
                      self$sugarlevel <- self$sugarlevel - sugaramt
                      if(self$sugarlevel == 0){
                          self$alive <- FALSE
                      }
                      invisible(self)
                  }, ## withdraw

                  withdrawal_possible = function(amt){
                      if(amt <= self$sugarlevel ){
                          return(TRUE)
                      }
                      return(FALSE)
                  }, ## withdrawal_possible
                  
                  get_ledger = function(){
                      return(self$transaction.log)
                  }, ## get_ledger

                  append_agent_id = function(agid){
                      self$v.member.ids <- c(self$v.member.ids, agid)
                  },

                  get_member_agent_ids = function(){
                      self$v.member.ids
                  }
                  ) ##end of public

)## end of Institution

