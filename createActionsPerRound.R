# Create steps in the game. To do that will cycle through each 
# betting timeframe (preflop, flop, etc..) and split out actions into multiple
# columns. After deleting out no-action and stacking columns underneath each
# other, it should be in order.
createActionPerRound <- function(round) {
  
  # split up actions during the round
  roundVar <- paste0(round, 'action')
  
  splitActs <- 
    stringr::str_split_fixed(allDat[[roundVar]], pattern = "", n = Inf)
  colnames(splitActs) <- paste0('round', seq(1, ncol(splitActs)))
  playActDt <- data.frame('timestamp' = allDat$timestamp,
                          'players' = allDat$players,
                          'playername' = allDat$playername, 
                          'position' = allDat$position,
                          splitActs, 
                          stringsAsFactors = FALSE)
  # order the actions
  byAct <- 
    playActDt %>% 
    melt(id.vars = c('timestamp', 'players', 'playername', 'position'), 
         value.name = 'action') %>% 
    filter(!(action %in% c("", "-"))) %>% 
    mutate(ord = (as.numeric(substr(variable, 6, 6)) * players) 
           + position) %>% 
    select(-variable) %>% 
    mutate('round' = round,
           'action' = actions[action]) %>%
    mutate('action' = case_when(
      position == 1 & action == 'blind' ~ 'small_blind',
      position == 2 & action == 'blind' ~ 'big_blind',
      TRUE                              ~ action
    )) %>% 
    select(-position) %>% 
    unique
  
  return(byAct)
}
