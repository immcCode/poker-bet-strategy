library(caret)
library(Hmisc)
library(data.table)
library(tidyverse)
library(corrplot)
library(Boruta)

# Raw Data Directory
dataDir <- file.path('data', 'holdem')

# Source in other R files
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if (grepl('DATS6450', nm)) next
    if (trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if (trace) cat("\n")
  }
}
sourceDir('R')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in Files ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Will use the year of 2000 as our raw data for modeling
allDat <- pullDataForYear('2000')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter to only certain players ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# check to see the distribution of hands player by player and games won by player
playerStats <- 
  allDat %>% 
  select(timestamp, playername, won) %>% 
  mutate('winner' = ifelse(won > 0, 1, 0)) %>% 
  group_by(playername) %>% 
  summarize(handsPlayed = n(),
            gamesWon = sum(winner))

keepAmt <- .01

ggplot(playerStats, aes(handsPlayed)) +
  geom_histogram(binwidth = 1000) + 
  ggtitle('Distribution of Number of Hands Played by Player') +
  geom_vline(xintercept = (quantile(playerStats$handsPlayed, 1 - keepAmt)),
             col = 'red')

ggplot(playerStats, aes(gamesWon)) +
  geom_histogram(binwidth = 100) + 
  ggtitle('Distribution of Number of Hands Won by Player') + 
  geom_vline(xintercept = (quantile(playerStats$gamesWon, 1 - keepAmt)),
             col = 'red')

cat('Number of players:', nrow(playerStats), 
    '\nWill minimize to', keepAmt, '% of players:', 
    floor(nrow(playerStats) * keepAmt), 'players.\n')

# Given that there are The majority of players won less than 100 hands. Will keep only
# those with over 100 wins to make sure there is a good sample of 
# wins amongst the remaining players
keepPlayers <- 
  playerStats %>% 
  filter(handsPlayed > quantile(playerStats$handsPlayed, 1 - keepAmt)) %>% 
  pull(playername)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get player actions per hand ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# In this section, will get thes sequence of player actions per hand between 
# whom bet, checked, folded, etc

# Create a vector to substitute actions
actionNames <- c('-', 'B', 'f', 'k', 'b', 'c', 'r', 'A', 'Q', 'K')
actions <- c('no_action', 'blind', 'fold', 'check', 'bet', 'call',
                 'raise', 'all_in', 'quit', 'kicked')
names(actions) <- actionNames

rounds <- c('preflop', 'flop', 'turn', 'river')
rndNums <- length(rounds)

# Get list of actions per betting round, combine into 1 table, update
# order variable to account for differences due to per round and sort
allActs <- 
  lapply(rounds, createActionPerRound) %>% 
  bind_rows %>%
  mutate(ord = case_when(
            round == 'flop' ~ ord + 100,
            round == 'turn' ~ ord + 200,
            round == 'river' ~ ord + 300,
            TRUE             ~ ord
            )) %>% 
  arrange(timestamp, ord) %>% 
  as.tbl

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Find walks, aka when everyone folds to big blind ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
potenAllFoldHands <-
  allDat %>%
  filter(playersflop == 0) %>% 
  pull(timestamp) %>% 
  unique

allFoldHands <-
  allActs %>%
  filter(timestamp %in% potenAllFoldHands,
         round == 'preflop') %>% 
  group_by(timestamp, players) %>% 
  summarise(n = n()) %>% 
  filter(players + 1 == n) %>% 
  pull(timestamp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter data, remove walks and unnecessary players ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

allActs <- 
  allActs %>% 
  filter(!(timestamp %in% allFoldHands) & playername %in% keepPlayers)

allDat <-
  allDat %>% 
  filter(!(timestamp %in% allFoldHands) & playername %in% keepPlayers)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create extra columns that can be used for analysis ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Whom bet first during the round ----
initBet <- 
  allActs %>% 
  filter(action == 'bet') %>% 
  mutate('flopInitBet' = ifelse(round == 'flop', 1, 0),
         'turnInitBet' = ifelse(round == 'turn', 1, 0),
         'riverInitBet' = ifelse(round == 'river', 1, 0)) %>% 
  select(-action, -round, -players, -ord)

# Average Num time raises the pot and voluntarily puts chips into the pot ----
allPlayerActs <-
  allDat %>% 
  select(timestamp, playername) %>% 
  mutate('timeplayer' = paste(timestamp, playername, sep = "|")) %>%
  pull(timeplayer) %>% 
  expand.grid('round' = rounds, stringsAsFactors = FALSE) %>% 
  tidyr::separate(Var1, c('timestamp', 'playername'), sep = '\\|') %>% 
  mutate(timestamp = as.numeric(timestamp))
  
betStats <- 
  allPlayerActs %>% 
  left_join(allActs, by = c('timestamp', 'playername', 'round')) %>% 
  
  mutate(pip = ifelse(action %in% c('bet', 'raise', 'call'), 1, 0),
         raise = ifelse(action %in% c('bet', 'raise'), 1, 0)) %>% 
  group_by(timestamp, playername, round) %>% 
  summarize(numBets = sum(raise),
            numPip = sum(pip)) %>% 
  group_by(timestamp, playername) %>% 
  summarize(avgRaises = mean(numBets),
            avgPip = mean(numPip))

# Num time check raises ----

### Create obj for later processing
checkRaiseDat <- 
  allActs %>% 
  filter(action %in% c('check', 'raise'))

### Find players with two actions in a round
twoActs <- 
  checkRaiseDat %>% 
  group_by(timestamp, round, playername) %>% 
  summarize(n = n()) %>% 
  filter(n > 1) %>% 
  select(-n)

finalCheckRaises <- 
  checkRaiseDat %>% 
  inner_join(twoActs, by = c('timestamp', 'playername', 'round')) %>% 
  group_by(timestamp, playername, round) %>% 
  mutate(prevAct = dplyr::lag(action),
         chkRaise = ifelse(action == 'raise' & prevAct == 'check', 1, 0)) %>% 
  group_by(timestamp, playername) %>% 
  summarize(num_checkRaises = sum(chkRaise))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create final dataset for analysis ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Will join over calculated columns, convert mising data to 0, create 
# dependent variable winner, and drop several variables:
#   - the cards since looking for betting strategy independent of cards
#   - action variables since accounted for in calculated variables
#   - extra player count variables --- all of the variables will be highly 
#     correlated so will only keep percentage of players at the flop
#   - money variables, since strategy should be independent of money
#   - extra unnecessary variables 

cleanedDat <-
  allDat %>%
  select(timestamp, playername, position, won, playersflop, players) %>% 
  mutate('winner' = ifelse(won > 0, 1, 0)) %>% 
  left_join(initBet, by = c('timestamp', 'playername')) %>% 
  left_join(betStats, by = c('timestamp', 'playername')) %>% 
  left_join(finalCheckRaises, by = c('timestamp', 'playername')) %>%
  mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
  select(-won, -playersflop) %>% # remove variables
  as.tbl

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exploratory Analysis ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create factor var object for later
factorVars <- c('flopInitBet', 'turnInitBet', 'riverInitBet', 'winner')
charVars <- c('playername')

filteredDat <- 
  cleanedDat %>% 
  filter(playername %in% keepPlayers) %>% 
  select(-one_of(charVars), -timestamp, -players, -avgPip)

# Look at correlation plot for numeric variables. Will leave binary variables in
filteredDat %>% 
  as.matrix %>% 
  cor %>% 
  corrplot(method = 'number')

# Avg Raises is highly correlated with winner (makes sense)

#~~~~~~~~~~~~~~
# Check what if drop raises instead of pip
filteredDat <- 
  cleanedDat %>% 
  filter(playername %in% keepPlayers) %>% 
  select(-one_of(charVars), -timestamp, -players, -avgRaises)

# Look at correlation plot for numeric variables. Will leave binary variables in
filteredDat %>% 
  as.matrix %>% 
  cor %>% 
  corrplot(method = 'number')

# And check if include both
filteredDat <- 
  cleanedDat %>% 
  filter(playername %in% keepPlayers) %>% 
  select(-one_of(charVars), -timestamp, -players)

# Look at correlation plot for numeric variables. Will leave binary variables in
filteredDat %>% 
  as.matrix %>% 
  cor %>% 
  corrplot(method = 'number', type = 'lower')




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Feature Selection ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Check normal logistic regression
checkLogRegForVars <- function(dat) {
  # Original model with all variables
  logRegMod <- glm(winner ~ ., family = 'binomial',
                   data = dat)
  # select sig. variables
  sumLogReg <- summary(logRegMod)
  print(sumLogReg)
  
  toselect.x <- sumLogReg$coeff[-1,4] < 0.01
  relevant.x <- names(toselect.x)[toselect.x == TRUE] 
  
  # formula with only sig variables
  sig.formula <- as.formula(paste("winner ~",paste(relevant.x, collapse = "+")))
  sig.model <- glm(formula = sig.formula, family = 'binomial',
                   data = dat)
  
  print(summary(sig.model))
  return(relevant.x)
}

fixedFactorDat <- 
  filteredDat %>% 
  mutate_at(vars(one_of(factorVars)), as.factor)

sigVars1 <- checkLogRegForVars(filteredDat)

# Feature selection through Boruta
borutaOut <- Boruta(winner ~ .,
                    data = fixedFactorDat,
                    doTrace = 2)

plot(borutaOut)

# collect Confirmed and Tentative variables
borutaSignif <- 
  borutaOut$finalDecision[borutaOut$finalDecision %in% 
                                c("Confirmed", "Tentative")] %>% 
  names

print(borutaSignif) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add heiarchical columns ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# breaks <- c(0, .14, .23, .32, .4, 1)
# vpipStats <- c('very-tight', 'tight', 'semi-loose', 'loose', 'very-loose')
vpipStats <- c('tight', 'semi-loose', 'loose')

# Calculate VPIP percentage, number of hands voluntarily put chips in pot
realVPIP <- 
  allActs %>% 
  filter(round == 'preflop') %>% 
  mutate(pip = ifelse(action %in% c('bet', 'raise', 'call'), 1, 0),
         raise = ifelse(action %in% c('bet', 'raise'), 1, 0)) %>% 
  group_by(playername, timestamp) %>% 
  summarize(pip = sum(pip),
            raise = sum(raise)) %>% 
  mutate(handpip = ifelse(pip > 0, 1, 0),
         handraise = ifelse(raise > 0, 1, 0)) %>% 
  filter(!(timestamp %in% allFoldHands)) %>% 
  group_by(playername) %>% 
  summarize(allpip = sum(handpip),
            allraise = sum(handraise),
            n = n()) %>% 
  mutate(vpip = allpip / n,
         pfr = allraise / n,
         pfrVpip = pfr / vpip,
         vpipCut = Hmisc::cut2(vpip, g = 3),
         vpipFinal = factor(vpipCut, levels(vpipCut), vpipStats))

ggplot(realVPIP, aes(pfrVpip)) + 
  geom_histogram(binwidth = .01) +
  ggtitle('Distribution of PFR/VPIP')

ggplot(realVPIP, aes(pfr)) + 
  geom_histogram(binwidth = .01) +
  ggtitle('Distribution of PFR')

ggplot(realVPIP, aes(vpip)) + 
  geom_histogram(binwidth = .01) +
  ggtitle('Distribution of VPIP')

finalRealVPIP <- 
  realVPIP %>% 
  select(playername, vpipFinal) %>% 
  rename('vpipType' = vpipFinal)

# Create experience column (this may need to be run separately to save memory)
rawDataYears <- getListOfYearsFromRawFiles(dataDir)
rawDataYears <- rawDataYears[as.numeric(rawDataYears) < 2000]

prior_exp <- 
  lapply(rawDataYears, getExperienceByYear, keepPlayers) %>% bind_rows

priorExp <-
  prior_exp %>% 
  group_by(playername) %>% 
  summarize(years = n(),
            handsPlayed = sum(n),
            totalWins = sum(wins)) %>% 
  arrange(-handsPlayed)

# Need to drop the top two players since they are so much bigger
ggplot(filter(priorExp, !(playername %in% c('r00lbot', 'kfish'))),
              aes(handsPlayed)) + 
  geom_histogram(binwidth = 1000) +
  ggtitle('Distribution of Total Hands Played by Players Prior to 2000')

# Need to drop the top two players since they are so much bigger
ggplot(priorExp, aes(totalWins)) + 
  geom_histogram(binwidth = 1000) +
  ggtitle('Distribution of Total Hands Won by Players Prior to 2000')

# Will use following cutoffs: 
#   7,500 as expert cutoff
#       0 as mediocre cutoff since graph peaks at 2,000
# no data as beginner

finalPriorExp <-
  tibble(playername = keepPlayers) %>% 
  left_join(priorExp, by = 'playername') %>% 
  mutate(priorExp = case_when(
          handsPlayed > 7500   ~ 'expert',
          handsPlayed >    0   ~ 'mediocre',
          TRUE                 ~ 'beginner')) %>% 
  select(playername, priorExp)

table(finalPriorExp$priorExp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create final dataset ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Based on correlation plot, will not include average Raises
# all predictors are okay
# possibly remove playername too?
finalDat <-
  cleanedDat %>% 
  select(-timestamp, -players) %>% 
  left_join(finalRealVPIP, by = 'playername') %>% 
  left_join(finalPriorExp, by = 'playername')

data.table::fwrite(finalDat, file.path('data', 'finalOutput.csv'))

# Uncomment to read in data instead of re-running all previous data
finalDat <-
  data.table::fread(file.path('data', 'finalOutput.csv'))
#factorVars <- c('flopInitBet', 'turnInitBet', 'riverInitBet', 'winner',
#                  'vpipType', 'priorExp')
factorVars <- c('winner','vpipType', 'priorExp')
charVars <- c('playername')


finalDat <-
  finalDat %>% 
  select(-playername) %>% 
  mutate_at(vars(one_of(factorVars)), as.factor)

# Update levels of heiararchy variables
exp <- c('expert', 'mediocre', 'beginner')
vpipStats <- c('tight', 'semi-loose', 'loose')

finalDat$vpipType <- factor(finalDat$vpipType, levels = vpipStats)
finalDat$priorExp <- factor(finalDat$priorExp, levels = exp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Split into train and test datasets ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setting seed for reproducibility
# set.seed(1234)

# trIndx <-
#   caret::createDataPartition(finalDat$winner, p = .8, list = FALSE, times = 1)
# 
# trainData <- finalDat[trIndx, ]
# testData <- finalDat[-trIndx, ]
# 
# table(trainData$winner)
# table(testData$winner)

# Account for massive class disparity in the training data by down-weighting
# since we have ample number of observations in train dataset
# set.seed(3456)
# downTrainData <- 
#   caret::downSample(x = select(trainData, -winner),
#                     y = trainData$winner, yname = 'winner') %>% 
#   as.tbl

test <- TRUE
test <- FALSE

if (test) {
  table(finalDat$priorExp, finalDat$vpipType, finalDat$winner)
  
  lngth <- 
    length(unique(finalDat$vpipType)) * length(unique(finalDat$priorExp)) * 2
  testDatList <- vector('list', lngth)
  
  set.seed(3456)
  ctr <- 1
  for (i in levels(finalDat$vpipType)) {
    for (n in levels(finalDat$priorExp)) {
      for (w in levels(finalDat$winner)) {
        tmp <- 
          finalDat %>% 
          filter(vpipType == i & priorExp == n & winner == w)
        
        indx <- sample(1:nrow(tmp), 20)
        
        testDatList[[ctr]] <- tmp[indx, ]
        ctr <- ctr + 1
      }
    } 
  }
  
  modelDat <- bind_rows(testDatList)
  
  table(modelDat$priorExp, modelDat$vpipType, modelDat$winner)

} else {
  set.seed(3456)
  modelDat <- 
    caret::downSample(x = select(finalDat, -winner),
                      y = finalDat$winner, yname = 'winner') %>% 
    as.tbl
  
  # Much better
  table(modelDat$winner)
}

# Make sure y variable is numeric, not factor for model to work
if (class(modelDat$winner) == 'factor') {
  modelDat$winner <- as.numeric(as.character(modelDat$winner))
}

# Make sure y variable is between 0 and 1
stopifnot(modelDat$winner >= 0 & modelDat$winner <= 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Build Model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

yName = "winner"
xName = c("avgRaises", 'avgPip', 'num_checkRaises',
          'flopInitBet', 'turnInitBet', 'riverInitBet')
sName = c('priorExp', 'vpipType')
fileNameRoot = "Project-JAGS-Model-TEST1-" 
numSavedSteps=15000 ; 
thinSteps=2
#.............................................................................
graphFileType = "eps" 
#------------------------------------------------------------------------------- 

# Generate the MCMC chain:
startTime = proc.time()
mcmcCoda = genMCMC2Heirarchy( data = as.data.frame(modelDat) , 
                    xName = xName, yName = yName, sName = sName,
                    numSavedSteps = numSavedSteps , thinSteps = thinSteps , 
                    saveName = fileNameRoot )
endTime = proc.time()

show(endTime - startTime)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=NULL , saveType=NULL )
            # saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC(mcmcCoda , data = as.data.frame(modelDat), 
                   xName = xName, yName = yName, sName = sName, 
                   compVal = NULL, pairsPlot = TRUE, showCurve = FALSE,
                   saveName = fileNameRoot, saveType = graphFileType)
#------------------------------------------------------------------------------- 

sName = c('vpipType')
startTime = proc.time()
mcmcCoda = genMCMC1Heirarchy( data = as.data.frame(modelDat) , 
                              xName = xName, yName = yName, sName = sName,
                              numSavedSteps = numSavedSteps, thinSteps = thinSteps, 
                              saveName = fileNameRoot )
endTime = proc.time()
show(endTime - startTime)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=NULL , saveType=NULL )
  # saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC(mcmcCoda , data = as.data.frame(modelDat), 
                   xName = xName, yName = yName, sName = sName, 
                   compVal = NULL, pairsPlot = TRUE, showCurve = FALSE,
                   saveName = fileNameRoot, saveType = graphFileType)
#------------------------------------------------------------------------------- 


