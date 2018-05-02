# Wrapper function to read in each month
readInRawFiles <- function(dataDir) {
  months <- paste(dataDir, list.files(dataDir), sep = '/')
  return(lapply(months, readInMonthRawFiles) %>% bind_rows)
}

# Get a list of months that exist in the data directory
getListofMonths <- function(dataDir) {
  paste(dataDir, list.files(dataDir), sep = '/')
}

# Function that will read in each month
readInMonthRawFiles <- function(monthDir, debug = FALSE) {
  cat('Reading in month dir:', monthDir, '\n')

  # Hand Database ----
  hdbFile <- file.path(monthDir, 'hdb')
  hdbColNames <- c('timestamp', 'gameset', 'game', 'players',
                   'flopAct', 'turnAct', 'riverAct', 'showdownAct', 'flop1', 
                   'flop2', 'flop3' ,'turn', 'river')
  
  # Read in HDB table
  cat('Reading in HDB file\n')
  # hdbRaw <- data.table::fread(hdbFile, header = FALSE, fill = TRUE, quote = "",
                              # col.names = hdbColNames, na.strings = "")
  hdbRaw <- scan(hdbFile, what = as.list(rep("", 13)), fill = TRUE,
                 na.strings = " ") %>% as.data.table
  names(hdbRaw) <- hdbColNames
  
  # Splits up flop, turn, river, and showdown 
  cleanHDB <- function(hdbRaw) {
    
    hdbTmp <- data.table::copy(hdbRaw)
    newColNames <- c('playersflop', 'potflop', 'playersturn',
                     'potturn', 'playersriver', 'potriver', 
                     'playersshowdown', 'potshowdown')
    
    vars <- c('flop', 'turn', 'river', 'showdown')

    for (var in vars) {
      newCols <- grep(var, newColNames, ignore.case = TRUE, value = TRUE)
      hdbTmp <- tidyr::separate_(hdbTmp, paste0(var, 'Act'), newCols, sep = '/')
    }
    
    return(hdbTmp)
  }
  
  hdbClean <- 
    cleanHDB(hdbRaw) 
  
  charVars <- c('flopAct', 'turnAct', 'riverAct', 'showdownAct', 'flop1', 
                'flop2', 'flop3' ,'turn', 'river')
  numVars <- setdiff(names(hdbClean), charVars)
 
  hdb <- 
    hdbClean %>% 
    mutate_at(numVars, as.numeric)
  
  
  # Roster Database ----
  
  # rosterFile <- file.path(monthDir, 'hroster')
  # roster <- data.table::fread(rosterFile, header = FALSE, fill = TRUE)

  
  # Player Database ----
  
  # Gets a list of all player files since each player's data is in a separate file
  pdbDir <- file.path(monthDir, 'pdb')
  playerFiles <- list.files(pdbDir)
  
  # Read in player data
  readPlayerDatabase <- function(player, skipNull = FALSE) {
    infile <- file.path(pdbDir, player)
    scan(infile, what = as.list(rep("", 13)), sep = '', fill = TRUE, 
         na.strings = "", quiet = TRUE, skipNul = skipNull) %>% as.data.table
  }
  
  cat('Reading in player database\n')
  player_data <- 
    lapply(playerFiles, readPlayerDatabase) %>% 
    data.table::rbindlist(use.names = TRUE)
  
  # Update column names
  playerColNames <- c('playername', 'timestamp', 'players', 'position', 
                      'preflopaction', 'flopaction', 'turnaction', 
                      'riveraction', 'bankroll', 'totalaction', 'won', 
                      'player_hand1', 'player_hand2')
  
  names(player_data) <- playerColNames
  
  # Check player files for issues --- used when initially reading in raw files
  # to clean up the files
  if (debug) {
    debugPlayerFiles(playerFiles, playerColNames, player_data)  
  }
  

  finalPlayer <- 
    player_data %>% 
    select(-players) # Drop number of players (exists in hdb)

  # Fix numeric columns and make sure player name is character
  playerCharVars <- c('playername', 'preflopaction', 'flopaction', 'turnaction', 
                      'riveraction', 'player_hand1', 'player_hand2')
  playerNumVars <- setdiff(names(finalPlayer), playerCharVars)
  
  finalPlayer <-
    finalPlayer %>% 
    mutate_at(playerNumVars, as.numeric)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Merge Files ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  allDat <- merge(hdb, finalPlayer, by = 'timestamp', all.x = TRUE)
  setorder(allDat, timestamp, position)
  return(allDat)
}


#~~~~~~~~~~~~~~~~
# Check to see if there are errors in the raw files
#~~~~~~
debugPlayerFiles <- function(playerFiles, playerColNames, player_data) {

  cat('Reading in player database skipped \n')
  player_dataSkip <-
    lapply(playerFiles, readPlayerDatabase, skipNull = TRUE) %>% 
    data.table::rbindlist(use.names = TRUE)
  
  cat('Finished reading in player database skipped \n')
  
  names(player_dataSkip) <- playerColNames
  
  summarizeData <- function(dat) {
    dat %>% 
      group_by(playername) %>% 
      summarize(n = n())
  }
  
  # Find the errors and stop processing if there are any
  findIssFiles <-
    summarizeData(player_data) %>% 
    left_join(summarizeData(player_dataSkip), by = 'playername') %>% 
    mutate(diff = n.x - n.y) %>% 
    filter(diff != 0 | is.na(diff))
  
  if (nrow(findIssFiles) > 0 | (nrow(player_data) != nrow(player_dataSkip))) {
    browser()
  }
}
