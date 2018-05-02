# This function will return the years corresponding to the raw data monthly files
getListOfYearsFromRawFiles <- function(dataDir) {
  months <- getListofMonths(dataDir)
  yearmon <- substring(months, 13)
  years <- substr(yearmon, 1, 4) %>% unique
  years  
}

# Get total number of hands played and won for a list of players in a year
getExperienceByYear <- function(yr, playersInFile) {
  pullDataForYear(yr) %>% 
    select(timestamp, playername, won) %>% 
    filter(playername %in% playersInFile) %>% 
    mutate('year' = yr) %>% 
    sumDatByYear
}


pullDataForYear <- function(yr) {
  months <- getListofMonths(dataDir)
  monthsInYear <- grep(yr, months, value = TRUE)
  yearDat <- lapply(monthsInYear, readInMonthRawFiles) %>% bind_rows
  yearDat
}

sumDatByYear <- function(yearDat) {
  yearDat %>% 
    mutate('winner' = ifelse(won > 0, 1, 0)) %>% 
    group_by(playername, year) %>% 
    summarize(n = n(),
              wins = sum(winner)) %>% 
    arrange(-wins)
}
