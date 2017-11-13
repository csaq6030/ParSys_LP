inputFile <- "/home/lee/Documents/ParSys_LP/exercise05/monteCarlo/output.dat"
inputFile <- "/home/lee/Documents/ParSys_LP/exercise05/nQueens/output.dat"
con  <- file(inputFile, open = "r")

means <- c()
stDe <- c()
count = 0;
t <- c()
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  line = unlist(strsplit(oneLine, ";"))
  if(line[1]=="Time"){
    count = count + 1
    t <- c(t, as.numeric(line[2]))
    if(count == 10){
      means <- c(means, mean(t))
      stDe <- c(stDe, sd(t))
      count = 0
      t <- c()
    }
  }
} 

close(con)

