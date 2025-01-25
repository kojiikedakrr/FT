loadDTSMdata <- function(flag.jpsmacro=FALSE) {
    ## Globals: creates Y, M.o, dates, mats, n.per, W
    n.per <<- 12

    start.date <- 19850101; start.month <- 198501
    end.date <- 20071231; end.month <- 200712

    ## yield data
    ## Anh Le's data
    load("data/le_data_monthly.RData") ## Y, dates
    mats <- c(3, 6, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120)
    Y <- Y[,mats]/n.per
    sel.sample <- dates>=start.date & dates<=end.date
    Y <- Y[sel.sample,]
    dates <- dates[sel.sample]
    cat("Start date:", min(dates), "\n")
    cat("End date:", max(dates), "\n")
    dates <<- dates
    mats <<- mats
    Y <<- Y

    ## yield PCs
    eig <- eigen(cov(Y))
    W <- t(eig$vectors)  ## rows contain loadings, first row is first PC
    W[1,] <- W[1,]/sum(W[1,])
    W <- W*1200  # yield factors are in percent per year
    W <<- W

    ## macro data
    if (flag.jpsmacro) {
        macro <- getJPSmacro()
        macro$GRO <- macro$GRO*100
        macro$INF <- macro$INF*100
        macro <- subset(macro, yyyymm<=end.month & yyyymm>=start.month)
    } else {
        macro <- read.csv("data/macro_monthly_corecpi_ugap.csv")
        macro <- subset(macro, DATE<=end.date & DATE>=floor(start.date/100)*100)
    }
    if (nrow(macro)!=nrow(Y))
        stop("yield and macro data do not have the same length")
    M.o <<- as.matrix(macro[,2:3])  # construct macro factors
}

getJPSmacro <- function()
    setNames(read.csv("data/jpsData_macro.csv", header=FALSE), c("yyyymm", "GRO", "INF"))

calcReturns <- function(yields, mats, h=12) {
    ## construct annual excess returns
    ## - expects columns y1, y2, etc in data.frame "yields"
    T <- nrow(yields)
    for (n in mats[mats>h]) {
        y.n <- yields[[paste('y', as.character(n/12), sep="")]]
        y.nmh <- yields[[paste('y', as.character(n/12-1), sep="")]]
        xr <- c(-(n-h)*y.nmh[(h+1):T] + n*y.n[1:(T-h)] - h*yields$y1[1:(T-h)], rep(NA, h))
        yields[paste('xr', as.character(n/12), sep="")] <- xr
    }
    yields$xr.avg <- rowMeans(as.matrix(yields[,grep("xr", names(yields))]))
    return(yields)
}

addPCs <- function(data, yield.cols, N=5) {
    ## calculate PCs of yields
    if (missing(yield.cols))
        stop("Need to provide yield.cols")
    Y <- data.matrix(data[,yield.cols])
    W <- eigen(cov(Y))$vectors[,1:N]
    PC <- data.frame(Y %*% W)
    names(PC) <- sapply(1:N, function(i) paste("PC", i, sep=""))
    return(data.frame(data, PC))
}

constructDataSet <- function(start = 198501, end = 200712, source.yields = "Le") {
    cat("Start date:", start, " End date:", end, "\n")

    macro <- getMacroData()

    jpsmacro <- getJPSmacro() ## this only available up to 2007-12
    macro <- merge(macro, jpsmacro) #, by="yyyymm", all=TRUE)

    altmeas <- getAltMeasSlackInfl()
    macro <- merge(macro, altmeas)

    ## yields
    if (source.yields == "Le") {
        yields <- getLeYields()
    } else if (source.yields == "GSW") {
        yields <- getGSWyields()
    } else if (source.yields == "JPS") {
        yields <- getJPSyields()
    }
    yields <- calcReturns(yields, mats)

    ## merge
    data <- merge(yields, macro)
    data <- data[do.call(order, data),]

    ## select subsample
    data <- subset(data, yyyymm>=start & yyyymm<=end)

    ## add principal components
    yield.cols <- min(which(names(data)=="y1"),
                      which(names(data)=="y0.5"),
                      which(names(data)=="y0.25")):which(names(data)=="y10")
    data <- addPCs(data, yield.cols)

    return(data)
}

getMacroData <- function() {
    ## read macro data
    ## read in macro variables, fix dates
    macro <- read.csv("data/data_macro.csv", na.strings="#N/A")
    macro <- subset(macro, select=-date_first)
    ## macro$date_first<- as.Date(macro$date_first, format="%m/%d/%Y")
    macro$yyyymm <- as.numeric(macro$yyyymm)

    ## fill in missing monthly nairu, calculate ugap
    require(zoo)
    macro <- transform(macro, nairu = na.locf(nairu))
    macro$ugap <- macro$ur - macro$nairu

    return(macro)
}

getAltMeasSlackInfl <- function() {
    altmeas <- read.csv("data/alt_measures_slack_inflation_data.csv", na.strings=c("#N/A", "#NA"))
    colnames(altmeas)[1] <- "date_first"
    altmeas$date_first<- as.Date(altmeas$date_first, format="%m/%d/%Y")
    altmeas$yyyymm <- as.numeric(format(altmeas$date_first, "%Y%m"))
    return(altmeas)
}

convertDate <- function(date)
    ## convert numeric date YYYYMMDD to Date object
    return(as.Date(as.character(date), format="%Y%m%d"))

getLeYields <- function() {
    load("data/le_data_monthly.RData")
    mats <<- c(3, 6, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120)
    yields <- data.frame(dates, 100*Y[, mats])
    colnames(yields) <- c('dates', sapply(mats/12, function(n) paste('y', as.character(n), sep="")))
    yields$date <- convertDate(yields$dates)
    yields$yyyymm <- floor(yields$dates/100)
    return(yields)
}
