# ======================================
# Read the data
# ======================================

readFRRF <- function(fn, dir="", specs=NULL, attr=NULL, ...){
  
  if (dir == "") Fn <- fn else Fn <- paste(dir, fn, sep="/")
  lines <- readLines(Fn, warn = FALSE)
  close(file(Fn))
  
  act2  <- grep("Act2", lines) 
  
  if (length(act2))
    ff <- parseFRRF.act2(lines, fn, specs, attr, ...)
  else 
    stop ("Can only read act2 files")
  ff
}

parseFRRF.act2 <- function(lines, fn, specs=NULL, attr=NULL, ...){
  
  getFRRFdate <- function(lines){
    idate <- grep("File date", lines)
    itime <- grep("File time", lines)
    date  <-  paste(strsplit(lines[idate], ",")[[1]][2], 
                    strsplit(lines[itime], ",")[[1]][2])
    lines <<- lines[-(1:(itime+1))]
    return(date)
  }    
  getFRRFfit <- function(lines){
    ip <- grep("Alpha:", lines)-1
    il <- grep("SErP:", lines)
    fits  <-  read.csv(textConnection(lines[ip:il]), header=TRUE)
    lines <<- lines[-((ip-1):(il+1))]
    return(fits)
  }    
  
  date  <- getFRRFdate(lines)
  fit   <- getFRRFfit(lines)
  skip  <- grep("Saq", lines) - 1
  
  LL <- lines[1:(skip-1)]
  LL <- LL[-which(LL=="")]
  pars <- read.csv(textConnection(LL), header=FALSE)
  lines <- lines[-(1:skip)]
  
  Dat <- read.csv(textConnection(lines), header=TRUE)
  Dat[,1] <- fn
  colnames(Dat)[1] <- "file"  
  Dat$date    <- date
  Dat$file    <- fn
  cn <- colnames(Dat)
  oldcolnames <- c("X.Chl.", "rP", "rP.1", "F.", "Fm.", "Fq..Fm.")
  newcolnames <- c("Chl", "rP_measured", "rP_fitted", "Fo", "Fm", "Fv/Fm")

  ij <- match(oldcolnames, cn)
  
  cn[ij] <- newcolnames
  colnames(Dat) <- cn 
  
  attributes(Dat)$station <- fn
  attributes(Dat)$date <- date 
  attributes(Dat)$pars <- pars
  attributes(Dat)$fitted.pars <- fit 
  attributes(ff)$processing <- c(attributes(ff)$processing, paste("Created at", Sys.time()))
  attributes(ff)$standardized <- FALSE
  
  if (length(specs)){
    ELL <- as.data.frame(matrix(nrow=nrow(Dat), ncol=length(specs), data=specs, byrow=TRUE))
    names(ELL) <- names(specs)
    Dat <- cbind(Dat, ELL)
  }
  ell <- list(...)
  if (length(ell)){
    ELL <- as.data.frame(matrix(nrow=nrow(Dat), ncol=length(ell), data=ell, byrow=TRUE))
    names(ELL) <- names(ell)
    Dat <- cbind(Dat, ELL)
  }
  if (length(attr))
    attr(Dat, which=names(attr)) <- attr
  Dat
}

# ======================================
# Standardize the data
# ======================================

standardizeFRRF <- function(frrf, Fblanc=0, keep=NULL, 
                            col = c("station", "date", "E", "Chl", "Fo", "Fm", "RSigma", "QR", "ADC"),
                            Ka=11800){  # instrument-specific constant, units of /m
  # save the attributes
  ATT <- attributes(frrf)
  ATT <- ATT[!names(ATT) %in% c("names", "row.names", "class")]
  
  # remove failed data points
  keep <- unique(c(col, keep))
  ff         <- frrf[, keep]
  ff$Fblanc  <- Fblanc
  ii <- which (is.na(ff$Fo))
  if (length(ii)) ff <- ff[-ii, ]
  
  # correct for the blanc (fluorescence in the absence of Chl)
  ff$Fo      <- with(ff,  Fo - Fblanc) 
  ff$Fm      <- with(ff,  Fm - Fblanc)
  
  # calculate values using the absorption method
  ff$Fq      <- with(ff,  Fm - Fo    )
  ff$a_LHII  <- with(ff,  Fo * Fm / Fq * (Ka / 10^6))     # cross-sectional surface of PSII system
  ff$FqFm    <- with(ff,  Fq/ Fm  )
  ff$JVPII   <- with(ff,  Fq/ Fm* a_LHII[1] * E * 3600/1000) # volumetric e-flux, [mmol m-3 h-1]
  
  # save attributes
  attributes(ff) <- c(attributes(ff), ATT)
  attributes(ff)$Fblanc <- Fblanc
  attributes(ff)$processing <- c(attributes(ff)$processing, paste("Standardized with Fblanc = ", Fblanc, "at", Sys.time()))
  attributes(ff)$standardized <- TRUE
  attributes(ff)$check <- checkFRRF(ff)  
  return(ff)
}  

standardizeFRRFsigma <- function(frrf, Fblanc=0, keep=NULL, # NOT WORKING
                            col = c("station", "date", "E", "Chl", "Fo", "Fm", "RSigma", "QR", "ADC"),
                            Ka=11800){  # instrument-specific constant, units of /m
  # save the attributes
  ATT <- attributes(frrf)
  ATT <- ATT[!names(ATT) %in% c("names", "row.names", "class")]
  
  # remove failed data points
  keep <- unique(c(col, keep))
  ff         <- frrf #[, keep]
  ff$Fblanc  <- Fblanc
  ii <- which (is.na(ff$Fo))
  if (length(ii)) ff <- ff[-ii, ]
  
  # correct for the blanc (fluorescence in the absence of Chl)
  ff$RCII    <- with(ff,  Fo/Sigma/1e-18*Ka/1e6*1/6.02e23) 
  ff$npsII   <- with(ff,  RCII/(ff$Chl/893.5/1000) )
  ff$Fm      <- with(ff,  Fm - Fblanc)
  
  # calculate values using the absorption method
  ff$Fq      <- with(ff,  Fm - Fo    )
  ff$a_LHII  <- with(ff,  Fo * Fm / Fq * (Ka / 10^6))     # cross-sectional surface of PSII system
  ff$FqFm    <- with(ff,  Fq/ Fm  )
  fac  <- 6.02e23/1e18/893.5
  ff$JVPII   <- with(ff,  -npsII * Fq/ Fm* Sigma[1] *fac* E * 3600/1000) # volumetric e-flux, [mmol m-3 h-1]
  
  # save attributes
  attributes(ff) <- c(attributes(ff), ATT)
  attributes(ff)$Fblanc <- Fblanc
  attributes(ff)$processing <- c(attributes(ff)$processing, paste("Standardized with Fblanc = ", Fblanc, "at", Sys.time()))
  attributes(ff)$standardized <- TRUE
  attributes(ff)$check <- checkFRRF(ff)  
  return(ff)
}  

# ======================================
# Check the data
# ======================================

checkFRRF <- function(frrf, 
                      QA_ADC = c(30, 80), 
                      QA_RSigma = c(0.035, 0.07),
                      QA_QR = c(6, NA),
                      plotit=FALSE)
{
  check <- data.frame()
  
  out <- sum(c(frrf$ADC < QA_ADC[1], 
               frrf$ADC > QA_ADC[2]))
  
  check <- data.frame(name = "ADC", OK = out==0,
                      minValue = min(frrf$ADC), maxValue = max(frrf$ADC),
                      minQA= QA_ADC[1], maxQA=QA_ADC[2], exceed = out)
  
  
  out <- sum(c(frrf$RSigma < QA_RSigma[1], 
               frrf$RSigma > QA_RSigma[2]))
  check <- rbind(check, data.frame(name = "Rsigma", OK = out==0,
                                   minValue = min(frrf$RSigma), maxValue = max(frrf$RSigma),
                                   minQA= QA_RSigma[1], maxQA=QA_RSigma[2], exceed = out))
  
  out <- sum(frrf$QR < QA_QR[1])
  check <- rbind(check, data.frame(name = "QR", OK = out==0,
                                   minValue = min(frrf$QR), maxValue = max(frrf$QR),
                                   minQA= QA_QR[1], maxQA=NA, exceed = out))
  
  if (plotit){
    labels <- NULL	
    out <- check[1,"exceed"]
    ylim <- range(c(frrf$ADC, QA_ADC))
    plot(frrf$E, frrf$ADC, ylim = ylim)
    abline(h = QA_ADC, lty = "dashed", lwd = 2, col = "red")
    
    if (out > 0) labels <- c(labels, paste(out, "values out of range for ADC: [", 
                                           paste(range(frrf$ADC), collapse=","), 
                                           "], should be in [", paste(QA_ADC, collapse=","), "]"))
    
    out <- check[2,"exceed"]
    ylim <- range(c(frrf$RSigma, QA_RSigma))
    plot(frrf$E, frrf$RSigma, ylim = ylim)
    abline(h = QA_RSigma, lty = "dashed", lwd = 2, col = "red")
    if (out > 0) labels <- c(labels, paste(out, "values out of range for RSigma: [", 
                                           paste(range(frrf$RSigma), collapse=","), 
                                           "], should be in [", paste(QA_RSigma, collapse=","), "]"))
    
    out <- check[3,"exceed"]
    ylim <- range(c(frrf$QR, QA_QR), na.rm = TRUE)
    plot(frrf$E, frrf$QR, ylim = ylim)
    abline(h = QA_QR, lty = "dashed", lwd = 2, col = "red")
    if (out > 0) labels <- c(labels, paste(out, "values out of range for QR: [", 
                                           paste(range(frrf$QR), collapse=","), 
                                           "], should be in [", paste(QA_QR, collapse=","), "]"))
    
    # plot message on the plot
    plot.new()
    for (i in 1:length(labels)){
      text(0, 1-(i*0.25), labels[i], adj = 0, cex = 0.75)
    }
  }
  check
}

# ======================================
# Fit the data
# ======================================

fitFRRF <- function(frrf, PAR="E", response="JVPII", model="EP", 
                    PARmax = NULL, main=NULL, plotit=TRUE, station=NULL){
  
  ATT <- attributes(frrf)
  if (is.null(station)) station <- frrf$station[1]
  
  notused <- NULL
  FRRF <- frrf
  if (! is.null(PARmax)){
    ii <- which(frrf[,PAR] < PARmax)
    FRRF <- frrf[ii,]
    notused <- frrf[-ii,]
  }
  
  # use the Eilers-Peeters model to fit the data
  FIT            <- PEfit(model=model, PAR = FRRF[,PAR], 
                          response = FRRF[, response])
  FitEP          <- data.frame(station=station, model=model, as.list(FIT$par))
  FitEP$nrow     <- nrow(FRRF)
  FitEP$maxE     <- max(FRRF$E)
  FitEP$meanChl  <- mean(FRRF$Chl) 
  FitEP$rsq      <- summary(FIT)$sigma
  
  if (is.null(main)) main <- station
   if (plotit){
     ylim <- range(frrf[, response])
     xlim <- range(frrf[, PAR])
     with(FRRF, plot(E, FRRF[,response], main=main, ylim=ylim, xlim=xlim))
    lines(EP(FIT$par, PAR=seq(0, max(frrf$E), length.out=100)))
    if (! is.null(notused))
      with(notused, points(E, JVPII, col="red"))
   }
  FitEP
} 

