# MSL Tools

# EuclideanDistance ----
euclideanDistance <- function(x1,y1,x2,y2) {
  df <- data.frame(x1,y1,x2,y2) %>%
    dplyr::mutate(EucDist = sqrt((x2 - x1)^2 + (y2 - y1)^2)) %>%
    dplyr::select(EucDist)
  
  return(df)
}

# Formant Analysis ----
formantAnalysis <- function(time,formant) {
  df <- data.frame(time,formant)
  df <- df %>%
    dplyr::mutate(time_incre = NA,
                  f_change = NA,
                  formant = na_if(formant,0))
  if (all(is.na(df$formant))){
    onset <- NA
    offset <- NA
    extent <- NA
    slope <- NA
    df <- df %>%
      dplyr::mutate(f_active = NA,
                    (f_steady = NA))
    percentActive_F <- NA
    f_change <- NA
    slope_2020 <- NA
    steady_2020 <- NA
    steadyDuration <- NA
    active_F <- NA
  } else {
    firstLast <- df %>%
      dplyr::select(formant) %>%
      na.omit(formant)
    onset <- as.numeric(data.table::first(firstLast))
    offset <- as.numeric(data.table::last(firstLast))
    extent <- base::abs(offset-onset)#*1000
    slope <- extent/((data.table::last(df$time) - data.table::first(df$time)))
    df <- df %>%
      dplyr::mutate(f_change = c(0,base::abs(diff(df$formant))),
                    time_change = c(0,base::abs(diff(df$time))),
                    time_incre = cumsum(time_change)
                    )
    
    # This is the loop to determine active transitions
    c <- 1
    c_1 <- c + 1
    fl <- as.numeric(length(df$formant))
    d <- fl
    while(c <= fl){
      if (!is.na(df$f_change[c])){
        if (df$f_change[c] >= .02){
          df$f_active[c] <- TRUE
        } else {
          df$f_active[c] <- FALSE
        }
      } else {
        df$f_active[c] <- NA
      } # If there is an NA in the formants, this will make the active status NA
      c <- c + 1
    }
    df$f_steady <- !df$f_active
    
    # Active: This is the loop to correct for times where slopes are not active for only 1 time sequence ----
    c <- 2
    b <- c - 1
    a <- c + 1
    df_na <- df %>%
      dplyr::filter(!is.na(f_active))
    if (NROW(df_na) != 0) {
      while (c <= (NROW(df_na)-1)) {
        if (df_na$f_active[c] == FALSE && df_na$f_active[b] == TRUE && df_na$f_active[a] == TRUE){
          df_na$f_active[c] <- TRUE
        }
        c <- c + 1
        b <- c - 1
        a <- c + 1
      }
      df_na <- df_na %>%
        dplyr::select(time,f_active)
      df <- df %>%
        dplyr::select(!f_active) %>%
        base::merge(df_na,by = "time",all = TRUE)
    }
    
    # Active: Calculate slopes ----
    slopes <- rle(df$f_active)$lengths
    df_na <- df %>%
      dplyr::filter(!is.na(f_active))
    if (!all(df_na$f_active == FALSE)){
      TRUE_slopes <- slopes
      TRUE_slopes[!rle(df$f_active)$values] <- NA
      longestSlope <- base::max(as.numeric(TRUE_slopes), na.rm = TRUE)
    } else {
      longestSlope <- 0
    }
    if (longestSlope <= 1){
      slope_2020 <- NA
      active_F <- NA
    } else {
      inds <- as.numeric(first(which(TRUE_slopes == longestSlope)))
      beginning <- ifelse(inds == 1,1,sum(slopes[1:(first(inds)-1)]) + 1)
      if ((beginning-1) + slopes[inds] == as.numeric(NROW(df$f_active))){
        end <- as.numeric(NROW(df$f_active))
      } else {
        end <- as.numeric(NROW(df$f_active)) - sum(slopes[((inds[1])+1):as.numeric(NROW(slopes))])
      }
      active_F <- df[beginning:end,]
      timeMaxFormant <- active_F %>%
        dplyr::filter(formant == base::max(as.numeric(active_F$formant), na.rm = TRUE)) %>%
        dplyr::select(time) %>%
        data.table::first() %>%
        as.numeric
      
      active_F <- active_F %>%
        dplyr::filter(time <= timeMaxFormant) %>%
        dplyr::select(time, formant)
      
      timeMinFormant <- active_F %>%
        dplyr::filter(formant == base::min(as.numeric(active_F$formant), na.rm = TRUE)) %>%
        dplyr::select(time) %>%
        data.table::first() %>%
        as.numeric
      
      active_F <- active_F %>%
        dplyr::filter(time >= timeMinFormant)
      
      slope_2020 <- base::abs((dplyr::first(active_F$formant)) - (dplyr::last(active_F$formant))) / ((dplyr::last(active_F$time) - dplyr::first(active_F$time)))
    }
    percentActive_F <- (sum(df$f_active, na.rm = TRUE)/as.numeric(length(df$f_active)))*100
    
    # Steady: This is the loop to correct for times where slopes are not steady for only 1 time sequence ----
    c <- 2
    b <- c - 1
    a <- c + 1
    df_na <- df %>%
      dplyr::filter(!is.na(f_steady))
    if (!NROW(df_na) == 0) {
      while (c <= (NROW(df_na)-1)) {
        if (df_na$f_steady[c] == FALSE && df_na$f_steady[b] == TRUE && df_na$f_steady[a] == TRUE){
          df_na$f_steady[c] <- TRUE
        }
        c <- c + 1
        b <- c - 1
        a <- c + 1
      }
      df_na %>%
        dplyr::select(time,f_steady) ->
        df_na
      df %>%
        dplyr::select(-f_steady) %>%
        base::merge(df_na,by.x = "time",all = TRUE) ->
        df
    }
    
    # Steady: Calculate slopes ----
    slopes <- rle(df$f_steady)$lengths
    df_na <- df[!is.na(df$f_steady),]
    if (!all(df_na$f_steady == FALSE)){
      TRUE_slopes <- slopes
      TRUE_slopes[!rle(df$f_steady)$values] <- NA
      longestSlope <- base::max(as.numeric(TRUE_slopes), na.rm = TRUE)
    } else {
      longestSlope <- 0
    }
    if (longestSlope <= 1){
      steady_2020 <- NA
      steadyDuration <- 0
    } else {
      inds <- as.numeric(first(which(rle(df$f_steady)$lengths == longestSlope)))
      beginning <- ifelse(inds == 1,1,sum(slopes[1:(first(inds)-1)]) + 1)
      if ((beginning-1) + slopes[inds] == as.numeric(NROW(df$f_steady))){
        end <- as.numeric(NROW(df$f_steady))
      } else {
        end <- as.numeric(NROW(df$f_steady)) - sum(slopes[((inds[1])+1):as.numeric(NROW(slopes))])
      }
      steady_F <- df[beginning:end,]
      steady_2020 <- base::abs((dplyr::first(steady_F$formant)) - (dplyr::last(steady_F$formant))) / ((dplyr::last(steady_F$time) - dplyr::first(steady_F$time)))
      steadyDuration <- last(steady_F$time) - first(steady_F$time)
    }
    percentSteady_F <- (sum(df$f_steady, na.rm = TRUE)/as.numeric(length(df$f_steady)))*100
  }
  list <- list("onset" = onset,
               "offset" = offset,
               "extent" = extent,
               "slope" = slope,
               "f_active" = df$f_active,
               "f_steady" = df$f_steady,
               "percentActive" = percentActive_F,
               "f_change" = df$f_change,
               "f_slopeSegment" = active_F,
               "f_slope2020" = slope_2020,
               "f_steady2020" = steady_2020,
               "f_steadyDuration" = steadyDuration
  )
  return(list)
}

# 1D Movement Analysis ----
movementAnalysis_1D <- function(time, movement_1D) {
  if (sum(!is.na(movement_1D)) < (.50 * NROW(movement_1D))) {
    df <- data.frame(time,movement_1D)
    df %>%
      tibble::add_column(active = NA) %>%
      tibble::add_column(steady = NA) ->
      df
    active <- df$active
    steady <- df$steady
    onset <- NA
    offset <- NA
    distance <- NA
    displacement <- NA
    speed_avg <- NA
    speed_max <- NA
    activeSpeed_avg <- NA
    activeSpeed_max <- NA
    steadySpeed_avg <- NA
    steadySpeed_max <- NA
    pathDist <- NA
    steadyDuration <- NA
  } else {
  df <- data.frame(time,movement_1D)
  noBlank <- subset(df, !is.na(movement_1D))
  onset <- dplyr::first(noBlank$movement_1D)
  offset <- dplyr::last(noBlank$movement_1D)
  distance <- sum(base::abs(diff(noBlank$movement_1D,lag=1, differences = 1)))
  displacement <- base::abs(onset-offset)
  L <- as.numeric(NROW(movement_1D))
  df$timeDiff[2:L] <- diff(df$time,lag=1, differences = 1)
  df$pathDist[2:L] <- diff(df$movement_1D,lag=1, differences = 1)
  
  # Set up the active column
  c <- 2
  df %>%
    tibble::add_column(active = NA) ->
    df
  while (c <= NROW(df$pathDist)){
    if (!is.na(df$pathDist[c])){
      if (base::abs(df$pathDist[c]) <= .18){
        df$active[c] <- FALSE
      } else {
        df$active[c] <- TRUE
      }
      c <- c + 1
    } else {
      df$active[c] <- NA
      c <- c + 1
    }
  }
  df$steady <- !df$active
  
  # Active: This is the loop to correct for times where slopes are not active for only 1 time sequence ----
  c <- 2
  b <- c - 1
  a <- c + 1
  df %>%
    base::subset(!is.na(df$active)==TRUE) ->
    df_na
  if (!NROW(df_na) == 0) {
    while (c <= (NROW(df_na)-1)) {
      if (df_na$active[c] == FALSE && df_na$active[b] == TRUE && df_na$active[a] == TRUE){
        df_na$active[c] <- TRUE
      }
      c <- c + 1
      b <- c - 1
      a <- c + 1
    }
    df_na %>%
      dplyr::select(time,active) ->
      df_na
    df %>%
      dplyr::select(-active) %>%
      base::merge(df_na,by.x = "time",all = TRUE) ->
      df
  }
  
  # Active: Calculate slopes ----
  slopes <- rle(df$active)$lengths
  df_na <- df[!is.na(df$active),]
  if (!all(df_na$active == FALSE)){
    TRUE_slopes <- slopes
    TRUE_slopes[!rle(df$active)$values] <- NA
    longestSlope <- base::max(as.numeric(TRUE_slopes), na.rm = TRUE)
  } else {
    longestSlope <- 0
  }
  if (longestSlope <= 1){
    activeSpeed_avg <- NA
    activeSpeed_max <- NA
  } else {
    inds <- as.numeric(first(which(TRUE_slopes == longestSlope)))
    beginning <- ifelse(inds == 1,1,sum(slopes[1:(first(inds)-1)]) + 1)
    if ((beginning-1) + slopes[inds] == as.numeric(NROW(df$active))){
      end <- as.numeric(NROW(df$active))
    } else {
      end <- as.numeric(NROW(df$active)) - sum(slopes[((inds[1])+1):as.numeric(NROW(slopes))])
    }
    active <- df[beginning:end,]
    active$activeSpeed <- active$pathDist/active$timeDiff
    activeSpeed_avg <- base::mean(base::abs(active$activeSpeed), na.rm = TRUE)
    activeSpeed_max <- base::max(base::abs(active$activeSpeed), na.rm = TRUE)
  }
  
  # Steady: This is the loop to correct for times where slopes are not steady for only 1 time sequence ----
  c <- 2
  b <- c - 1
  a <- c + 1
  df %>%
    base::subset(!is.na(df$steady)==TRUE) ->
    df_na
  if (!NROW(df_na) == 0) {
    while (c <= (NROW(df_na)-1)) {
      if (df_na$steady[c] == FALSE && df_na$steady[b] == TRUE && df_na$steady[a] == TRUE){
        df_na$steady[c] <- TRUE
      }
      c <- c + 1
      b <- c - 1
      a <- c + 1
    }
    df_na %>%
      dplyr::select(time,steady) ->
      df_na
    df %>%
      dplyr::select(-steady) %>%
      base::merge(df_na,by.x = "time",all = TRUE) ->
      df
  }
  
  # Steady: Calculate slopes ----
  slopes <- rle(df$steady)$lengths
  df_na <- df[!is.na(df$steady),]
  if (!all(df_na$steady == FALSE)){
    TRUE_slopes <- slopes
    TRUE_slopes[!rle(df$steady)$values] <- NA
    longestSlope <- base::max(as.numeric(TRUE_slopes), na.rm = TRUE)
  } else {
    longestSlope <- 0
  }
  if (longestSlope <= 1){
    steadySpeed_avg <- NA
    steadySpeed_max <- NA
    steadyDuration <- 0
  } else {
    inds <- as.numeric(first(which(TRUE_slopes == longestSlope)))
    beginning <- ifelse(inds == 1,1,sum(slopes[1:(first(inds)-1)]) + 1)
    if ((beginning-1) + slopes[inds] == as.numeric(NROW(df$steady))){
      end <- as.numeric(NROW(df$steady))
    } else {
      end <- as.numeric(NROW(df$steady)) - sum(slopes[((inds[1])+1):as.numeric(NROW(slopes))])
    }
    steady <- df[beginning:end,]
    steady$steadySpeed <- steady$pathDist/steady$timeDiff
    steadySpeed_avg <- base::mean(base::abs(steady$steadySpeed), na.rm = TRUE)
    steadySpeed_max <- base::max(base::abs(steady$steadySpeed), na.rm = TRUE)
    steadyDuration <- last(steady$time*1000) - first(steady$time*1000)
  }
  
  
  pathDist <- df$pathDist
  active <- df$active
  steady <- df$steady
  speed <- df$pathDist / df$timeDiff
  speed_avg <- base::mean(base::abs(speed), na.rm = TRUE)
  speed_max <- base::max(base::abs(speed), na.rm = TRUE)
  }
  
  list <- list("onset" = onset,
               "offset" = offset,
               "distance" = distance,
               "displacement" = displacement,
               "speed_avg" = speed_avg,
               "speed_max" = speed_max,
               "activeSpeed_avg" = activeSpeed_avg,
               "activeSpeed_max" = activeSpeed_max,
               "steadySpeed_avg" = steadySpeed_avg,
               "steadySpeed_max" = steadySpeed_max,
               "pathDist" = pathDist,
               "active" = active,
               "steady" = steady,
               "steadyDuration" = steadyDuration)
  return(list)
}

# 2D Movement Analysis ----
movementAnalysis_2D <- function(time,movement_x,movement_y) {
  if (sum(!is.na(movement_x)) < (.50 * NROW(movement_x))||sum(!is.na(movement_y)) < (.50 * NROW(movement_y))){
    onset_x <- NA
    offset_x <- NA
    onset_y <- NA
    offset_y <- NA
    distance <- NA
    displacement <- NA
    speed_avg <- NA
    speed_max <- NA
  } else {
  xy <- data.frame(time,movement_x,movement_y)
  xy <- subset(xy,!is.na(movement_x) & !is.na(movement_y))
  onset_x <- dplyr::first(xy$movement_x)
  offset_x <- dplyr::last(xy$movement_x)
  onset_y <- dplyr::first(xy$movement_y)
  offset_y <- dplyr::last(xy$movement_y)
  
  distCount <- NROW(xy) #Counter for Distance Measure
  c <- 1 # counter for distance measure
  distance <- 0 #Base for net distance
  p1 <- 1
  p2 <- p1 + 1
  xy$pathDist <- NA
  xy$time_diff <- NA
  while (c < distCount) {
    new_distance <- raster::pointDistance(c(xy$movement_x[p1],xy$movement_y[p1]),c(xy$movement_x[p2],xy$movement_y[p2]), lonlat=FALSE)
    distance <- distance + new_distance
    time_diff <- xy$time[p2] - xy$time[p1]
    p1 <- p1 + 1
    p2 <- p1 + 1
    c <- c + 1
    xy$pathDist[c] <- new_distance
    xy$time_diff[c] <- time_diff
  }
  distance <- sum(xy$pathDist,na.rm = TRUE)
  
  displacement <- sqrt((offset_x - onset_x)^2 + 
                         (offset_y - onset_y)^2)
  xy$speed <- xy$pathDist / xy$time_diff
  speed_avg <- base::mean(base::abs(xy$speed), na.rm = TRUE)
  speed_max <- base::max(base::abs(xy$speed), na.rm = TRUE)
  }
  list <- list("onset_x" = onset_x,
               "offset_x" = offset_x,
               "onset_y" = onset_y,
               "offset_y" = offset_y,
               "distance" = distance,
               "displacement" = displacement,
               "speed_avg" = speed_avg,
               "speed_max" = speed_max)
  return(list)
}

# Convex Hull ----
library(sp)
cHull <- function(x,y) {
  x1 <- na.omit(x)
  y1 <- na.omit(y)
  hpts <- chull(x = x1, y = y1)
  hpts <- c(hpts, hpts[1])
  xy.coords <- cbind(x1, y1)
  chull.coords <- xy.coords[hpts,]
  if (NROW(chull.coords) < 4) {
    cHull <- NA
  } else {
    chull.poly <- sp::Polygon(chull.coords, hole=F)
    cHull <- chull.poly@area
  }
  return(cHull)
}

# Plot Convex Hull
Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
}  

# Hz to Bark ----
hz2bark <- function(f) {
  z <- 26.81/(1 + 1960/f) - 0.53
  if (z < 2) {
    z <- z + (0.15 * (2 - z))
  } else if (z > 20.1) {
    z <- z + (0.22 * (z - 20.1))
  }
  return(z)
}

# AAVS ----
aavs <- function(f1_bark, f2_bark) {
  sgv <- base::sqrt(base::det(stats::cov(as.matrix(f1_bark, f2_bark))))
  return(sgv)
}

# Temp Test ----
#lblFilePath <- "Data/PD/Sub11/Sub11_1/Sub11_1_005_sync.lbl"
#fbwFilePath <- "Data/PD/Sub11/Sub11_1/SUB11_1_005_SYNC.FBW"
#F0FilePath <- "Data/PD/Sub11/Sub11_1/SUB11_1_005_SYNC.F0"
#decoupledMovementFilePath <- "Data/PD/Sub11/Decoupled Movement/Sub11_1_005_Decoupled.csv"
#includeF0Data = TRUE
#Group <- "PD"
#Database.ID <- "Sub11"
#Study.ID <- "PD01"

# test <- KASA(Database.ID = "Sub11",
#             Study.ID = "PD01",
#             Group = "PD",
#             saveData = TRUE,
#             lblFilePath = lblFilePath,
#             fbwFilePath = fbwFilePath,
#             include_F0Data = TRUE,
#             F0FilePath = F0FilePath,
#             decoupledMovementFilePath = decoupledMovementFilePath,
#             include_Duration_s = TRUE,
#             include_F2.slope_2020.Hz.ms = TRUE,
#             AS.include_F2.Active = TRUE)

library(tidyverse)

# Kinematic and Acoustic Segment Analysis (KASA) ----
KASA <- function(Database.ID,
                 Study.ID,
                 Group,
                 saveData = FALSE,
                 lblFilePath,
                 fbwFilePath,
                 include_F0Data = FALSE,
                 F0FilePath,
                 decoupledMovementFilePath,
                 # Temporal Measures
                 include_Duration_s = FALSE,
                 include_Duration_ms = FALSE,
                 include_Onset_ms = FALSE,
                 include_Offset_ms = FALSE,
                 include_Onset_s = FALSE,
                 include_Offset_s = FALSE,
                 include_TempMid_ms.Aco = FALSE,
                 include_TempMid_ms.Kin = FALSE,
                 include_F1.tempMid = FALSE,
                 include_F2.tempMid = FALSE,
                 include_F3.tempMid = FALSE,
                 include_dTF.x_tempMid = FALSE,
                 include_dTF.y_tempMid = FALSE,
                 include_dTB.x_tempMid = FALSE,
                 include_dTB.y_tempMid = FALSE,
                 # Max Constriction Info
                 include_TimeMaxCons = FALSE,
                 include_NormTimeMaxCons = FALSE,
                 include_F1_maxConstriction = FALSE,
                 include_F2_maxConstriction = FALSE,
                 include_F3_maxConstriction = FALSE,
                 include_dTF.x_maxConstriction = FALSE,
                 include_dTF.y_maxConstriction = FALSE,
                 # Pitch Measures (F0)
                 include_F0.mean = FALSE,
                 include_F0.min = FALSE,
                 include_F0.max = FALSE,
                 include_F0.sd = FALSE,
                 # Intensity Measures (dB)
                 include_dB.mean = FALSE,
                 include_dB.min = FALSE,
                 include_dB.max = FALSE,
                 include_dB.sd = FALSE,
                 # F1 measures
                 include_F1.onset = FALSE,
                 include_F1.offset = FALSE,
                 include_F1.extent = FALSE,
                 include_F1.slope_Hz.ms = FALSE,
                 include_F1.slope_2020.Hz.ms = FALSE,
                 # F2 measures
                 include_F2.onset = FALSE,
                 include_F2.offset = FALSE,
                 include_F2.extent = FALSE,
                 include_F2.slope_Hz.ms = FALSE,
                 include_F2.slope_2020.Hz.ms = FALSE,
                 include_F2.steadyDuration = FALSE,
                 # F3 measures
                 include_F3.onset = FALSE,
                 include_F3.offset = FALSE,
                 include_F3.extent = FALSE,
                 include_F3.slope_Hz.ms = FALSE,
                 include_F3.slope_2020.Hz.ms = FALSE,
                 # Formant Measures
                 include_minF2_F3 = FALSE,
                 # TF measures
                   # TF - X
                   include_TF.x_Onset = FALSE,
                   include_TF.x_Offset = FALSE,
                   include_TF.x_Distance = FALSE,
                   include_TF.x_Displacement = FALSE,
                   include_TF.x_SpeedAvg = FALSE,
                   include_TF.x_SpeedMax = FALSE,
                   include_TF.x_steadyDuration = FALSE,
                   # TF - Y
                   include_TF.y_Onset = FALSE,
                   include_TF.y_Offset = FALSE,
                   include_TF.y_Distance = FALSE,
                   include_TF.y_Displacement = FALSE,
                   include_TF.y_SpeedAvg = FALSE,
                   include_TF.y_SpeedMax = FALSE,
                   include_TF.y_steadyDuration = FALSE,
                   # TF - XY
                   include_TF.xy_Distance = FALSE,
                   include_TF.xy_Displacement = FALSE,
                   include_TF.xy_SpeedAvg = FALSE,
                   include_TF.xy_SpeedMax = FALSE,
                   include_TF.cHull = FALSE,
                 # Decoupled TF measures
                   # dTF - X
                   include_dTF.x_Onset = FALSE,
                   include_dTF.x_Offset = FALSE,
                   include_dTF.x_Distance = FALSE,
                   include_dTF.x_Displacement = FALSE,
                   include_dTF.x_SpeedAvg = FALSE,
                   include_dTF.x_SpeedMax = FALSE,
                   include_dTF.x_steadyDuration = FALSE,
                   # dTF - Y
                   include_dTF.y_Onset = FALSE,
                   include_dTF.y_Offset = FALSE,
                   include_dTF.y_Distance = FALSE,
                   include_dTF.y_Displacement = FALSE,
                   include_dTF.y_SpeedAvg = FALSE,
                   include_dTF.y_SpeedMax = FALSE,
                   include_dTF.y_steadyDuration = FALSE,
                   # dTF - XY
                   include_dTF.xy_Distance = FALSE,
                   include_dTF.xy_Displacement = FALSE,
                   include_dTF.xy_SpeedAvg = FALSE,
                   include_dTF.xy_SpeedMax = FALSE,
                   include_dTF.cHull = FALSE,
                 # TB measures
                   # TB - X
                   include_TB.x_Onset = FALSE,
                   include_TB.x_Offset = FALSE,
                   include_TB.x_Distance = FALSE,
                   include_TB.x_Displacement = FALSE,
                   include_TB.x_SpeedAvg = FALSE,
                   include_TB.x_SpeedMax = FALSE,
                   include_TB.x_steadyDuration = FALSE,
                   # TB - Y
                   include_TB.y_Onset = FALSE,
                   include_TB.y_Offset = FALSE,
                   include_TB.y_Distance = FALSE,
                   include_TB.y_Displacement = FALSE,
                   include_TB.y_SpeedAvg = FALSE,
                   include_TB.y_SpeedMax = FALSE,
                   include_TB.y_steadyDuration = FALSE,
                   # TB - XY
                   include_TB.xy_Distance = FALSE,
                   include_TB.xy_Displacement = FALSE,
                   include_TB.xy_SpeedAvg = FALSE,
                   include_TB.xy_SpeedMax = FALSE,
                   include_TB.cHull = FALSE,
                   # Decoupled TB measures
                   # dTB - X
                   include_dTB.x_Onset = FALSE,
                   include_dTB.x_Offset = FALSE,
                   include_dTB.x_Distance = FALSE,
                   include_dTB.x_Displacement = FALSE,
                   include_dTB.x_SpeedAvg = FALSE,
                   include_dTB.x_SpeedMax = FALSE,
                   include_dTB.x_steadyDuration = FALSE,
                   # dTB - Y
                   include_dTB.y_Onset = FALSE,
                   include_dTB.y_Offset = FALSE,
                   include_dTB.y_Distance = FALSE,
                   include_dTB.y_Displacement = FALSE,
                   include_dTB.y_SpeedAvg = FALSE,
                   include_dTB.y_SpeedMax = FALSE,
                   include_dTB.y_steadyDuration = FALSE,
                   # dTB - XY
                   include_dTB.xy_Distance = FALSE,
                   include_dTB.xy_Displacement = FALSE,
                   include_dTB.xy_SpeedAvg = FALSE,
                   include_dTB.xy_SpeedMax = FALSE,
                   include_dTB.cHull = FALSE,
                 # UL measures
                   # UL - X
                   include_UL.x_Onset = FALSE,
                   include_UL.x_Offset = FALSE,
                   include_UL.x_Distance = FALSE,
                   include_UL.x_Displacement = FALSE,
                   include_UL.x_SpeedAvg = FALSE,
                   include_UL.x_SpeedMax = FALSE,
                   # UL - Y
                   include_UL.y_Onset = FALSE,
                   include_UL.y_Offset = FALSE,
                   include_UL.y_Distance = FALSE,
                   include_UL.y_Displacement = FALSE,
                   include_UL.y_SpeedAvg = FALSE,
                   include_UL.y_SpeedMax = FALSE,
                   # UL - XY
                   include_UL.xy_Distance = FALSE,
                   include_UL.xy_Displacement = FALSE,
                   include_UL.xy_SpeedAvg = FALSE,
                   include_UL.xy_SpeedMax = FALSE,
                   include_UL.cHull = FALSE,
                 # LL measures
                   # LL - X
                   include_LL.x_Onset = FALSE,
                   include_LL.x_Offset = FALSE,
                   include_LL.x_Distance = FALSE,
                   include_LL.x_Displacement = FALSE,
                   include_LL.x_SpeedAvg = FALSE,
                   include_LL.x_SpeedMax = FALSE,
                   # LL - Y
                   include_LL.y_Onset = FALSE,
                   include_LL.y_Offset = FALSE,
                   include_LL.y_Distance = FALSE,
                   include_LL.y_Displacement = FALSE,
                   include_LL.y_SpeedAvg = FALSE,
                   include_LL.y_SpeedMax = FALSE,
                   # LL - XY
                   include_LL.xy_Distance = FALSE,
                   include_LL.xy_Displacement = FALSE,
                   include_LL.xy_SpeedAvg = FALSE,
                   include_LL.xy_SpeedMax = FALSE,
                   include_LL.cHull = FALSE,
                   # Decoupled LL measures
                   # dLL - X
                   include_dLL.x_Onset = FALSE,
                   include_dLL.x_Offset = FALSE,
                   include_dLL.x_Distance = FALSE,
                   include_dLL.x_Displacement = FALSE,
                   include_dLL.x_SpeedAvg = FALSE,
                   include_dLL.x_SpeedMax = FALSE,
                   # dLL - Y
                   include_dLL.y_Onset = FALSE,
                   include_dLL.y_Offset = FALSE,
                   include_dLL.y_Distance = FALSE,
                   include_dLL.y_Displacement = FALSE,
                   include_dLL.y_SpeedAvg = FALSE,
                   include_dLL.y_SpeedMax = FALSE,
                   # dLL - XY
                   include_dLL.xy_Distance = FALSE,
                   include_dLL.xy_Displacement = FALSE,
                   include_dLL.xy_SpeedAvg = FALSE,
                   include_dLL.xy_SpeedMax = FALSE,
                   include_dLL.cHull = FALSE,
                   # Lip Aperture
                   include_lipAperture.min = FALSE,
                 # Jaw measures
                 # Jaw - X
                 include_Jaw.x_Onset = FALSE,
                 include_Jaw.x_Offset = FALSE,
                 include_Jaw.x_Distance = FALSE,
                 include_Jaw.x_Displacement = FALSE,
                 include_Jaw.x_SpeedAvg = FALSE,
                 include_Jaw.x_SpeedMax = FALSE,
                 # Jaw - Y
                 include_Jaw.y_Onset = FALSE,
                 include_Jaw.y_Offset = FALSE,
                 include_Jaw.y_Distance = FALSE,
                 include_Jaw.y_Displacement = FALSE,
                 include_Jaw.y_SpeedAvg = FALSE,
                 include_Jaw.y_SpeedMax = FALSE,
                 # Jaw - XY
                 include_Jaw.xy_Distance = FALSE,
                 include_Jaw.xy_Displacement = FALSE,
                 include_Jaw.xy_SpeedAvg = FALSE,
                 include_Jaw.xy_SpeedMax = FALSE,
                 include_Jaw.cHull = FALSE,
                 
                 # Acoustic ROM Measures
                 include_AAVS = FALSE,
                 
                 # Acoustic Segments
                 AS.include_NormalizeTime = FALSE,
                 AS.include_F1.Bark = FALSE,
                 AS.include_F1.Active = FALSE,
                 AS.include_F1.Steady = FALSE,
                 AS.include_F2.Bark = FALSE,
                 AS.include_F2.Active = FALSE,
                 AS.include_F2.Steady = FALSE,
                 AS.include_F3.Bark = FALSE,
                 AS.include_F3.Active = FALSE,
                 AS.include_F3.Steady = FALSE,
                 AS.include_F2_F3constriction = FALSE,
                 
                 # Movement Segments
                 MS.include_TF_x.PathDist = FALSE,
                 MS.include_TF_x.Active = FALSE,
                 MS.include_TF_x.Steady = FALSE,
                 MS.include_dTF_x.PathDist = FALSE,
                 MS.include_dTF_x.Active = FALSE,
                 MS.include_dTF_x.Steady = FALSE,
                 MS.include_TF_y.PathDist = FALSE,
                 MS.include_TF_y.Active = FALSE,
                 MS.include_TF_y.Steady = FALSE,
                 MS.include_dTF_y.PathDist = FALSE,
                 MS.include_dTF_y.Active = FALSE,
                 MS.include_dTF_y.Steady = FALSE,
                 MS.include_TB_x.PathDist = FALSE,
                 MS.include_TB_x.Active = FALSE,
                 MS.include_TB_x.Steady = FALSE,
                 MS.include_dTB_x.PathDist = FALSE,
                 MS.include_dTB_x.Active = FALSE,
                 MS.include_dTB_x.Steady = FALSE,
                 MS.include_TB_y.PathDist = FALSE,
                 MS.include_TB_y.Active = FALSE,
                 MS.include_TB_y.Steady = FALSE,
                 MS.include_dTB_y.PathDist = FALSE,
                 MS.include_dTB_y.Active = FALSE,
                 MS.include_dTB_y.Steady = FALSE,
                 
                 MS.include_LipAperture = FALSE
                 ) {
  
  # Loading the acoustic data into R 
  acousticData <- read.delim(fbwFilePath,
                             header = FALSE) # Upload the acoustic FBW file.
  colnames(acousticData) <- c("Time", "F1","F2","F3","dB","1","2")
  
  acousticSegments <- openxlsx::createWorkbook()
  
  # Creating Bark Variables 
  acousticData$F1_bark <- emuR::bark(acousticData$F1*1000)
  acousticData$F2_bark <- emuR::bark(acousticData$F2*1000)
  acousticData$F3_bark <- emuR::bark(acousticData$F3*1000)
  
  # Loading F0 file into R, if specified to include 
  if (includeF0Data == TRUE) {
    F0Data <- read.delim(F0FilePath,
                         header = FALSE)
    F0Data %>%
      mutate(V2 = ifelse(V2 < 50,NA,V2)) ->
      F0Data
  }
  
  # Loading the decoupled movement data into R 
  movementData <- read.csv(file = decoupledMovementFilePath, header=TRUE, stringsAsFactors=FALSE)
  movementSegments <- openxlsx::createWorkbook()
  
  # Loading the target segments (lbl file) into R 
  targetIntervals <- read.table(lblFilePath,
                                header = FALSE,
                                sep = "\t") # Upload the segmented LBL file.
  targetIntervals %>%
    dplyr::rename(onsets = V1,
                  offsets = V2,
                  targetNames = V3) %>%
    dplyr::mutate(targetNames = base::make.unique(as.character(targetNames))) ->
    targetIntervals

  # Set-up & Defining On/Offsets 
  onsets_ms <- targetIntervals$onsets
  offsets_ms <- targetIntervals$offsets
  onsets_s <- onsets_ms/1000
  offsets_s <- offsets_ms/1000
  targetNames <- targetIntervals$targetNames
  
  # Analysis Loop 
      # Calculations and recording of each segment identified in the Target Intervals
      k <- 1
      mastersheet <- data.frame(k)
  
      while (k <= NROW(targetIntervals$targetNames)) {
        # * Sub and Target info reported to the mastersheet 
        mastersheet$k[k] <- k
        mastersheet$Group[k] <- paste(Group)
        mastersheet$Sub[k] <- Database.ID
        mastersheet$Speaker[k] <- Study.ID
        mastersheet$Target[k] <- paste(targetIntervals$targetNames[k])
        
        # * Finding the target segments 
        segOnset_s <- onsets_s[k] #Segment onset in Seconds
        segOffset_s <- offsets_s[k] #Segment offset in Seconds
        segOnset_ms <- onsets_ms[k] #Segment onset in Milliseconds
        segOffset_ms <- offsets_ms[k] #Segment offset in Milliseconds
        segmentMovement <- base::subset(movementData, movementData$time > segOnset_s & movementData$time < segOffset_s) #ID segment from kinematic data
        segmentAcoustic <- base::subset(acousticData, acousticData$Time > segOnset_ms & acousticData$Time < segOffset_ms) #ID segement from acoustic data
        segmentAcoustic %>%
          naniar::replace_with_na(replace = list(
            F1 = c(0),
            F2 = c(0),
            F3 = c(0))) ->
          segmentAcoustic
        
        if (include_F0Data == TRUE) {segmentF0 <- base::subset(F0Data, F0Data$V1 > segOnset_ms & F0Data$V1 < segOffset_ms)}
        
        # * Duration Measures 
        if (AS.include_NormalizeTime == TRUE) {segmentAcoustic$normTime <- as.numeric(scales::rescale(segmentAcoustic$Time,to = c(0,100)))}
        
        mastersheet$duration_ms[k] <- targetIntervals$offsets[k] - targetIntervals$onsets[k] # Segment duration (ms)
        mastersheet$duration_s[k] <- mastersheet$duration_ms[k]/1000
        mastersheet$onset_ms[k] <- segOnset_ms
        mastersheet$offset_ms[k] <- segOffset_ms
        mastersheet$onset_s[k] <- segOnset_s
        mastersheet$offset_s[k] <- segOffset_s
        mastersheet$midValue_ms[k] <- (mastersheet$duration_ms[k]/2) + targetIntervals$onsets[k] # Temporal midpoint of segment (ms)
        mastersheet$midValue_s[k] <- mastersheet$midValue_ms[k]/1000 # Temporal midpoint of segment (s)
        tempMid.aco <- segmentAcoustic$Time[base::which.min(base::abs((mastersheet$midValue_ms[k]) - segmentAcoustic$Time))]
        tempMid.kin <- base::which.min(base::abs((segmentAcoustic$Time[segmentAcoustic$Time == tempMid.aco] / 1000) - segmentMovement$time))
        mastersheet$tempMid_ms.Aco <- tempMid.aco
        mastersheet$tempMid_ms.Kin <- tempMid.kin
        # Temporal Midpoints
        if (include_F1.tempMid == TRUE) {mastersheet$F1_tempMid[k] <- segmentAcoustic$F1[segmentAcoustic$Time == tempMid.aco]}
        if (include_F2.tempMid == TRUE) {mastersheet$F2_tempMid[k] <- segmentAcoustic$F2[segmentAcoustic$Time == tempMid.aco]}
        if (include_F3.tempMid == TRUE) {mastersheet$F3_tempMid[k] <- segmentAcoustic$F3[segmentAcoustic$Time == tempMid.aco]}
        
        if (include_dTF.x_tempMid == TRUE) {mastersheet$dTF.x_tempMid[k] <- segmentMovement$decoupled_tf_x[tempMid.kin]}
        if (include_dTF.y_tempMid == TRUE) {mastersheet$dTF.y_tempMid[k] <- segmentMovement$decoupled_tf_y[tempMid.kin]}
        
        if (include_dTB.x_tempMid == TRUE) {mastersheet$dTB.x_tempMid[k] <- segmentMovement$decoupled_tb_x[tempMid.kin]}
        if (include_dTB.y_tempMid == TRUE) {mastersheet$dTB.y_tempMid[k] <- segmentMovement$decoupled_tb_y[tempMid.kin]}
        
        # AAVS
        if (include_AAVS == TRUE) {mastersheet$aavs[k] <- aavs(segmentAcoustic$F1_bark,segmentAcoustic$F2_bark)}
        
        # * Constriction Information
        segmentAcoustic$F2_F3constriction <- segmentAcoustic$F3 - segmentAcoustic$F2
        if (!all(is.na(segmentAcoustic$F2_F3constriction))){
          if(include_minF2_F3 == TRUE) {mastersheet$minF2_F3[k] <- ifelse(NROW(segmentAcoustic$F2_F3constriction) == 0 ,NA,base::min(segmentAcoustic$F2_F3constriction))}
          maxConstriction.aco <- base::which.min(segmentAcoustic$F3 - segmentAcoustic$F2)
          maxConstriction.kin <- base::which.min(base::abs((segmentAcoustic$Time[maxConstriction.aco] / 1000) - segmentMovement$time))
          if (include_TimeMaxCons == TRUE) {mastersheet$timeMaxCons[k] <- segmentAcoustic$Time[maxConstriction.aco]}
          if (include_NormTimeMaxCons == TRUE) {mastersheet$NormTimeMaxCons[k] <- segmentAcoustic$normTime[maxConstriction.aco]}
          
          # * Formants at Max Constriction
          if(include_F1_maxConstriction == TRUE) {mastersheet$F1_maxConstriction[k] <- segmentAcoustic$F1[maxConstriction.aco]}
          if(include_F2_maxConstriction == TRUE) {mastersheet$F2_maxConstriction[k] <- segmentAcoustic$F2[maxConstriction.aco]}
          if(include_F3_maxConstriction == TRUE) {mastersheet$F3_maxConstriction[k] <- segmentAcoustic$F3[maxConstriction.aco]}
          
          if (include_dTF.x_maxConstriction == TRUE) {mastersheet$dTF.x_maxConstriction[k] <- segmentMovement$decoupled_tf_x[maxConstriction.kin]}
          if (include_dTF.y_maxConstriction == TRUE) {mastersheet$dTF.y_maxConstriction[k] <- segmentMovement$decoupled_tf_y[maxConstriction.kin]}
        } else {
          if(include_minF2_F3 == TRUE) {mastersheet$minF2_F3[k] <- NA}
          maxConstriction.aco <- NA
          maxConstriction.kin <- NA
          if (include_TimeMaxCons == TRUE) {mastersheet$timeMaxCons[k] <- NA}
          if (include_NormTimeMaxCons == TRUE) {mastersheet$NormTimeMaxCons[k] <- NA}
          
          # * Formants at Max Constriction 
          mastersheet$F1_maxConstriction[k] <- NA
          mastersheet$F2_maxConstriction[k] <- NA
          mastersheet$F3_maxConstriction[k] <- NA
          
          mastersheet$dTF.x_maxConstriction[k] <- NA
          mastersheet$dTF.y_maxConstriction[k] <- NA
        }
        
        
        # Acoustic Measures
        # If the acoustic segment is too short, this segment makes the measures NA
        if (include_F0Data == TRUE) {
         if (is.null(segmentF0)) {
          # Pitch
           if (include_F0.mean == TRUE) {mastersheet$f0_mean[k] <- NA}
           if (include_F0.max == TRUE) {mastersheet$f0_max[k] <- NA}
           if (include_F0.min == TRUE) {mastersheet$f0_min[k] <- NA}
           if (include_F0.sd == TRUE) {mastersheet$f0_sd[k] <- NA}
           } else {
          # Pitch
          if (include_F0.mean == TRUE) {mastersheet$f0_mean[k] <- base::mean(segmentF0$V2,na.rm = TRUE)}
          if (include_F0.max == TRUE) {mastersheet$f0_max[k] <- base::ifelse(NROW(segmentF0) == 0,NA,max(segmentF0$V2, na.rm = TRUE))}
          if (include_F0.min == TRUE) {mastersheet$f0_min[k] <- base::ifelse(NROW(segmentF0) == 0,NA,min(segmentF0$V2, na.rm = TRUE))}
          if (include_F0.sd == TRUE) {mastersheet$f0_sd[k] <- stats::sd(segmentF0$V2, na.rm = TRUE)}
        }
       }
        if (NROW(segmentAcoustic) < 1) {
          # * F1 Measures
          if (include_F1.onset == TRUE) {mastersheet$F1_Onset.kHz[k] <- NA}
          if (include_F1.offset == TRUE) {mastersheet$F1_Offset.kHz[k] <- NA}
          if (include_F1.extent == TRUE) {mastersheet$F1_Extent.kHz[k] <- NA}
          if (include_F1.slope_Hz.ms == TRUE) {mastersheet$F1_Slope.Hz.ms[k] <- NA}
          if (include_F1.slope_2020.Hz.ms == TRUE) {mastersheet$F1_Slope_2020.Hz.ms[k] <- NA}
          
          # * F2 Measures
          if (include_F2.onset == TRUE) {mastersheet$F2_Onset.kHz[k] <- NA}
          if (include_F2.offset == TRUE) {mastersheet$F2_Offset.kHz[k] <- NA}
          if (include_F2.extent == TRUE) {mastersheet$F2_Extent.kHz[k] <- NA}
          if (include_F2.slope_Hz.ms == TRUE) {mastersheet$F2_Slope.Hz.ms[k] <- NA}
          if (include_F2.slope_2020.Hz.ms == TRUE) {mastersheet$F2_Slope_2020.Hz.ms[k] <- NA}
          if (include_F2.steadyDuration == TRUE) {mastersheet$F2_steadyDuration[k] <- NA}
          
          # * F3 Measures
          if (include_F3.onset == TRUE) {mastersheet$F3_Onset.kHz[k] <- NA}
          if (include_F3.offset == TRUE) {mastersheet$F3_Offset.kHz[k] <- NA}
          if (include_F3.extent == TRUE) {mastersheet$F3_Extent.kHz[k] <- NA}
          if (include_F3.slope_Hz.ms == TRUE) {mastersheet$F3_Slope.Hz.ms[k] <- NA}
          if (include_F3.slope_2020.Hz.ms == TRUE) {mastersheet$F3_Slope_2020.Hz.ms[k] <- NA}
          
          # * dB
          if (include_dB.mean == TRUE) {mastersheet$dB_mean[k] <- NA}
          if (include_dB.max == TRUE) {mastersheet$dB_max[k] <- NA}
          if (include_dB.min == TRUE) {mastersheet$dB_min[k] <- NA}
          if (include_dB.sd == TRUE) {mastersheet$dB_sd[k] <- NA}
          
        } else {
          # * F1 Measures
          F1 <- formantAnalysis(segmentAcoustic$Time,segmentAcoustic$F1)
          if (include_F1.onset == TRUE) {mastersheet$F1_Onset.kHz[k] <- F1$onset}
          if (include_F1.offset == TRUE) {mastersheet$F1_Offset.kHz[k] <- F1$offset}
          if (include_F1.extent == TRUE) {mastersheet$F1_Extent.kHz[k] <- F1$extent}
          if (include_F1.slope_Hz.ms == TRUE) {mastersheet$F1_Slope.Hz.ms[k] <- F1$slope}
          if (include_F1.slope_2020.Hz.ms == TRUE) {mastersheet$F1_Slope_2020.Hz.ms[k] <- F1$f_slope2020}
          
          if (NROW(segmentAcoustic) == NROW(F1$f_active)) {
            segmentAcoustic$F1_active <- F1$f_active
            segmentAcoustic$F1_steady <- F1$f_steady
          } else {
            segmentAcoustic$F1_active <- NA
            segmentAcoustic$F1_steady <- NA        
          }
          rm(F1)
          
          # * F2 Measures
          F2 <- formantAnalysis(segmentAcoustic$Time,segmentAcoustic$F2)
          if (include_F2.onset == TRUE) {mastersheet$F2_Onset.kHz[k] <- F2$onset}
          if (include_F2.offset == TRUE) {mastersheet$F2_Offset.kHz[k] <- F2$offset}
          if (include_F2.extent == TRUE) {mastersheet$F2_Extent.kHz[k] <- F2$extent}
          if (include_F2.slope_Hz.ms == TRUE) {mastersheet$F2_Slope.Hz.ms[k] <- F2$slope}
          if (include_F2.slope_2020.Hz.ms == TRUE) {mastersheet$F2_Slope_2020.Hz.ms[k] <- F2$f_slope2020}
          
          if (NROW(segmentAcoustic) == NROW(F2$f_active)) {
            segmentAcoustic$F2_active <- F2$f_active
            segmentAcoustic$F2_steady <- F2$f_steady
          } else {
            segmentAcoustic$F2_active <- NA
            segmentAcoustic$F2_steady <- NA        
          }
          if (include_F2.steadyDuration == TRUE) {mastersheet$F2_steadyDuration[k] <- F2$f_steadyDuration}
          rm(F2)
          
          # * F3 Measures 
          F3 <- formantAnalysis(segmentAcoustic$Time,segmentAcoustic$F3)
          if (include_F3.onset == TRUE) {mastersheet$F3_Onset.kHz[k] <- F3$onset}
          if (include_F3.offset == TRUE) {mastersheet$F3_Offset.kHz[k] <- F3$offset}
          if (include_F3.extent == TRUE) {mastersheet$F3_Extent.kHz[k] <- F3$extent}
          if (include_F3.slope_Hz.ms == TRUE) {mastersheet$F3_Slope.Hz.ms[k] <- F3$slope}
          if (include_F3.slope_2020.Hz.ms == TRUE) {mastersheet$F3_Slope_2020.Hz.ms[k] <- F3$f_slope2020}
          
          if (NROW(segmentAcoustic) == NROW(F3$f_active)) {
            segmentAcoustic$F3_active <- F3$f_active
            segmentAcoustic$F3_steady <- F3$f_steady
          } else {
            segmentAcoustic$F3_active <- NA
            segmentAcoustic$F3_steady <- NA        
          }
          rm(F3)
          
          # * dB 
          dBData <- subset(segmentAcoustic,segmentAcoustic$dB > -100)
          dBData$dB <- dBData$dB + 100
          segmentAcoustic$dB <- segmentAcoustic$dB + 100
          if (include_dB.mean == TRUE) {mastersheet$dB_mean[k] <- base::mean(dBData$dB, na.rm = TRUE)}
          if (include_dB.max == TRUE) {mastersheet$dB_max[k] <- base::max(dBData$dB, na.rm = TRUE)}
          if (include_dB.min == TRUE) {mastersheet$dB_min[k] <- base::min(dBData$dB, na.rm = TRUE)}
          if (include_dB.sd == TRUE) {mastersheet$dB_sd[k] <- stats::sd(dBData$dB, na.rm = TRUE)}
          rm(dBData)
          
        } # For segments of normal length, these are the acoustic measurement calculations
        
        # Acoustic Segments
        if (AS.include_F1.Bark == FALSE) {segmentAcoustic %>% dplyr::select(!F1_bark) -> segmentAcoustic}
        if (AS.include_F1.Active == FALSE) {segmentAcoustic %>% dplyr::select(!F1_active) -> segmentAcoustic}
        if (AS.include_F1.Steady == FALSE) {segmentAcoustic %>% dplyr::select(!F1_steady) -> segmentAcoustic}
        if (AS.include_F2.Bark == FALSE) {segmentAcoustic %>% dplyr::select(!F2_bark) -> segmentAcoustic}
        if (AS.include_F2.Active == FALSE) {segmentAcoustic %>% dplyr::select(!F2_active) -> segmentAcoustic}
        if (AS.include_F2.Steady == FALSE) {segmentAcoustic %>% dplyr::select(!F2_steady) -> segmentAcoustic}
        if (AS.include_F3.Bark == FALSE) {segmentAcoustic %>% dplyr::select(!F3_bark) -> segmentAcoustic}
        if (AS.include_F3.Active == FALSE) {segmentAcoustic %>% dplyr::select(!F3_active) -> segmentAcoustic}
        if (AS.include_F3.Steady == FALSE) {segmentAcoustic %>% dplyr::select(!F3_steady) -> segmentAcoustic}
        if (AS.include_F2_F3constriction == FALSE) {segmentAcoustic %>% dplyr::select(!F2_F3constriction) -> segmentAcoustic}
        
        # Saving Acoustic Segment to List
        sheetname <- as.character(targetIntervals$targetNames)[k]
        openxlsx::addWorksheet(acousticSegments, sheetname)
        openxlsx::writeData(acousticSegments, sheet = sheetname, x = segmentAcoustic)
        
        # Kinematic Measures 
        # If the movement segment is too short, this segment makes the measures NA
        if (NROW(segmentMovement) <= 1) {
          # NA TF Measures
          {
            # * X
            if (include_TF.x_Onset == TRUE) {mastersheet$TF.x_Onset[k] <- NA}
            if (include_TF.x_Offset == TRUE) {mastersheet$TF.x_Offset[k] <- NA}
            if (include_TF.x_Distance == TRUE) {mastersheet$TF.x_Distance[k] <- NA}
            if (include_TF.x_Displacement == TRUE) {mastersheet$TF.x_Displacement[k] <- NA}
            if (include_TF.x_SpeedAvg == TRUE) {mastersheet$TF.x_SpeedAvg[k] <- NA}
            if (include_TF.x_SpeedMax == TRUE) {mastersheet$TF.x_SpeedMax[k] <- NA}
            if (include_TF.x_steadyDuration == TRUE) {mastersheet$TF.x_steadyDuration[k] <- NA}
            
            # * X (Decoupled)
            if (include_dTF.x_Onset == TRUE) {mastersheet$dTF.x_Onset[k] <- NA}
            if (include_dTF.x_Offset == TRUE) {mastersheet$dTF.x_Offset[k] <- NA}
            if (include_dTF.x_Distance == TRUE) {mastersheet$dTF.x_Distance[k] <- NA}
            if (include_dTF.x_Displacement == TRUE) {mastersheet$dTF.x_Displacement[k] <- NA}
            if (include_dTF.x_SpeedAvg == TRUE) {mastersheet$dTF.x_SpeedAvg[k] <- NA}
            if (include_dTF.x_SpeedMax == TRUE) {mastersheet$dTF.x_SpeedMax[k] <- NA}
            if (include_dTF.x_steadyDuration == TRUE) {mastersheet$dTF.x_steadyDuration[k] <- NA}
            
            # * Y
            if (include_TF.y_Onset == TRUE) {mastersheet$TF.y_Onset[k] <- NA}
            if (include_TF.y_Offset == TRUE) {mastersheet$TF.y_Offset[k] <- NA}
            if (include_TF.y_Distance == TRUE) {mastersheet$TF.y_Distance[k] <- NA}
            if (include_TF.y_Displacement == TRUE) {mastersheet$TF.y_Displacement[k] <- NA}
            if (include_TF.y_SpeedAvg == TRUE) {mastersheet$TF.y_SpeedAvg[k] <- NA}
            if (include_TF.y_SpeedMax == TRUE) {mastersheet$TF.y_SpeedMax[k] <- NA}
            if (include_TF.y_steadyDuration == TRUE) {mastersheet$TF.y_steadyDuration[k] <- NA}
            
            # * Y (Decoupled)
            if (include_dTF.y_Onset == TRUE) {mastersheet$dTF.y_Onset[k] <- NA}
            if (include_dTF.y_Offset == TRUE) {mastersheet$dTF.y_Offset[k] <- NA}
            if (include_dTF.y_Distance == TRUE) {mastersheet$dTF.y_Distance[k] <- NA}
            if (include_dTF.y_Displacement == TRUE) {mastersheet$dTF.y_Displacement[k] <- NA}
            if (include_dTF.y_SpeedAvg == TRUE) {mastersheet$dTF.y_SpeedAvg[k] <- NA}
            if (include_dTF.y_SpeedMax == TRUE) {mastersheet$dTF.y_SpeedMax[k] <- NA}
            if (include_dTF.y_steadyDuration == TRUE) {mastersheet$dTF.y_steadyDuration[k] <- NA}
            
            # * XY
            if (include_TF.xy_Distance == TRUE) {mastersheet$TF.xy_Distance[k] <- NA}
            if (include_TF.xy_Displacement == TRUE) {mastersheet$TF.xy_Displacement[k] <- NA}
            if (include_TF.xy_SpeedAvg == TRUE) {mastersheet$TF.xy_SpeedAvg[k] <- NA}
            if (include_TF.xy_SpeedMax == TRUE) {mastersheet$TF.xy_SpeedMax[k] <- NA}
            
            # * XY (Decoupled)
            if (include_dTF.xy_Distance == TRUE) {mastersheet$dTF.xy_Distance[k] <- NA}
            if (include_dTF.xy_Displacement == TRUE) {mastersheet$dTF.xy_Displacement[k] <- NA}
            if (include_dTF.xy_SpeedAvg == TRUE) {mastersheet$dTF.xy_SpeedAvg[k] <- NA}
            if (include_dTF.xy_SpeedMax == TRUE) {mastersheet$dTF.xy_SpeedMax[k] <- NA}
            
            # * Convex Hull - TF
            if (include_TF.cHull == TRUE) {mastersheet$TF.cHull[k] <- NA}
            if (include_dTF.cHull == TRUE) {mastersheet$dTF.cHull[k] <- NA}
            
            mastersheet$dTF.x_maxConstriction[k] <- NA
            mastersheet$dTF.y_maxConstriction[k] <- NA
          }
          # NA TB Measures
          {
            # * X
            if (include_TB.x_Onset == TRUE) {mastersheet$TB.x_Onset[k] <- NA}
            if (include_TB.x_Offset == TRUE) {mastersheet$TB.x_Offset[k] <- NA}
            if (include_TB.x_Distance == TRUE) {mastersheet$TB.x_Distance[k] <- NA}
            if (include_TB.x_Displacement == TRUE) {mastersheet$TB.x_Displacement[k] <- NA}
            if (include_TB.x_SpeedAvg == TRUE) {mastersheet$TB.x_SpeedAvg[k] <- NA}
            if (include_TB.x_SpeedMax == TRUE) {mastersheet$TB.x_SpeedMax[k] <- NA}
            if (include_TB.x_steadyDuration == TRUE) {mastersheet$TB.x_steadyDuration[k] <- NA}
            
            # * X (Decoupled)
            if (include_dTB.x_Onset == TRUE) {mastersheet$dTB.x_Onset[k] <- NA}
            if (include_dTB.x_Offset == TRUE) {mastersheet$dTB.x_Offset[k] <- NA}
            if (include_dTB.x_Distance == TRUE) {mastersheet$dTB.x_Distance[k] <- NA}
            if (include_dTB.x_Displacement == TRUE) {mastersheet$dTB.x_Displacement[k] <- NA}
            if (include_dTB.x_SpeedAvg == TRUE) {mastersheet$dTB.x_SpeedAvg[k] <- NA}
            if (include_dTB.x_SpeedMax == TRUE) {mastersheet$dTB.x_SpeedMax[k] <- NA}
            if (include_dTB.x_steadyDuration == TRUE) {mastersheet$dTB.x_steadyDuration[k] <- NA}
            
            # * Y
            if (include_TB.y_Onset == TRUE) {mastersheet$TB.y_Onset[k] <- NA}
            if (include_TB.y_Offset == TRUE) {mastersheet$TB.y_Offset[k] <- NA}
            if (include_TB.y_Distance == TRUE) {mastersheet$TB.y_Distance[k] <- NA}
            if (include_TB.y_Displacement == TRUE) {mastersheet$TB.y_Displacement[k] <- NA}
            if (include_TB.y_SpeedAvg == TRUE) {mastersheet$TB.y_SpeedAvg[k] <- NA}
            if (include_TB.y_SpeedMax == TRUE) {mastersheet$TB.y_SpeedMax[k] <- NA}
            if (include_TB.y_steadyDuration == TRUE) {mastersheet$TB.y_steadyDuration[k] <- NA}
            
            # * Y (Decoupled)
            if (include_dTB.y_Onset == TRUE) {mastersheet$dTB.y_Onset[k] <- NA}
            if (include_dTB.y_Offset == TRUE) {mastersheet$dTB.y_Offset[k] <- NA}
            if (include_dTB.y_Distance == TRUE) {mastersheet$dTB.y_Distance[k] <- NA}
            if (include_dTB.y_Displacement == TRUE) {mastersheet$dTB.y_Displacement[k] <- NA}
            if (include_dTB.y_SpeedAvg == TRUE) {mastersheet$dTB.y_SpeedAvg[k] <- NA}
            if (include_dTB.y_SpeedMax == TRUE) {mastersheet$dTB.y_SpeedMax[k] <- NA}
            if (include_dTB.y_steadyDuration == TRUE) {mastersheet$dTB.y_steadyDuration[k] <- NA}
            
            # * XY
            if (include_TB.xy_Distance == TRUE) {mastersheet$TB.xy_Distance[k] <- NA}
            if (include_TB.xy_Displacement == TRUE) {mastersheet$TB.xy_Displacement[k] <- NA}
            if (include_TB.xy_SpeedAvg == TRUE) {mastersheet$TB.xy_SpeedAvg[k] <- NA}
            if (include_TB.xy_SpeedMax == TRUE) {mastersheet$TB.xy_SpeedMax[k] <- NA}
            
            # * XY (Decoupled)
            if (include_dTB.xy_Distance == TRUE) {mastersheet$dTB.xy_Distance[k] <- NA}
            if (include_dTB.xy_Displacement == TRUE) {mastersheet$dTB.xy_Displacement[k] <- NA}
            if (include_dTB.xy_SpeedAvg == TRUE) {mastersheet$dTB.xy_SpeedAvg[k] <- NA}
            if (include_dTB.xy_SpeedMax == TRUE) {mastersheet$dTB.xy_SpeedMax[k] <- NA}
            
            # * Convex Hull - TB
            if (include_TB.cHull == TRUE) {mastersheet$TB.cHull[k] <- NA}
            if (include_dTB.cHull == TRUE) {mastersheet$dTB.cHull[k] <- NA}
            
            mastersheet$dTB.x_maxConstriction[k] <- NA
            mastersheet$dTB.y_maxConstriction[k] <- NA
          }
          # NA UL Measures
          {
            # * X
            if (include_UL.x_Onset == TRUE) {mastersheet$UL.x_Onset[k] <- NA}
            if (include_UL.x_Offset == TRUE) {mastersheet$UL.x_Offset[k] <- NA}
            if (include_UL.x_Distance == TRUE) {mastersheet$UL.x_Distance[k] <- NA}
            if (include_UL.x_Displacement == TRUE) {mastersheet$UL.x_Displacement[k] <- NA}
            if (include_UL.x_SpeedAvg == TRUE) {mastersheet$UL.x_SpeedAvg[k] <- NA}
            if (include_UL.x_SpeedMax == TRUE) {mastersheet$UL.x_SpeedMax[k] <- NA}
            
            # * Y
            if (include_UL.y_Onset == TRUE) {mastersheet$UL.y_Onset[k] <- NA}
            if (include_UL.y_Offset == TRUE) {mastersheet$UL.y_Offset[k] <- NA}
            if (include_UL.y_Distance == TRUE) {mastersheet$UL.y_Distance[k] <- NA}
            if (include_UL.y_Displacement == TRUE) {mastersheet$UL.y_Displacement[k] <- NA}
            if (include_UL.y_SpeedAvg == TRUE) {mastersheet$UL.y_SpeedAvg[k] <- NA}
            if (include_UL.y_SpeedMax == TRUE) {mastersheet$UL.y_SpeedMax[k] <- NA}
            
            # * XY
            if (include_UL.xy_Distance == TRUE) {mastersheet$UL.xy_Distance[k] <- NA}
            if (include_UL.xy_Displacement == TRUE) {mastersheet$UL.xy_Displacement[k] <- NA}
            if (include_UL.xy_SpeedAvg == TRUE) {mastersheet$UL.xy_SpeedAvg[k] <- NA}
            if (include_UL.xy_SpeedMax == TRUE) {mastersheet$UL.xy_SpeedMax[k] <- NA}
            
            # * Convex Hull - UL
            if (include_UL.cHull == TRUE) {mastersheet$UL.cHull[k] <- NA}
            
          }
          # NA LL Measures
          {
            # * X
            if (include_LL.x_Onset == TRUE) {mastersheet$LL.x_Onset[k] <- NA}
            if (include_LL.x_Offset == TRUE) {mastersheet$LL.x_Offset[k] <- NA}
            if (include_LL.x_Distance == TRUE) {mastersheet$LL.x_Distance[k] <- NA}
            if (include_LL.x_Displacement == TRUE) {mastersheet$LL.x_Displacement[k] <- NA}
            if (include_LL.x_SpeedAvg == TRUE) {mastersheet$LL.x_SpeedAvg[k] <- NA}
            if (include_LL.x_SpeedMax == TRUE) {mastersheet$LL.x_SpeedMax[k] <- NA}
            
            # * X (Decoupled)
            if (include_dLL.x_Onset == TRUE) {mastersheet$dLL.x_Onset[k] <- NA}
            if (include_dLL.x_Offset == TRUE) {mastersheet$dLL.x_Offset[k] <- NA}
            if (include_dLL.x_Distance == TRUE) {mastersheet$dLL.x_Distance[k] <- NA}
            if (include_dLL.x_Displacement == TRUE) {mastersheet$dLL.x_Displacement[k] <- NA}
            if (include_dLL.x_SpeedAvg == TRUE) {mastersheet$dLL.x_SpeedAvg[k] <- NA}
            if (include_dLL.x_SpeedMax == TRUE) {mastersheet$dLL.x_SpeedMax[k] <- NA}
            
            # * Y
            if (include_LL.y_Onset == TRUE) {mastersheet$LL.y_Onset[k] <- NA}
            if (include_LL.y_Offset == TRUE) {mastersheet$LL.y_Offset[k] <- NA}
            if (include_LL.y_Distance == TRUE) {mastersheet$LL.y_Distance[k] <- NA}
            if (include_LL.y_Displacement == TRUE) {mastersheet$LL.y_Displacement[k] <- NA}
            if (include_LL.y_SpeedAvg == TRUE) {mastersheet$LL.y_SpeedAvg[k] <- NA}
            if (include_LL.y_SpeedMax == TRUE) {mastersheet$LL.y_SpeedMax[k] <- NA}
            
            # * Y (Decoupled)
            if (include_dLL.y_Onset == TRUE) {mastersheet$dLL.y_Onset[k] <- NA}
            if (include_dLL.y_Offset == TRUE) {mastersheet$dLL.y_Offset[k] <- NA}
            if (include_dLL.y_Distance == TRUE) {mastersheet$dLL.y_Distance[k] <- NA}
            if (include_dLL.y_Displacement == TRUE) {mastersheet$dLL.y_Displacement[k] <- NA}
            if (include_dLL.y_SpeedAvg == TRUE) {mastersheet$dLL.y_SpeedAvg[k] <- NA}
            if (include_dLL.y_SpeedMax == TRUE) {mastersheet$dLL.y_SpeedMax[k] <- NA}
            
            # * XY
            if (include_LL.xy_Distance == TRUE) {mastersheet$LL.xy_Distance[k] <- NA}
            if (include_LL.xy_Displacement == TRUE) {mastersheet$LL.xy_Displacement[k] <- NA}
            if (include_LL.xy_SpeedAvg == TRUE) {mastersheet$LL.xy_SpeedAvg[k] <- NA}
            if (include_LL.xy_SpeedMax == TRUE) {mastersheet$LL.xy_SpeedMax[k] <- NA}
            
            # * XY (Decoupled)
            if (include_dLL.xy_Distance == TRUE) {mastersheet$dLL.xy_Distance[k] <- NA}
            if (include_dLL.xy_Displacement == TRUE) {mastersheet$dLL.xy_Displacement[k] <- NA}
            if (include_dLL.xy_SpeedAvg == TRUE) {mastersheet$dLL.xy_SpeedAvg[k] <- NA}
            if (include_dLL.xy_SpeedMax == TRUE) {mastersheet$dLL.xy_SpeedMax[k] <- NA}
            
            # * Convex Hull - LL
            if (include_LL.cHull == TRUE) {mastersheet$LL.cHull[k] <- NA}
            if (include_dLL.cHull == TRUE) {mastersheet$dLL.cHull[k] <- NA}
            
            mastersheet$dLL.x_maxConstriction[k] <- NA
            mastersheet$dLL.y_maxConstriction[k] <- NA
          }
          # NA Jaw Measures
          {
            # * X
            mastersheet$Jaw.x_Onset[k] <- NA
            mastersheet$Jaw.x_Offset[k] <- NA
            mastersheet$Jaw.x_Distance[k] <- NA
            mastersheet$Jaw.x_Displacement[k] <- NA
            mastersheet$Jaw.x_SpeedAvg[k] <- NA
            mastersheet$Jaw.x_SpeedMax[k] <- NA
            
            # * Y
            mastersheet$Jaw.y_Onset[k] <- NA
            mastersheet$Jaw.y_Offset[k] <- NA
            mastersheet$Jaw.y_Distance[k] <- NA
            mastersheet$Jaw.y_Displacement[k] <- NA
            mastersheet$Jaw.y_SpeedAvg[k] <- NA
            mastersheet$Jaw.y_SpeedMax[k] <- NA
            
            # * XY
            mastersheet$Jaw.xy_Distance[k] <- NA
            mastersheet$Jaw.xy_Displacement[k] <- NA
            mastersheet$Jaw.xy_SpeedAvg[k] <- NA
            mastersheet$Jaw.xy_SpeedMax[k] <- NA
            
            # * Convex Hull - Jaw
            mastersheet$Jaw.cHull[k] <- NA
          }
          
        } else {
          # TF Measures 
          {
            # * X 
            tf.x <- movementAnalysis_1D(segmentMovement$time,segmentMovement$tf_x)
            mastersheet$TF.x_Onset[k] <- tf.x$onset
            mastersheet$TF.x_Offset[k] <- tf.x$offset
            mastersheet$TF.x_Distance[k] <- tf.x$distance
            mastersheet$TF.x_Displacement[k] <- tf.x$displacement
            mastersheet$TF.x_SpeedAvg[k] <- tf.x$speed_avg
            mastersheet$TF.x_SpeedMax[k] <- tf.x$speed_max
            
            segmentMovement$TF_x.PathDist <- tf.x$pathDist
            segmentMovement$TF_x.Active <- tf.x$active
            segmentMovement$TF_x.Steady <- tf.x$steady
            mastersheet$TF.x_steadyDuration[k] <- tf.x$steadyDuration
            
            rm(tf.x)
            
            # * X (Decoupled) 
            dtf.x <- movementAnalysis_1D(segmentMovement$time,segmentMovement$decoupled_tf_x)
            mastersheet$dTF.x_Onset[k] <- dtf.x$onset
            mastersheet$dTF.x_Offset[k] <- dtf.x$offset
            mastersheet$dTF.x_Distance[k] <- dtf.x$distance
            mastersheet$dTF.x_Displacement[k] <- dtf.x$displacement
            mastersheet$dTF.x_SpeedAvg[k] <- dtf.x$speed_avg
            mastersheet$dTF.x_SpeedMax[k] <- dtf.x$speed_max
            
            segmentMovement$dTF_x.PathDist <- dtf.x$pathDist
            segmentMovement$dTF_x.Active <- dtf.x$active
            segmentMovement$dTF_x.Steady <- dtf.x$steady
            mastersheet$dTF.x_steadyDuration[k] <- dtf.x$steadyDuration
            
            rm(dtf.x)
            
            # * Y 
            tf.y <- movementAnalysis_1D(segmentMovement$time,segmentMovement$tf_y)
            mastersheet$TF.y_Onset[k] <- tf.y$onset
            mastersheet$TF.y_Offset[k] <- tf.y$offset
            mastersheet$TF.y_Distance[k] <- tf.y$distance
            mastersheet$TF.y_Displacement[k] <- tf.y$displacement
            mastersheet$TF.y_SpeedAvg[k] <- tf.y$speed_avg
            mastersheet$TF.y_SpeedMax[k] <- tf.y$speed_max
            
            segmentMovement$TF_y.PathDist <- tf.y$pathDist
            segmentMovement$TF_y.Active <- tf.y$active
            segmentMovement$TF_y.Steady <- tf.y$steady
            mastersheet$TF.y_steadyDuration[k] <- tf.y$steadyDuration
            
            rm(tf.y)
            
            # * Y (Decoupled) 
            dtf.y <- movementAnalysis_1D(segmentMovement$time,segmentMovement$decoupled_tf_y)
            mastersheet$dTF.y_Onset[k] <- dtf.y$onset
            mastersheet$dTF.y_Offset[k] <- dtf.y$offset
            mastersheet$dTF.y_Distance[k] <- dtf.y$distance
            mastersheet$dTF.y_Displacement[k] <- dtf.y$displacement
            mastersheet$dTF.y_SpeedAvg[k] <- dtf.y$speed_avg
            mastersheet$dTF.y_SpeedMax[k] <- dtf.y$speed_max
            
            segmentMovement$dTF_y.PathDist <- dtf.y$pathDist
            segmentMovement$dTF_y.Active <- dtf.y$active
            segmentMovement$dTF_y.Steady <- dtf.y$steady
            mastersheet$dTF.y_steadyDuration[k] <- dtf.y$steadyDuration
            
            rm(dtf.y)
            
            # * XY  -
            tf.xy <- movementAnalysis_2D(segmentMovement$time,segmentMovement$tf_x,segmentMovement$tf_y)
            mastersheet$TF.xy_Distance[k] <- tf.xy$distance
            mastersheet$TF.xy_Displacement[k] <- tf.xy$displacement
            mastersheet$TF.xy_SpeedAvg[k] <- tf.xy$speed_avg
            mastersheet$TF.xy_SpeedMax[k] <- tf.xy$speed_max
            
            rm(tf.xy)
            
            # * XY (Decoupled)  -
            dtf.xy <- movementAnalysis_2D(segmentMovement$time,segmentMovement$decoupled_tf_x,segmentMovement$decoupled_tf_y)
            mastersheet$dTF.xy_Distance[k] <- dtf.xy$distance
            mastersheet$dTF.xy_Displacement[k] <- dtf.xy$displacement
            mastersheet$dTF.xy_SpeedAvg[k] <- dtf.xy$speed_avg
            mastersheet$dTF.xy_SpeedMax[k] <- dtf.xy$speed_max
            
            rm(dtf.xy)
            
            # * Convex Hull - TF 
            if (include_TF.cHull == TRUE) {mastersheet$TF.cHull[k] <- cHull(segmentMovement$tf_x,segmentMovement$tf_y)}
            if (include_dTF.cHull == TRUE) {mastersheet$dTF.cHull[k] <- cHull(segmentMovement$decoupled_tf_x,segmentMovement$decoupled_tf_y)}
          }
          # TB Measures 
          {
            # * X 
            tb.x <- movementAnalysis_1D(segmentMovement$time,segmentMovement$tb_x)
            mastersheet$TB.x_Onset[k] <- tb.x$onset
            mastersheet$TB.x_Offset[k] <- tb.x$offset
            mastersheet$TB.x_Distance[k] <- tb.x$distance
            mastersheet$TB.x_Displacement[k] <- tb.x$displacement
            mastersheet$TB.x_SpeedAvg[k] <- tb.x$speed_avg
            mastersheet$TB.x_SpeedMax[k] <- tb.x$speed_max
            
            segmentMovement$TB_x.PathDist <- tb.x$pathDist
            segmentMovement$TB_x.Active <- tb.x$active
            segmentMovement$TB_x.Steady <- tb.x$steady
            mastersheet$TB.x_steadyDuration[k] <- tb.x$steadyDuration
            
            rm(tb.x)
            
            # * X (Decoupled) 
            dtb.x <- movementAnalysis_1D(segmentMovement$time,segmentMovement$decoupled_tb_x)
            mastersheet$dTB.x_Onset[k] <- dtb.x$onset
            mastersheet$dTB.x_Offset[k] <- dtb.x$offset
            mastersheet$dTB.x_Distance[k] <- dtb.x$distance
            mastersheet$dTB.x_Displacement[k] <- dtb.x$displacement
            mastersheet$dTB.x_SpeedAvg[k] <- dtb.x$speed_avg
            mastersheet$dTB.x_SpeedMax[k] <- dtb.x$speed_max
            
            segmentMovement$dTB_x.PathDist <- dtb.x$pathDist
            segmentMovement$dTB_x.Active <- dtb.x$active
            segmentMovement$dTB_x.Steady <- dtb.x$steady
            mastersheet$dTB.x_steadyDuration[k] <- dtb.x$steadyDuration
            
            rm(dtb.x)
            
            # * Y 
            tb.y <- movementAnalysis_1D(segmentMovement$time,segmentMovement$tb_y)
            mastersheet$TB.y_Onset[k] <- tb.y$onset
            mastersheet$TB.y_Offset[k] <- tb.y$offset
            mastersheet$TB.y_Distance[k] <- tb.y$distance
            mastersheet$TB.y_Displacement[k] <- tb.y$displacement
            mastersheet$TB.y_SpeedAvg[k] <- tb.y$speed_avg
            mastersheet$TB.y_SpeedMax[k] <- tb.y$speed_max
            
            segmentMovement$TB_y.PathDist <- tb.y$pathDist
            segmentMovement$TB_y.Active <- tb.y$active
            segmentMovement$TB_y.Steady <- tb.y$steady
            mastersheet$TB.y_steadyDuration[k] <- tb.y$steadyDuration
            
            rm(tb.y)
            
            # * Y (Decoupled) 
            dtb.y <- movementAnalysis_1D(segmentMovement$time,segmentMovement$decoupled_tb_y)
            mastersheet$dTB.y_Onset[k] <- dtb.y$onset
            mastersheet$dTB.y_Offset[k] <- dtb.y$offset
            mastersheet$dTB.y_Distance[k] <- dtb.y$distance
            mastersheet$dTB.y_Displacement[k] <- dtb.y$displacement
            mastersheet$dTB.y_SpeedAvg[k] <- dtb.y$speed_avg
            mastersheet$dTB.y_SpeedMax[k] <- dtb.y$speed_max
            
            segmentMovement$dTB_y.PathDist <- dtb.y$pathDist
            segmentMovement$dTB_y.Active <- dtb.y$active
            segmentMovement$dTB_y.Steady <- dtb.y$steady
            mastersheet$dTB.y_steadyDuration[k] <- dtb.y$steadyDuration
            
            rm(dtb.y)
            
            # * XY  -
            tb.xy <- movementAnalysis_2D(segmentMovement$time,segmentMovement$tb_x,segmentMovement$tb_y)
            mastersheet$TB.xy_Distance[k] <- tb.xy$distance
            mastersheet$TB.xy_Displacement[k] <- tb.xy$displacement
            mastersheet$TB.xy_SpeedAvg[k] <- tb.xy$speed_avg
            mastersheet$TB.xy_SpeedMax[k] <- tb.xy$speed_max
            
            rm(tb.xy)
            
            # * XY (Decoupled)  -
            dtb.xy <- movementAnalysis_2D(segmentMovement$time,segmentMovement$decoupled_tb_x,segmentMovement$decoupled_tb_y)
            mastersheet$dTB.xy_Distance[k] <- dtb.xy$distance
            mastersheet$dTB.xy_Displacement[k] <- dtb.xy$displacement
            mastersheet$dTB.xy_SpeedAvg[k] <- dtb.xy$speed_avg
            mastersheet$dTB.xy_SpeedMax[k] <- dtb.xy$speed_max
            
            rm(dtb.xy)
            
            # * Convex Hull - TB 
            if (include_TB.cHull == TRUE) {mastersheet$TB.cHull[k] <- cHull(segmentMovement$tb_x,segmentMovement$tb_y)}
            if (include_dTB.cHull == TRUE) {mastersheet$dTB.cHull[k] <- cHull(segmentMovement$decoupled_tb_x,segmentMovement$decoupled_tb_y)}
          }
          # UL Measures 
          {
            # * X 
            ul.x <- movementAnalysis_1D(segmentMovement$time,segmentMovement$ul_x)
            mastersheet$UL.x_Onset[k] <- ul.x$onset
            mastersheet$UL.x_Offset[k] <- ul.x$offset
            mastersheet$UL.x_Distance[k] <- ul.x$distance
            mastersheet$UL.x_Displacement[k] <- ul.x$displacement
            mastersheet$UL.x_SpeedAvg[k] <- ul.x$speed_avg
            mastersheet$UL.x_SpeedMax[k] <- ul.x$speed_max
            
            rm(ul.x)
            
            # * Y 
            ul.y <- movementAnalysis_1D(segmentMovement$time,segmentMovement$ul_y)
            mastersheet$UL.y_Onset[k] <- ul.y$onset
            mastersheet$UL.y_Offset[k] <- ul.y$offset
            mastersheet$UL.y_Distance[k] <- ul.y$distance
            mastersheet$UL.y_Displacement[k] <- ul.y$displacement
            mastersheet$UL.y_SpeedAvg[k] <- ul.y$speed_avg
            mastersheet$UL.y_SpeedMax[k] <- ul.y$speed_max
            
            rm(ul.y)
            
            # * XY  -
            ul.xy <- movementAnalysis_2D(segmentMovement$time,segmentMovement$ul_x,segmentMovement$ul_y)
            mastersheet$UL.xy_Distance[k] <- ul.xy$distance
            mastersheet$UL.xy_Displacement[k] <- ul.xy$displacement
            mastersheet$UL.xy_SpeedAvg[k] <- ul.xy$speed_avg
            mastersheet$UL.xy_SpeedMax[k] <- ul.xy$speed_max
            
            rm(ul.xy)
            
            # * Convex Hull - UL 
            if (include_UL.cHull == TRUE) {mastersheet$UL.cHull[k] <- cHull(segmentMovement$ul_x,segmentMovement$ul_y)}
          }
          # LL Measures 
          {
            # * X 
            ll.x <- movementAnalysis_1D(segmentMovement$time,segmentMovement$ll_x)
            mastersheet$LL.x_Onset[k] <- ll.x$onset
            mastersheet$LL.x_Offset[k] <- ll.x$offset
            mastersheet$LL.x_Distance[k] <- ll.x$distance
            mastersheet$LL.x_Displacement[k] <- ll.x$displacement
            mastersheet$LL.x_SpeedAvg[k] <- ll.x$speed_avg
            mastersheet$LL.x_SpeedMax[k] <- ll.x$speed_max
            
            rm(ll.x)
            
            # * X (Decoupled) 
            dll.x <- movementAnalysis_1D(segmentMovement$time,segmentMovement$decoupled_ll_x)
            mastersheet$dLL.x_Onset[k] <- dll.x$onset
            mastersheet$dLL.x_Offset[k] <- dll.x$offset
            mastersheet$dLL.x_Distance[k] <- dll.x$distance
            mastersheet$dLL.x_Displacement[k] <- dll.x$displacement
            mastersheet$dLL.x_SpeedAvg[k] <- dll.x$speed_avg
            mastersheet$dLL.x_SpeedMax[k] <- dll.x$speed_max
            
            rm(dll.x)
            
            # * Y 
            ll.y <- movementAnalysis_1D(segmentMovement$time,segmentMovement$ll_y)
            mastersheet$LL.y_Onset[k] <- ll.y$onset
            mastersheet$LL.y_Offset[k] <- ll.y$offset
            mastersheet$LL.y_Distance[k] <- ll.y$distance
            mastersheet$LL.y_Displacement[k] <- ll.y$displacement
            mastersheet$LL.y_SpeedAvg[k] <- ll.y$speed_avg
            mastersheet$LL.y_SpeedMax[k] <- ll.y$speed_max
            
            rm(ll.y)
            
            # * Y (Decoupled) 
            dll.y <- movementAnalysis_1D(segmentMovement$time,segmentMovement$decoupled_ll_y)
            mastersheet$dLL.y_Onset[k] <- dll.y$onset
            mastersheet$dLL.y_Offset[k] <- dll.y$offset
            mastersheet$dLL.y_Distance[k] <- dll.y$distance
            mastersheet$dLL.y_Displacement[k] <- dll.y$displacement
            mastersheet$dLL.y_SpeedAvg[k] <- dll.y$speed_avg
            mastersheet$dLL.y_SpeedMax[k] <- dll.y$speed_max
            
            rm(dll.y)
            
            # * XY  -
            ll.xy <- movementAnalysis_2D(segmentMovement$time,segmentMovement$ll_x,segmentMovement$ll_y)
            mastersheet$LL.xy_Distance[k] <- ll.xy$distance
            mastersheet$LL.xy_Displacement[k] <- ll.xy$displacement
            mastersheet$LL.xy_SpeedAvg[k] <- ll.xy$speed_avg
            mastersheet$LL.xy_SpeedMax[k] <- ll.xy$speed_max
            
            rm(ll.xy)
            
            # * XY (Decoupled)  -
            dll.xy <- movementAnalysis_2D(segmentMovement$time,segmentMovement$decoupled_ll_x,segmentMovement$decoupled_ll_y)
            mastersheet$dLL.xy_Distance[k] <- dll.xy$distance
            mastersheet$dLL.xy_Displacement[k] <- dll.xy$displacement
            mastersheet$dLL.xy_SpeedAvg[k] <- dll.xy$speed_avg
            mastersheet$dLL.xy_SpeedMax[k] <- dll.xy$speed_max
            
            rm(dll.xy)
            
            # * Convex Hull - LL 
            if (include_LL.cHull == TRUE) {mastersheet$LL.cHull[k] <- cHull(segmentMovement$ll_x,segmentMovement$ll_y)}
            if (include_dLL.cHull == TRUE) {mastersheet$dLL.cHull[k] <- cHull(segmentMovement$decoupled_ll_x,segmentMovement$decoupled_ll_y)}
            
            # * Lip Aperture 
            segmentMovement$lipAperture <- segmentMovement$ul_y - segmentMovement$ll_y
            if (include_lipAperture.min) {mastersheet$lipAperture.min[k] <- base::min(segmentMovement$lipAperture, na.rm = TRUE)}
          }
          # Jaw Measures 
          {
            # * X 
            jaw.x <- movementAnalysis_1D(segmentMovement$time,segmentMovement$jaw_x)
            mastersheet$Jaw.x_Onset[k] <- jaw.x$onset
            mastersheet$Jaw.x_Offset[k] <- jaw.x$offset
            mastersheet$Jaw.x_Distance[k] <- jaw.x$distance
            mastersheet$Jaw.x_Displacement[k] <- jaw.x$displacement
            mastersheet$Jaw.x_SpeedAvg[k] <- jaw.x$speed_avg
            mastersheet$Jaw.x_SpeedMax[k] <- jaw.x$speed_max
            
            rm(jaw.x)
            
            # * Y 
            jaw.y <- movementAnalysis_1D(segmentMovement$time,segmentMovement$jaw_y)
            mastersheet$Jaw.y_Onset[k] <- jaw.y$onset
            mastersheet$Jaw.y_Offset[k] <- jaw.y$offset
            mastersheet$Jaw.y_Distance[k] <- jaw.y$distance
            mastersheet$Jaw.y_Displacement[k] <- jaw.y$displacement
            mastersheet$Jaw.y_SpeedAvg[k] <- jaw.y$speed_avg
            mastersheet$Jaw.y_SpeedMax[k] <- jaw.y$speed_max
            
            rm(jaw.y)
            
            # * XY  -
            jaw.xy <- movementAnalysis_2D(segmentMovement$time,segmentMovement$jaw_x,segmentMovement$jaw_y)
            mastersheet$Jaw.xy_Distance[k] <- jaw.xy$distance
            mastersheet$Jaw.xy_Displacement[k] <- jaw.xy$displacement
            mastersheet$Jaw.xy_SpeedAvg[k] <- jaw.xy$speed_avg
            mastersheet$Jaw.xy_SpeedMax[k] <- jaw.xy$speed_max
            
            rm(jaw.xy)
            
            # * Convex Hull - Jaw 
            if (include_Jaw.cHull == TRUE) {mastersheet$Jaw.cHull[k] <- cHull(segmentMovement$jaw_x,segmentMovement$jaw_y)}
          }
        }
        
        # Removing unwanted Movement Segment Measures
        # TF - X
        if (MS.include_TF_x.PathDist == FALSE) {segmentMovement %>% dplyr::select(!TF_x.PathDist) -> segmentMovement}
        if (MS.include_TF_x.Active == FALSE) {segmentMovement %>% dplyr::select(!TF_x.Active) -> segmentMovement}
        if (MS.include_TF_x.Steady == FALSE) {segmentMovement %>% dplyr::select(!TF_x.Steady) -> segmentMovement}
        if (MS.include_dTF_x.PathDist == FALSE) {segmentMovement %>% dplyr::select(!dTF_x.PathDist) -> segmentMovement}
        if (MS.include_dTF_x.Active == FALSE) {segmentMovement %>% dplyr::select(!dTF_x.Active) -> segmentMovement}
        if (MS.include_dTF_x.Steady == FALSE) {segmentMovement %>% dplyr::select(!dTF_x.Steady) -> segmentMovement}
        # TF - Y
        if (MS.include_TF_y.PathDist == FALSE) {segmentMovement %>% dplyr::select(!TF_y.PathDist) -> segmentMovement}
        if (MS.include_TF_y.Active == FALSE) {segmentMovement %>% dplyr::select(!TF_y.Active) -> segmentMovement}
        if (MS.include_TF_y.Steady == FALSE) {segmentMovement %>% dplyr::select(!TF_y.Steady) -> segmentMovement}
        if (MS.include_dTF_y.PathDist == FALSE) {segmentMovement %>% dplyr::select(!dTF_y.PathDist) -> segmentMovement}
        if (MS.include_dTF_y.Active == FALSE) {segmentMovement %>% dplyr::select(!dTF_y.Active) -> segmentMovement}
        if (MS.include_dTF_y.Steady == FALSE) {segmentMovement %>% dplyr::select(!dTF_y.Steady) -> segmentMovement}
        # TB - X
        if (MS.include_TB_x.PathDist == FALSE) {segmentMovement %>% dplyr::select(!TB_x.PathDist) -> segmentMovement}
        if (MS.include_TB_x.Active == FALSE) {segmentMovement %>% dplyr::select(!TB_x.Active) -> segmentMovement}
        if (MS.include_TB_x.Steady == FALSE) {segmentMovement %>% dplyr::select(!TB_x.Steady) -> segmentMovement}
        if (MS.include_dTB_x.PathDist == FALSE) {segmentMovement %>% dplyr::select(!dTB_x.PathDist) -> segmentMovement}
        if (MS.include_dTB_x.Active == FALSE) {segmentMovement %>% dplyr::select(!dTB_x.Active) -> segmentMovement}
        if (MS.include_dTB_x.Steady == FALSE) {segmentMovement %>% dplyr::select(!dTB_x.Steady) -> segmentMovement}
        # TB - Y
        if (MS.include_TB_y.PathDist == FALSE) {segmentMovement %>% dplyr::select(!TB_y.PathDist) -> segmentMovement}
        if (MS.include_TB_y.Active == FALSE) {segmentMovement %>% dplyr::select(!TB_y.Active) -> segmentMovement}
        if (MS.include_TB_y.Steady == FALSE) {segmentMovement %>% dplyr::select(!TB_y.Steady) -> segmentMovement}
        if (MS.include_dTB_y.PathDist == FALSE) {segmentMovement %>% dplyr::select(!dTB_y.PathDist) -> segmentMovement}
        if (MS.include_dTB_y.Active == FALSE) {segmentMovement %>% dplyr::select(!dTB_y.Active) -> segmentMovement}
        if (MS.include_dTB_y.Steady == FALSE) {segmentMovement %>% dplyr::select(!dTB_y.Steady) -> segmentMovement}
        
        if (MS.include_LipAperture == FALSE) {segmentMovement %>% dplyr::select(!lipAperture) -> segmentMovement}

        # Movement Segments
        sheetname <- as.character(targetIntervals$targetNames)[k]
        openxlsx::addWorksheet(movementSegments, sheetname)
        openxlsx::writeData(movementSegments, sheet = sheetname, x = segmentMovement)
        
        # Remove unwanted measures from Mastersheet
        # TF measures 
        # TF - X
        if (include_TF.x_Onset == FALSE) {mastersheet %>% dplyr::select(!TF.x_Onset) -> mastersheet}
        if (include_TF.x_Offset == FALSE) {mastersheet %>% dplyr::select(!TF.x_Offset) -> mastersheet}
        if (include_TF.x_Distance == FALSE) {mastersheet %>% dplyr::select(!TF.x_Distance) -> mastersheet}
        if (include_TF.x_Displacement == FALSE) {mastersheet %>% dplyr::select(!TF.x_Displacement) -> mastersheet}
        if (include_TF.x_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!TF.x_SpeedAvg) -> mastersheet}
        if (include_TF.x_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!TF.x_SpeedMax) -> mastersheet}
        if (include_TF.x_steadyDuration == FALSE) {mastersheet %>% dplyr::select(!TF.x_steadyDuration) -> mastersheet}
        # TF - Y
        if (include_TF.y_Onset == FALSE) {mastersheet %>% dplyr::select(!TF.y_Onset) -> mastersheet}
        if (include_TF.y_Offset == FALSE) {mastersheet %>% dplyr::select(!TF.y_Offset) -> mastersheet}
        if (include_TF.y_Distance == FALSE) {mastersheet %>% dplyr::select(!TF.y_Distance) -> mastersheet}
        if (include_TF.y_Displacement == FALSE) {mastersheet %>% dplyr::select(!TF.y_Displacement) -> mastersheet}
        if (include_TF.y_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!TF.y_SpeedAvg) -> mastersheet}
        if (include_TF.y_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!TF.y_SpeedMax) -> mastersheet}
        if (include_TF.y_steadyDuration == FALSE) {mastersheet %>% dplyr::select(!TF.y_steadyDuration) -> mastersheet}
        # TF - XY
        if (include_TF.xy_Distance == FALSE) {mastersheet %>% dplyr::select(!TF.xy_Distance) -> mastersheet}
        if (include_TF.xy_Displacement == FALSE) {mastersheet %>% dplyr::select(!TF.xy_Displacement) -> mastersheet}
        if (include_TF.xy_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!TF.xy_SpeedAvg) -> mastersheet}
        if (include_TF.xy_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!TF.xy_SpeedMax) -> mastersheet}
        # Decoupled TF measures
        # dTF - X
        if (include_dTF.x_Onset == FALSE) {mastersheet %>% dplyr::select(!dTF.x_Onset) -> mastersheet}
        if (include_dTF.x_Offset == FALSE) {mastersheet %>% dplyr::select(!dTF.x_Offset) -> mastersheet}
        if (include_dTF.x_Distance == FALSE) {mastersheet %>% dplyr::select(!dTF.x_Distance) -> mastersheet}
        if (include_dTF.x_Displacement == FALSE) {mastersheet %>% dplyr::select(!dTF.x_Displacement) -> mastersheet}
        if (include_dTF.x_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!dTF.x_SpeedAvg) -> mastersheet}
        if (include_dTF.x_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!dTF.x_SpeedMax) -> mastersheet}
        if (include_dTF.x_steadyDuration == FALSE) {mastersheet %>% dplyr::select(!dTF.x_steadyDuration) -> mastersheet}
        # dTF - Y
        if (include_dTF.y_Onset == FALSE) {mastersheet %>% dplyr::select(!dTF.y_Onset) -> mastersheet}
        if (include_dTF.y_Offset == FALSE) {mastersheet %>% dplyr::select(!dTF.y_Offset) -> mastersheet}
        if (include_dTF.y_Distance == FALSE) {mastersheet %>% dplyr::select(!dTF.y_Distance) -> mastersheet}
        if (include_dTF.y_Displacement == FALSE) {mastersheet %>% dplyr::select(!dTF.y_Displacement) -> mastersheet}
        if (include_dTF.y_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!dTF.y_SpeedAvg) -> mastersheet}
        if (include_dTF.y_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!dTF.y_SpeedMax) -> mastersheet}
        if (include_dTF.y_steadyDuration == FALSE) {mastersheet %>% dplyr::select(!dTF.y_steadyDuration) -> mastersheet}
        # dTF - XY
        if (include_dTF.xy_Distance == FALSE) {mastersheet %>% dplyr::select(!dTF.xy_Distance) -> mastersheet}
        if (include_dTF.xy_Displacement == FALSE) {mastersheet %>% dplyr::select(!dTF.xy_Displacement) -> mastersheet}
        if (include_dTF.xy_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!dTF.xy_SpeedAvg) -> mastersheet}
        if (include_dTF.xy_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!dTF.xy_SpeedMax) -> mastersheet}
        # TB measures 
        # TB - X
        if (include_TB.x_Onset == FALSE) {mastersheet %>% dplyr::select(!TB.x_Onset) -> mastersheet}
        if (include_TB.x_Offset == FALSE) {mastersheet %>% dplyr::select(!TB.x_Offset) -> mastersheet}
        if (include_TB.x_Distance == FALSE) {mastersheet %>% dplyr::select(!TB.x_Distance) -> mastersheet}
        if (include_TB.x_Displacement == FALSE) {mastersheet %>% dplyr::select(!TB.x_Displacement) -> mastersheet}
        if (include_TB.x_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!TB.x_SpeedAvg) -> mastersheet}
        if (include_TB.x_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!TB.x_SpeedMax) -> mastersheet}
        if (include_TB.x_steadyDuration == FALSE) {mastersheet %>% dplyr::select(!TB.x_steadyDuration) -> mastersheet}
        # TB - Y
        if (include_TB.y_Onset == FALSE) {mastersheet %>% dplyr::select(!TB.y_Onset) -> mastersheet}
        if (include_TB.y_Offset == FALSE) {mastersheet %>% dplyr::select(!TB.y_Offset) -> mastersheet}
        if (include_TB.y_Distance == FALSE) {mastersheet %>% dplyr::select(!TB.y_Distance) -> mastersheet}
        if (include_TB.y_Displacement == FALSE) {mastersheet %>% dplyr::select(!TB.y_Displacement) -> mastersheet}
        if (include_TB.y_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!TB.y_SpeedAvg) -> mastersheet}
        if (include_TB.y_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!TB.y_SpeedMax) -> mastersheet}
        if (include_TB.y_steadyDuration == FALSE) {mastersheet %>% dplyr::select(!TB.y_steadyDuration) -> mastersheet}
        # TB - XY
        if (include_TB.xy_Distance == FALSE) {mastersheet %>% dplyr::select(!TB.xy_Distance) -> mastersheet}
        if (include_TB.xy_Displacement == FALSE) {mastersheet %>% dplyr::select(!TB.xy_Displacement) -> mastersheet}
        if (include_TB.xy_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!TB.xy_SpeedAvg) -> mastersheet}
        if (include_TB.xy_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!TB.xy_SpeedMax) -> mastersheet}
        # Decoupled TB measures
        # dTB - X
        if (include_dTB.x_Onset == FALSE) {mastersheet %>% dplyr::select(!dTB.x_Onset) -> mastersheet}
        if (include_dTB.x_Offset == FALSE) {mastersheet %>% dplyr::select(!dTB.x_Offset) -> mastersheet}
        if (include_dTB.x_Distance == FALSE) {mastersheet %>% dplyr::select(!dTB.x_Distance) -> mastersheet}
        if (include_dTB.x_Displacement == FALSE) {mastersheet %>% dplyr::select(!dTB.x_Displacement) -> mastersheet}
        if (include_dTB.x_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!dTB.x_SpeedAvg) -> mastersheet}
        if (include_dTB.x_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!dTB.x_SpeedMax) -> mastersheet}
        if (include_dTB.x_steadyDuration == FALSE) {mastersheet %>% dplyr::select(!dTB.x_steadyDuration) -> mastersheet}
        # dTB - Y
        if (include_dTB.y_Onset == FALSE) {mastersheet %>% dplyr::select(!dTB.y_Onset) -> mastersheet}
        if (include_dTB.y_Offset == FALSE) {mastersheet %>% dplyr::select(!dTB.y_Offset) -> mastersheet}
        if (include_dTB.y_Distance == FALSE) {mastersheet %>% dplyr::select(!dTB.y_Distance) -> mastersheet}
        if (include_dTB.y_Displacement == FALSE) {mastersheet %>% dplyr::select(!dTB.y_Displacement) -> mastersheet}
        if (include_dTB.y_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!dTB.y_SpeedAvg) -> mastersheet}
        if (include_dTB.y_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!dTB.y_SpeedMax) -> mastersheet}
        if (include_dTB.y_steadyDuration == FALSE) {mastersheet %>% dplyr::select(!dTB.y_steadyDuration) -> mastersheet}
        # dTB - XY
        if (include_dTB.xy_Distance == FALSE) {mastersheet %>% dplyr::select(!dTB.xy_Distance) -> mastersheet}
        if (include_dTB.xy_Displacement == FALSE) {mastersheet %>% dplyr::select(!dTB.xy_Displacement) -> mastersheet}
        if (include_dTB.xy_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!dTB.xy_SpeedAvg) -> mastersheet}
        if (include_dTB.xy_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!dTB.xy_SpeedMax) -> mastersheet}
        # UL measures 
        # UL - X
        if (include_UL.x_Onset == FALSE) {mastersheet %>% dplyr::select(!UL.x_Onset) -> mastersheet}
        if (include_UL.x_Offset == FALSE) {mastersheet %>% dplyr::select(!UL.x_Offset) -> mastersheet}
        if (include_UL.x_Distance == FALSE) {mastersheet %>% dplyr::select(!UL.x_Distance) -> mastersheet}
        if (include_UL.x_Displacement == FALSE) {mastersheet %>% dplyr::select(!UL.x_Displacement) -> mastersheet}
        if (include_UL.x_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!UL.x_SpeedAvg) -> mastersheet}
        if (include_UL.x_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!UL.x_SpeedMax) -> mastersheet}
        # UL - Y
        if (include_UL.y_Onset == FALSE) {mastersheet %>% dplyr::select(!UL.y_Onset) -> mastersheet}
        if (include_UL.y_Offset == FALSE) {mastersheet %>% dplyr::select(!UL.y_Offset) -> mastersheet}
        if (include_UL.y_Distance == FALSE) {mastersheet %>% dplyr::select(!UL.y_Distance) -> mastersheet}
        if (include_UL.y_Displacement == FALSE) {mastersheet %>% dplyr::select(!UL.y_Displacement) -> mastersheet}
        if (include_UL.y_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!UL.y_SpeedAvg) -> mastersheet}
        if (include_UL.y_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!UL.y_SpeedMax) -> mastersheet}
        # UL - XY
        if (include_UL.xy_Distance == FALSE) {mastersheet %>% dplyr::select(!UL.xy_Distance) -> mastersheet}
        if (include_UL.xy_Displacement == FALSE) {mastersheet %>% dplyr::select(!UL.xy_Displacement) -> mastersheet}
        if (include_UL.xy_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!UL.xy_SpeedAvg) -> mastersheet}
        if (include_UL.xy_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!UL.xy_SpeedMax) -> mastersheet}
        # LL measures 
        # LL - X
        if (include_LL.x_Onset == FALSE) {mastersheet %>% dplyr::select(!LL.x_Onset) -> mastersheet}
        if (include_LL.x_Offset == FALSE) {mastersheet %>% dplyr::select(!LL.x_Offset) -> mastersheet}
        if (include_LL.x_Distance == FALSE) {mastersheet %>% dplyr::select(!LL.x_Distance) -> mastersheet}
        if (include_LL.x_Displacement == FALSE) {mastersheet %>% dplyr::select(!LL.x_Displacement) -> mastersheet}
        if (include_LL.x_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!LL.x_SpeedAvg) -> mastersheet}
        if (include_LL.x_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!LL.x_SpeedMax) -> mastersheet}
        # LL - Y
        if (include_LL.y_Onset == FALSE) {mastersheet %>% dplyr::select(!LL.y_Onset) -> mastersheet}
        if (include_LL.y_Offset == FALSE) {mastersheet %>% dplyr::select(!LL.y_Offset) -> mastersheet}
        if (include_LL.y_Distance == FALSE) {mastersheet %>% dplyr::select(!LL.y_Distance) -> mastersheet}
        if (include_LL.y_Displacement == FALSE) {mastersheet %>% dplyr::select(!LL.y_Displacement) -> mastersheet}
        if (include_LL.y_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!LL.y_SpeedAvg) -> mastersheet}
        if (include_LL.y_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!LL.y_SpeedMax) -> mastersheet}
        # LL - XY
        if (include_LL.xy_Distance == FALSE) {mastersheet %>% dplyr::select(!LL.xy_Distance) -> mastersheet}
        if (include_LL.xy_Displacement == FALSE) {mastersheet %>% dplyr::select(!LL.xy_Displacement) -> mastersheet}
        if (include_LL.xy_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!LL.xy_SpeedAvg) -> mastersheet}
        if (include_LL.xy_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!LL.xy_SpeedMax) -> mastersheet}
        # Decoupled LL measures
        # dLL - X
        if (include_dLL.x_Onset == FALSE) {mastersheet %>% dplyr::select(!dLL.x_Onset) -> mastersheet}
        if (include_dLL.x_Offset == FALSE) {mastersheet %>% dplyr::select(!dLL.x_Offset) -> mastersheet}
        if (include_dLL.x_Distance == FALSE) {mastersheet %>% dplyr::select(!dLL.x_Distance) -> mastersheet}
        if (include_dLL.x_Displacement == FALSE) {mastersheet %>% dplyr::select(!dLL.x_Displacement) -> mastersheet}
        if (include_dLL.x_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!dLL.x_SpeedAvg) -> mastersheet}
        if (include_dLL.x_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!dLL.x_SpeedMax) -> mastersheet}
        # dLL - Y
        if (include_dLL.y_Onset == FALSE) {mastersheet %>% dplyr::select(!dLL.y_Onset) -> mastersheet}
        if (include_dLL.y_Offset == FALSE) {mastersheet %>% dplyr::select(!dLL.y_Offset) -> mastersheet}
        if (include_dLL.y_Distance == FALSE) {mastersheet %>% dplyr::select(!dLL.y_Distance) -> mastersheet}
        if (include_dLL.y_Displacement == FALSE) {mastersheet %>% dplyr::select(!dLL.y_Displacement) -> mastersheet}
        if (include_dLL.y_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!dLL.y_SpeedAvg) -> mastersheet}
        if (include_dLL.y_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!dLL.y_SpeedMax) -> mastersheet}
        # dLL - XY
        if (include_dLL.xy_Distance == FALSE) {mastersheet %>% dplyr::select(!dLL.xy_Distance) -> mastersheet}
        if (include_dLL.xy_Displacement == FALSE) {mastersheet %>% dplyr::select(!dLL.xy_Displacement) -> mastersheet}
        if (include_dLL.xy_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!dLL.xy_SpeedAvg) -> mastersheet}
        if (include_dLL.xy_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!dLL.xy_SpeedMax) -> mastersheet}
        # Jaw measures 
        # Jaw - X
        if (include_Jaw.x_Onset == FALSE) {mastersheet %>% dplyr::select(!Jaw.x_Onset) -> mastersheet}
        if (include_Jaw.x_Offset == FALSE) {mastersheet %>% dplyr::select(!Jaw.x_Offset) -> mastersheet}
        if (include_Jaw.x_Distance == FALSE) {mastersheet %>% dplyr::select(!Jaw.x_Distance) -> mastersheet}
        if (include_Jaw.x_Displacement == FALSE) {mastersheet %>% dplyr::select(!Jaw.x_Displacement) -> mastersheet}
        if (include_Jaw.x_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!Jaw.x_SpeedAvg) -> mastersheet}
        if (include_Jaw.x_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!Jaw.x_SpeedMax) -> mastersheet}
        # Jaw - Y
        if (include_Jaw.y_Onset == FALSE) {mastersheet %>% dplyr::select(!Jaw.y_Onset) -> mastersheet}
        if (include_Jaw.y_Offset == FALSE) {mastersheet %>% dplyr::select(!Jaw.y_Offset) -> mastersheet}
        if (include_Jaw.y_Distance == FALSE) {mastersheet %>% dplyr::select(!Jaw.y_Distance) -> mastersheet}
        if (include_Jaw.y_Displacement == FALSE) {mastersheet %>% dplyr::select(!Jaw.y_Displacement) -> mastersheet}
        if (include_Jaw.y_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!Jaw.y_SpeedAvg) -> mastersheet}
        if (include_Jaw.y_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!Jaw.y_SpeedMax) -> mastersheet}
        # Jaw - XY
        if (include_Jaw.xy_Distance == FALSE) {mastersheet %>% dplyr::select(!Jaw.xy_Distance) -> mastersheet}
        if (include_Jaw.xy_Displacement == FALSE) {mastersheet %>% dplyr::select(!Jaw.xy_Displacement) -> mastersheet}
        if (include_Jaw.xy_SpeedAvg == FALSE) {mastersheet %>% dplyr::select(!Jaw.xy_SpeedAvg) -> mastersheet}
        if (include_Jaw.xy_SpeedMax == FALSE) {mastersheet %>% dplyr::select(!Jaw.xy_SpeedMax) -> mastersheet}
        
        # Finishing the loop 
        k <- k + 1 # This is the counter for the analysis targets.
        mastersheet <- add_row(mastersheet)
        rm(segmentAcoustic,segmentF0,segmentMovement)
      }
      
      # Removing unwanted variables
      if (include_Duration_s == FALSE) {mastersheet %>% dplyr::select(!duration_s) -> mastersheet}
      if (include_Duration_ms == FALSE) {mastersheet %>% dplyr::select(!duration_ms) -> mastersheet}
      if (include_Onset_ms == FALSE) {mastersheet %>% dplyr::select(!onset_ms) -> mastersheet}
      if (include_Offset_ms == FALSE) {mastersheet %>% dplyr::select(!offset_ms) -> mastersheet}
      if (include_Onset_s == FALSE) {mastersheet %>% dplyr::select(!onset_s) -> mastersheet}
      if (include_Offset_s == FALSE) {mastersheet %>% dplyr::select(!offset_s) -> mastersheet}
      mastersheet %>% dplyr::select(!c(midValue_ms,midValue_s)) -> mastersheet
      if (include_TempMid_ms.Aco == FALSE) {mastersheet %>% dplyr::select(!tempMid_ms.Aco) -> mastersheet}
      if (include_TempMid_ms.Kin == FALSE) {mastersheet %>% dplyr::select(!tempMid_ms.Kin) -> mastersheet}
      
      # Saving the Data
      if (saveData == TRUE) {
      dir.create(paste("KASA"), showWarnings = FALSE)
      dir.create(paste("KASA/Analyzed Segments"), showWarnings = FALSE)
      dir.create(paste("KASA/Acoustic Segments"), showWarnings = FALSE)
      dir.create(paste("KASA/Movement Segments"), showWarnings = FALSE)
      
      # * Acoustic Segments
      AcoOutputFile <- paste("KASA/Acoustic Segments/",Study.ID,"_Acoustic Segments.xlsx",sep="")
      openxlsx::saveWorkbook(acousticSegments, AcoOutputFile, overwrite = TRUE)
      
      # * Movement Segments
      MovementOutputFile <- paste("KASA/Movement Segments/",Study.ID,"_Movement Segments.xlsx",sep="")
      openxlsx::saveWorkbook(movementSegments, MovementOutputFile, overwrite = TRUE)
      
      # * Data Analysis
      OutputFile <- paste("KASA/","Analyzed Segments/",Study.ID,"_Analyzed Segments.xlsx", sep="")
      openxlsx::write.xlsx(mastersheet, file = OutputFile,
                 sheetName = paste(Study.ID),
                 append = TRUE)
      }
      print(paste("Kinematic and Acoustic Analysis completed for ",Study.ID," (",Database.ID,").", sep = ""))
      list <- list("AnalyzedData" = mastersheet,
                   "AcousticSegments" = acousticSegments,
                   "MovementSegments" = movementSegments
      )
      return(list)
  }
