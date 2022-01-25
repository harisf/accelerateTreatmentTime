#' Simulate multistate data
#'
#' @param n Sample size
#'
#' @export
simulateData <- function(n){

  filename <- "sim_data"

  set.seed(1)

  # Generate multistate data with 5 states ----------------------------------

  # Transition rates
  a13 <- 0.22; a14 <- 0.7; a15 <- 0.08 # sum(a13+a14+a15)
  a23 <- 0.45; a24 <- 0.45; a25 <- 0.1 # sum(a23+a24+a25)
  a34 <- 0.17; a35 <- 0.33 # sum(a34+a35)
  a45 <- 0.02

  # Time grid
  endTime <- 10
  xx <- seq(0,endTime,length.out = 20)

  # It is possible to specify time-varying transition rates. Here, they are
  # constant in time
  spl13 <- smooth.spline(xx,seq(a13,a13,length.out=length(xx)))
  spl14 <- smooth.spline(xx,seq(a14,a14,length.out=length(xx)))
  spl15 <- smooth.spline(xx,seq(a15,a15,length.out=length(xx)))
  spl23 <- smooth.spline(xx,seq(a23,a23,length.out=length(xx)))
  spl24 <- smooth.spline(xx,seq(a24,a24,length.out=length(xx)))
  spl25 <- smooth.spline(xx,seq(a25,a25,length.out=length(xx)))
  spl34 <- smooth.spline(xx,seq(a34,a34,length.out=length(xx)))
  spl35 <- smooth.spline(xx,seq(a35,a35,length.out=length(xx)))
  spl45 <- smooth.spline(xx,seq(a45,a45,length.out=length(xx)))

  # Transition matrix
  transitionMatrix <- vector("list",length = 5)
  transitionMatrix[[1]]$smoothspline[[1]] <- list(NULL)
  transitionMatrix[[1]]$smoothspline[[2]] <- list(NULL)
  transitionMatrix[[1]]$smoothspline[[3]] <- list(spl13)
  transitionMatrix[[1]]$smoothspline[[4]] <- list(spl14)
  transitionMatrix[[1]]$smoothspline[[5]] <- list(spl15)
  transitionMatrix[[1]]$neighbours <- c(3,4,5)
  transitionMatrix[[2]]$smoothspline[[1]] <- list(NULL)
  transitionMatrix[[2]]$smoothspline[[2]] <- list(NULL)
  transitionMatrix[[2]]$smoothspline[[3]] <- list(spl23)
  transitionMatrix[[2]]$smoothspline[[4]] <- list(spl24)
  transitionMatrix[[2]]$smoothspline[[5]] <- list(spl25)
  transitionMatrix[[2]]$neighbours <- c(3,4,5)
  transitionMatrix[[3]]$smoothspline[[1]] <- list(NULL)
  transitionMatrix[[3]]$smoothspline[[2]] <- list(NULL)
  transitionMatrix[[3]]$smoothspline[[3]] <- list(NULL)
  transitionMatrix[[3]]$smoothspline[[4]] <- list(spl34)
  transitionMatrix[[3]]$smoothspline[[5]] <- list(spl35)
  transitionMatrix[[3]]$neighbours <- c(4,5)
  transitionMatrix[[4]]$smoothspline[[1]] <- list(NULL)
  transitionMatrix[[4]]$smoothspline[[2]] <- list(NULL)
  transitionMatrix[[4]]$smoothspline[[3]] <- list(NULL)
  transitionMatrix[[4]]$smoothspline[[4]] <- list(NULL)
  transitionMatrix[[4]]$smoothspline[[5]] <- list(spl45)
  transitionMatrix[[4]]$neighbours <- c(5)
  transitionMatrix[[5]]$smoothspline[[1]] <- list(NULL)
  transitionMatrix[[5]]$smoothspline[[2]] <- list(NULL)
  transitionMatrix[[5]]$smoothspline[[3]] <- list(NULL)
  transitionMatrix[[5]]$smoothspline[[4]] <- list(NULL)
  transitionMatrix[[5]]$smoothspline[[5]] <- list(NULL)
  transitionMatrix[[5]]$neighbours <- NULL

  # Starting states: 1 and 2
  dfr <- miscFunctions::generateFrame(n, endTime, transitionMatrix, c(1,2))

  # Terminating state: 5
  dfr <- dfr[!(dfr$from.state == 5),]


  # Create labels for each state --------------------------------------------

  # state --> label
  # 1 --> diag
  # 2 --> diag
  # 3 --> conf.t
  # 4 --> treat
  # 5 --> dead

  dt <- as.data.table(dfr)

  # X: baseline confounder
  dt[, X := ifelse(from.state[[1]] == 1, 0, 1), by = id]

  # L: time-varying confounder
  dt[, L := ifelse(to.state == 3, 1, 0)]
  dt[, L := cumsum(L), by = id]
  # L(t) = I(T_L < t) is left-continous, while Lcurrent(t) = I(T_L <= t) is
  # right-continuous
  dt[, Lcurrent := L]
  dt[to.state == 3, L := 0]

  # A: time-varying treatment
  dt[, A := ifelse(to.state == 4, 1, 0)]
  dt[, A := cumsum(A), by = id]
  dt[to.state == 4, A := A - 1]

  # D: death
  dt[, D := 0]
  dt[to.state == 5, D := 1]

  dt[from.state %in% c(1, 2), from.state.str := "diag"]
  dt[from.state == 3, from.state.str := "conf.t"]
  dt[from.state == 4, from.state.str := "treat"]

  dt[to.state == 3, to.state.str := "conf.t"]
  dt[to.state == 4, to.state.str := "treat"]
  dt[to.state == 5, to.state.str := "dead"]

  # Censored for death
  dt[to == 10, to.state := 0]
  dt[to == 10, to.state.str := "cens"]

  # Create a factor variable of states which are transitioned to. This is only
  # relevant for checking state occupancy probabilities
  dt[, to.state.fac := factor(to.state, levels = c(0, 3, 4, 5),
                              labels = c("cens", "conf.t", "treat", "dead"))]

  # Delete and rename variable
  dt[, c("from.state", "to.state") := NULL]
  setnames(dt, "from.state.str", "from.state")
  setnames(dt, "to.state.str", "to.state")
  setcolorder(dt, c("id", "from", "to", "from.state", "to.state", "to.state.fac"))

  # Save simulated data
  saveRDS(dt, file = paste0(filename, ".rds"))


  # Check transitions and state occupancies ---------------------------------
  #
  # # Transitions
  # survival::survcheck(survival::Surv(from, to, to.state.fac)~1, dt, id = id)
  #
  # # State occupancy probabilities
  # sfit <- survival::survfit(survival::Surv(from, to, to.state.fac)~X, dt, id = id)
  #
  # # Solid lines are for individuals with X = 0, dashed lines are for those with X = 1
  # plot(sfit, col = rep(c(1, 3, 2), each = 2), lty = 1:2)
  # text(rep(5, 3), c(0.1, 0.6, 0.3), c("conf.t", "treat", "death"), col = c(1, 3, 2))

}
