#' Random walk with restart algorithm
#' @description The function \code{rwr} is to calculate global proximity between
#'   drug seed and other drugs in complex network.
#' @details The equation of random walk with restart is as follow:
#'   \deqn{P_{t+1} = (1-\sigma)H * P_{t} + \sigma P_{0}}
#' @references Valdeolivas A, Tichit L, Navarro C, et al. Random walk with
#'   restart on multiplex and heterogeneous biological networks[J].
#'   Bioinformatics, 2019, 35(3): 497-505.
#' @param tm transition matrix generated from drug-gene/pathway interaction
#'   networks via \code{TransitionMatrix}
#' @param r numeric, global restart parameter
#' @param seeds_score drugseed-score \code{data.frame}, generated via
#'   \code{seedscore}
#' @return  global proximity between drug seed and other drug in
#'   drug-gene/pathway interaction network.
#' @export
#'
#' @examples
#' \dontrun{
#' tm = TransitionMatrix(...)
#' seeds = seedscore(seeds = "Sorafenib")
#' rwr(tm, r = 0.7, seeds_score = seeds)
#' }
#'


rwr = function(tm,
               r = 0.7,
               seeds_score){
  residue <- 1
  iter <- 1
  P0 <- matrix(0,nrow = ncol(tm),ncol=1)
  for(i in 1:nrow(seeds_score)){
    P0[which(rownames(tm) %in% seeds_score[i,1])] <- seeds_score[i,2]}
  P0  <- P0/sum(P0)
  Pt <-  as.matrix(P0)

  while(residue >= 1e-10){
    Pt1 <- (1-r) * (tm %*% Pt) + r * P0

    residue = sum(abs(Pt1 - Pt))
    Pt <- Pt1
    iter <- iter + 1
    # print(paste(iter,residue))
  }
  return(Pt1)
}
