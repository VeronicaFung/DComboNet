#' Generate drugseed-score file
#' @description The function \code{seedscore} is to generate the seed-score file
#' @param seeds strings, randow walk with restart process will start from the
#'   provided drug seed which should be included in drug network
#' @param eta numeric parameter to controls the probability of restarting in the
#'   corresponding network
#' @return \code{seeds_score} drugseed-score \code{data.frame}
#' @export
#' @examples
#' seeds = "Sorafenib"
#' seedscore(seeds = seeds)

seedscore <- function(seeds,
                      eta = 1){

  drugs = seeds[1]
  genes = seeds[-1]
  if(length(drugs)!= 0 && length(genes) ==0){
    seeds_score <- data.frame(Seeds_ID = seeds,
                              Score = eta/length(drugs))
  }else if(length(drugs)!= 0 && length(genes) !=0){
    drugseed_score = eta/length(drugs)
    geneseed_score = rep((1-eta)/length(genes),length(genes))
    seeds_score <- data.frame(Seeds_ID = seeds,
                              Score = c(drugseed_score,geneseed_score),
                              stringsAsFactors = FALSE)
  }
  return(seeds_score)
}


