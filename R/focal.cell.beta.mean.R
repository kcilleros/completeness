#' Mean beta diversity for focal cell
#'
#' This function computes the mean beta-diversity between the focal cell and its neighboring cells (defined according to the radius).
#' @param focal.cell rownames of the focal cell to estimate completeness.
#' @param grid.with.species a data.frame with long, lat, id.cell and species occurence data.
#' @param radius will define the number of neighbor cells. Default to 10000 (1 cell for 10x10km grid).
#' @param grid.with.species a data.frame with long, lat, id.cell and species occurence data.
#' @param X.long character, the name of the column that contains X or longitude.
#' @param Y.lat character, the name of the column that contains Y or latitude.
#' @param index.family family of dissimilarity indices, partial match of "sorensen" or "jaccard".
#' @param component which component (turnover, nestedness, total) of the beta-diversity is returned by the function: "sim", "sne" or "sor", or jtu", "jne" or "jac".
#' @keywords completeness
#' @export
#' @examples
#' focal.cell.beta.mean()

focal.cell.beta.mean <- function(focal.cell, 
                                 grid.with.species, 
                                 radius = 10000, 
                                 X.long, Y.lat, 
                                 index.family = "jac", 
                                 component = "dis"){
  
  #This function computes the mean beta-diversity between the focal cell and its neighboring cells (defined according to the radius)
  
  #Function arguments:
  ##focal.cell: rownames of the focal cell to estimate completeness
  ##grid.with.species: a data.frame with long, lat, id.cell and species occurence data
  ##radius: will define the number of neighbor cells
  ##X.long: the name of the column that contains X or longitude (character)
  ##Y.lat: the name of the column that contains Y or latitude (character)
  ##index.family: family of dissimilarity indices, partial match of "sorensen" or "jaccard"
  ##component: which component (turnover, nestedness, total) of the beta-diversity is returned by the function: "sim", "sne" or "sor", or jtu", "jne" or "jac" 
  
  require(vegan)
  require(CommEcol)
  require(magrittr)
  require(betapart) #beta.pair function
  
  #Select the neighboring cells based on the radius
  resu<-select.window(xf=grid.with.species[focal.cell,X.long],
                      yf=grid.with.species[focal.cell,Y.lat], 
                      radius=radius, 
                      xydata=grid.with.species)
  if(nrow(resu) < 2 | ncol(resu) < 4){ #if cell isolated
    return(NA)
  }else{
    if(sum(resu[1,-c(1:3)]) < 1){
      return(NA)
    }else{
      #Compute beta diversity with beta.pair function
      beta_dist <- beta.pair(resu[,-c(1:3)], index.family = index.family)
      
      #Extract the component
      beta_dist_comp <- beta_dist[[paste0("beta.", component)]]
      
      #Extract the value for the focal cell
      beta_dist_comp_focal <- beta_dist_comp[1:(nrow(resu)-1)]
      
      #mean value
      beta_dist_comp_focal_mean <- mean(beta_dist_comp_focal, na.rm=T)
      
      
      #percentage of the cell richness over Chao estimate
      return(beta_dist_comp_focal_mean)}
    
  }
  
  
  
}
