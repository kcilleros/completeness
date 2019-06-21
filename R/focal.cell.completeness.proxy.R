#' Proxy of completeness
#'
#' This function compares the richness of one focal cell to the richness estimated from the neighboring cells with one the following indices: Chao, jackknife (1 or 2), bootstrap.
#' @param focal.cell rownames of the focal cell to estimate completeness.
#' @param grid.with.species a data.frame with long, lat, id.cell and species occurence data.
#' @param radius will define the number of neighbor cells. Default to 10000 (1 cell for 10x10km grid).
#' @param focal.cell.in.estimate does the focal need to be taken into account to estimate diversity. default to TRUE
#' @param grid.with.species a data.frame with long, lat, id.cell and species occurence data.
#' @param X.long character, the name of the column that contains X or longitude.
#' @param Y.lat character, the name of the column that contains Y or latitude.
#' @param index character, one of the index available in specpool: chao, jack1, jack2, boot.
#' @keywords completeness
#' @export
#' @examples
#' focal.cell.completeness.proxy()

focal.cell.completeness.proxy <- function(focal.cell, 
                                    grid.with.species, 
                                    radius = 10000, 
                                    focal.cell.in.estimate = T,
                                    X.long, Y.lat, 
                                    index){
  
  #This function compares the richness of one focal cell to the richness estimated from the neighboring cells with the Chao's index 
  
  #Function arguments:
  ##focal.cell: rownames of the focal cell to estimate completeness
  ##grid.with.species: a data.frame with long, lat, id.cell and species occurence data
  ##radius: will define the number of neighbor cells
  ##focal.cell.in.estimate: does the focal need to be taken into account to estimate diversity
  ##X.long: the name of the column that contains X or longitude (character)
  ##Y.lat: the name of the column that contains Y or latitude (character)
  ##index: one of the index available in specpool: chao, jack1, jack2, boot (character)
  
  require(vegan)
  require(CommEcol)
  require(magrittr)
  
  #Select the neighboring cells based on the radius
  resu<-select.window(xf=grid.with.species[focal.cell,X.long],
                      yf=grid.with.species[focal.cell,Y.lat], 
                      radius=radius, 
                      xydata=grid.with.species)
  if(nrow(resu) < 2){ #if cell isolated
    return(NA)
  }else{
    
    #Computation of Chao estimate for occurence data
    if(index == "chao"){
      if(focal.cell.in.estimate == T){
        local.estimate <- specpool(resu[,-c(1:3)])$chao
      }else{
        local.estimate <- specpool(resu[-1,-c(1:3)])$chao
      }
    }
    
    #Computation of First order jackknife estimate for occurence data
    if(index == "jack1"){
      if(focal.cell.in.estimate == T){
        local.estimate <- specpool(resu[,-c(1:3)])$jack1
      }else{
        local.estimate <- specpool(resu[-1,-c(1:3)])$jack1
      }
    }	
    
    #Computation of Second order jackknife estimate for occurence data
    if(index == "jack2"){
      if(focal.cell.in.estimate == T){
        local.estimate <- specpool(resu[,-c(1:3)])$jack2
      }else{
        local.estimate <- specpool(resu[-1,-c(1:3)])$jack2
      }
    }	   
    
    #Computation of Bootstrap estimate for occurence data
    if(index == "boot"){
      if(focal.cell.in.estimate == T){
        local.estimate <- specpool(resu[,-c(1:3)])$boot
      }else{
        local.estimate <- specpool(resu[-1,-c(1:3)])$boot
      }
    }	
    
    #Species richness of the focal cell
    focal.cell.richness <- grid.with.species[focal.cell,-c(1:3)] %>% sum
    
    #percentage of the cell richness over Chao estimate
    return(100*focal.cell.richness/local.estimate)
  }
  
  
  
}
