#' Richness estimate for a cell based on cells from a smaller grid within the one of interest
#'
#' This function computes richness estimates for a focal cell and returns the percentage of observed richness compared to one estimate index. To estimate richness, sampling units used are cells from a smaller grid that are included in the focal cell.
#' @param focal.cell id of the focal cell for which richness estimates are needed.
#' @param data.long.format.grid.subgrid.species a data.frame in long format (1 row = 1 occurence) with species, cell id for a grid (e.g. 50x50km) and cell id for a smaller grid (e.g. 10x10km).
#' @param indice.estimate character, one of the index available in specpool: chao (default), jack1, jack2, boot.
#' @param colname_cell character, the name of the column that contains the cell id of the largest grid.
#' @param colname_subcell character, the name of the column that contains the cell id of the small grid that will be used as sampling unit.
#' @param colname_species character, the name of the column that contains the species.
#' @keywords completeness
#' @export
#' @examples
#' estimate.cell.from.subcells()
#' 
#' 

estimate.cell.from.subcells <- function(focal.cell,
                                        data.long.format.grid.subgrid.species,
                                        indice.estimate = "chao",
                                        colname_cell,
                                        colname_subcell,
                                        colname_species){
  
  require(vegan)
  require(reshape)
  require(reshape2)
  require(magrittr)
  
  
  
  data_grid_subset <- subset(data.long.format.grid.subgrid.species, 
                             data.long.format.grid.subgrid.species[,colname_cell] == focal.cell)
  data_grid_subset %>%
    cast(data = .,
         formula = as.formula(paste0(colname_subcell, " ~ ", colname_species)), 
         fun.aggregate = length, value = colname_cell) %>%
    decostand(x = .[,-1], method = "pa") -> matrice_comm
  
  estimates <- specpool(matrice_comm)
  estimate <- 100*estimates$Species/estimates[,indice.estimate]
  
  return(c(cell = focal.cell,
         n_subcell = estimates$n,
         richness_estimate = round(estimate, digits = 2)))
  
}
