devtools::load_all('../packages/twasAnalysis')



resistance_df_to_nested_list <- function(){}


test_df <- list(knockouts='',
                vals=c(1,2,3),
                children=list(
                  list(knockouts=c('SNQ2'),
                       vals=c(1,2,3),
                       children=NULL
                  )
                )
  
)