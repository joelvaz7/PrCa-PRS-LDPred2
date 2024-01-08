library(dplyr)


ethnicity_col <- fread("Analysis/UKBB/UKBB_PrCa_7296_ETHcol.txt")

ethnicity_col <- left_join(y, ethnicity_col, by = 'FID')

top_level_ethnicities = list("White" = c("British", "Irish", "Any other white background", "White"),
                             "Mixed" = c("Mixed", "White and Black African", "White and Black Caribbean", "Any other mixed background", "White and Asian"),
                             "Asian or Asian British" = c("Indian", "Pakistani", "Bangladeshi", "Any other Asian background"),
                             "Black or Black British" = c("African", "Caribbean"),
                             "Chinese" = "Chinese",
                             "Do not know" = "Do not know",
                             "Prefer not to answer" = "Prefer not to answer",
                             "Other ethnic group"  = "Other ethnic group")

# Create a new column based on the mapping
ethnicity_col$TopLevelEthnicity <- sapply(ethnicity_col$Ethnicity, function(ethnicity) {
  for (top_level_ethnicity in names(top_level_ethnicities)) {
    if (ethnicity %in% top_level_ethnicities[[top_level_ethnicity]]) {
      return(top_level_ethnicity)
    }
  }
  return(NA) # If no match is found
})
ggplot(ethnicity_col, aes(x=PC1,y=PC2,color=TopLevelEthnicity)) + geom_point(alpha = 0.9)
