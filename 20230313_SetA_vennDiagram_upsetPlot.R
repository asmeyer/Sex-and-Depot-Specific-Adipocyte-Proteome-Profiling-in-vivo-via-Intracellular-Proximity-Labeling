# 2023-03-13
# TLS
# Cytoplasmic 

# usage: venn diagram and upset plot for LC-MS/MS for enrichment score data for 
#INT versus POS

# last updated: 2023-03-14


################################################################## begin ####################################################################

######## venn diagram #######

library(ggvenn)

# variables
POS_M <-c ("P20152","Q6NSR3","Q9QZS7","E9Q616","Q9Z0F7","P97447")
POS_F <-c ("Q3V037","Q9QYR9")
INT_M <-c ("Q9CPZ6","E9Q616","Q9QZS7")
INT_F <-c ("Q9QYR9")


allT <- list ("SAT_Male"=POS_M, "SAT_Female"=POS_F, "BAT_Male"=INT_M, "BAT_Female"=INT_F)

allT  

# venn all tissues 
ggvenn(allT, fill_color = c("#257DCF", "#42853D", "#AA0A3C", "grey30"), show_percentage = F)



####### upset plot #######

library(UpSetR)
protein_list <- list(SAT_Male = c("P20152","Q6NSR3","Q9QZS7","E9Q616","Q9Z0F7","P97447"), SAT_Female =c("Q3V037","Q9QYR9"), BAT_Male = c("Q9CPZ6","E9Q616","Q9QZS7"), 
                  BAT_Female = c("Q9QYR9"))

protein_list

plot1 <- upset(fromList(protein_list), order.by = "freq")
str(plot1)
plot1


