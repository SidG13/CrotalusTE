###    ###
# This script is used to make a simplified figure to show percentages of re-wiring, based on BWP's schematic #
###    ###

library(tidyverse)

trin_mat <- read.table('../z_data/2SPECIES_Joint_Tip100_TE_TFBS_Ciiider_scan_footprint_bound_trinMat_May2025.txt', header = T)
re_sites <- data.frame(
      stringsAsFactors = FALSE,
             TFBS.Name = c("NFATC1","Creb3l2",
                           "TBX3","Irx2","Irx2.1","EHF","ELF5","NFATC1.1",
                           "TBX3.3","MEIS1.2","TBX3.1","CREB3L1","Creb3l2.1",
                           "TBX3.2","MEIS1.1","Creb3l2.2","MEIS1","MEIS1.3",
                           "Dlx4","POU6F1","PITX2"),
  TFBS.Alignment.Start = c(92L,175L,177L,213L,
                           214L,230L,232L,233L,245L,302L,410L,454L,475L,
                           476L,480L,513L,517L,530L,651L,691L,702L),
    TFBS.Alignment.End = c(112L,190L,198L,229L,
                           230L,244L,242L,242L,308L,308L,424L,508L,493L,
                           508L,508L,527L,530L,537L,665L,707L,718L),
   TFBS.Homology.Group = c("L2166","L250","L261",
                           "L358","L359","L411","L415","L418","L488",
                           "L674","L1039","L1186","L1236","L1244","L1283",
                           "L1352","L1371","L1415","L1779","L1888","L1936")
) %>% 
  arrange(TFBS.Alignment.Start)

trin_mat_re <- trin_mat[, re_sites$TFBS.Homology.Group]

trin_mat_re2 <- trin_mat_re %>%
  rownames_to_column('Tip100') %>%
  mutate(type = if_else(grepl('SVSP', Tip100), 'SVSP', 'nonSVSP'))

tfbs_map <- setNames(re_sites$TFBS.Name, re_sites$TFBS.Homology.Group)

colnames(trin_mat_re2) <- colnames(trin_mat_re2) %>%
  map_chr(~ ifelse(.x %in% names(tfbs_map), tfbs_map[.x], .x))

trin_long <- trin_mat_re2 %>%
  pivot_longer(cols = -c(Tip100, type), names_to = "TFBS", values_to = "Value") %>%
  mutate(Value = factor(Value, levels = c(0,1,2)))

tfbs_order <- re_sites %>%
  arrange(TFBS.Alignment.Start) %>%
  pull(TFBS.Name)

trin_long <- trin_long %>%
  mutate(TFBS = factor(TFBS, levels = tfbs_order))

trin_frac <- trin_long %>%
  group_by(type, TFBS, Value) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(type, TFBS) %>%
  mutate(fraction = count / sum(count))

