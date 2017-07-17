library(tidyverse)

aln <- read.delim("aln.table2", sep="\t")
reg <- read.delim("mergedSpecies.regions", sep="\t", header=T)
go <- read.delim("human_go.txt", sep="\t", header=T)

reg = reg %>% select(-GO)

names(aln) <- c("Name", "desc",
                "GID", "PID1", "PID2", "PID3",
                "Matches1", "Matches2", "Matches3",
                "alignPosition", "percentGg", "percentDr"
                )
go2 = go %>%
    select(Ensembl.Protein.ID, GO.Term.Accession) %>%
    group_by(Ensembl.Protein.ID) %>%
    summarise(goAcc = toString(GO.Term.Accession)) %>%
    ungroup
names(go2) = c("PID1", "GOterms") 
df1 <- aln %>% full_join(reg)
mainTable <- df1 %>% left_join(go2)
save(mainTable, file="mainTable.Rdata")
