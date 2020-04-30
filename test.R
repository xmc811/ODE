
get_mtx_from_dds(rnaseq[[1]], rnaseq[[2]]) %>%
    mtx_rescale()

test <- rlog(rnaseq[[1]], blind = FALSE)


get_mtx_from_dds(rnaseq[[1]], rnaseq[[2]], raw = T)

get_mtx_from_dds(rnaseq[[1]], rnaseq[[2]])

deseq_cluster(rnaseq[[1]], rnaseq[[2]], dir = TRUE, palette = "Spectral")

a <- get_nm_count_dds(rnaseq[[1]], rnaseq[[2]], "Grade")

rnaseq[[2]]@elementMetadata[2,2]

deseq_box(rnaseq[[1]], rnaseq[[2]], "Grade", "Set2")
