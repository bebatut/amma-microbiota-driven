get_perc = function(v, l){
    s = sum(v, na.rm=T)
    return(s)
    #return(round(s*100/l, digits=2))
}

get_sign_padj = function(dge_res) {
    return(get_perc(dge_res$padj < 0.05, sum(!is.na(dge_res$log2FoldChange))))
}

get_pos_sign_padj = function(dge_res) {
    return(get_perc(dge_res$padj < 0.05 & dge_res$log2FoldChange > 0, sum(!is.na(dge_res$log2FoldChange))))
}

get_neg_sign_padj = function(dge_res) {
    return(get_perc(dge_res$padj < 0.05 & dge_res$log2FoldChange < 0, sum(!is.na(dge_res$log2FoldChange))))
}

get_stats_padj  = function(dge_res) {
    v = c(get_sign_padj(dge_res), get_pos_sign_padj(dge_res), get_neg_sign_padj(dge_res))
    names(v) = c("Wald padj < 0.05", "LFC > 0 (Wald padj < 0.05)", "LFC < 0 (Wald padj < 0.05)")
    return(v)
}

get_dge_results = function(x, dge, constrasts, sign_adj_pvalue = 0.05, sign_fc = 1.5){
    contrast = constrasts %>%
        filter(Info == x) %>%
        select(-c(Info)) %>%
        unlist
    res = results(dge, contrast=as.numeric(contrast), alpha=0.05, test="Wald")
    return(as.data.frame(res) %>%
        rownames_to_column('genes') %>%
        mutate(sign_padj = padj < sign_adj_pvalue) %>%
        mutate(sign_padj_and_fc = (padj < sign_adj_pvalue & abs(log2FoldChange) >= log2(sign_fc))))
}

extract_DEG_log2FC = function(in_l, dir_path){
    # create dir if it does not exist
    full_dir_path = dir_path#paste("../results/dge/", dir_path, sep="")
    dir.create(full_dir_path, showWarnings = FALSE)
    l = list()
    # extract the log2FC of the genes with significant p-value and significant p-value + FC
    l$fc_deg = dplyr::data_frame(genes=character())
    l$sign_fc_deg = dplyr::data_frame(genes=character())
    for(i in names(in_l)){
        df = in_l[[i]]
        fc_deg = df %>%
            filter(sign_padj) %>%
            select(c(genes, log2FoldChange)) %>%
            rename(!!i:= log2FoldChange)
        l$fc_deg = l$fc_deg %>%
            full_join(fc_deg, by="genes")
        sign_fc_deg = df %>%
            filter(sign_padj_and_fc) %>%
            select(c(genes, log2FoldChange)) %>%
            rename(!!i:= log2FoldChange)
        l$sign_fc_deg = l$sign_fc_deg %>%
            full_join(sign_fc_deg, by="genes")
    }
    write.table(l$fc_deg, paste(full_dir_path, "fc_deg", sep=""), sep = "\t", quote = FALSE)
    write.table(l$sign_fc_deg, paste(full_dir_path, "sign_fc_deg", sep=""), sep = "\t", quote = FALSE)
    return(l)
}

nb_non_na = function(vec) return(sum(!is.na(vec)))
nb_pos = function(vec) return(sum(!is.na(vec) & vec > 0))
nb_neg = function(vec) return(sum(!is.na(vec) & vec < 0))

extract_DEG_stats = function(l, dir_path){
    l$stat = rbind(
        "All DEG (Wald padj < 0.05)"=l$fc_deg %>% select(-genes) %>% summarise_all(funs(nb_non_na)),
        "All over-expressed genes (Wald padj < 0.05 & FC > 0)"=l$fc_deg %>% select(-genes) %>% summarise_all(funs(nb_pos)),
        "All under-expressed genes (Wald padj < 0.05 & FC < 0)"=l$fc_deg %>% select(-genes) %>% summarise_all(funs(nb_neg)),
        "DEG (Wald padj < 0.05 & abs(FC) >= 1.5)"=l$sign_fc_deg %>% select(-genes) %>% summarise_all(funs(nb_non_na)),
        "Over-expressed genes (Wald padj < 0.05 & FC >= 1.5)"=l$sign_fc_deg %>% select(-genes) %>% summarise_all(funs(nb_pos)),
        "Under-expressed genes (Wald padj < 0.05 & FC <= -1.5)"=l$sign_fc_deg %>% select(-genes) %>% summarise_all(funs(nb_neg)))
    # plot stat
    mat = melt(as.data.frame(l$stat %>% rownames_to_column('type')))
    mat$variable = factor(mat$variable)
    mat$type = factor(mat$type)
    g = ggplot(mat, aes(x = reorder(variable, desc(variable)), y = reorder(type, desc(type)))) +
        labs(x = "", y = "") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y = element_text(size = rel(1.8)), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        geom_point(aes(size=value,col=value)) +
        scale_colour_gradient(low = "blue", high="red")
    print(g)
    # plot barplot of DEGs numbers
    s = l$stat %>%
        rownames_to_column('type') %>%
        filter(type %in% c("Over-expressed genes (Wald padj < 0.05 & FC >= 1.5)", "Under-expressed genes (Wald padj < 0.05 & FC <= -1.5)")) %>%
        mutate(type = gsub("genes \\(Wald padj < 0.05 & FC >= 1.5\\)", "", type)) %>%
        mutate(type = gsub("genes \\(Wald padj < 0.05 & FC <= -1.5\\)", "", type)) %>%
        mutate(type = gsub("expressed", "regulated", type)) %>%
        mutate(type = gsub("Over", "Up", type)) %>%
        mutate(type = gsub("Under", "Down", type)) #%>%
        #rename_all(funs(gsub("[A-Za-z\\-]+ vs [A-Za-z\\-]+ \\(", "", .))) %>%
        #rename_all(funs(gsub("\\)", "", .)))
    m = as.matrix(s %>% select(-c(type)))
    par(mar=c(10,3,3,3))
    a = barplot(m, main="", col=c("red", "darkblue"), beside=TRUE, las=1, border="white", axisnames = FALSE, ylim=range(pretty(c(0, m))))
    text(a, 0, srt = 60, adj=1.1, xpd = TRUE, labels = rep(colnames(m), each=2), cex=0.65)

    pdf(paste(dir_path, 'deg_nb.pdf', sep=""))
    par(mar=c(10,3,3,3))
    a = barplot(m, main="", col=c("red", "darkblue"), beside=TRUE, las=1, border="white", axisnames = FALSE, ylim=range(pretty(c(0, m))))
    text(a, 0, srt = 60, adj=1.1, xpd = TRUE, labels = rep(colnames(m), each=2), cex=0.65)
    dev.off()
    return(l)
}

get_boolean_df_wo_genes = function(df){
    return(1*as.data.frame(df %>% select(-genes) %>% transmute_all(funs(!is.na(.)))))
}

plot_sign_DEG_upset = function(l, nsets = 6){
    upset(get_boolean_df_wo_genes(l$fc_deg), nsets = nsets)
}

plot_sign_FC_DEG_upset = function(l, nsets = 6){
    upset(get_boolean_df_wo_genes(l$sign_fc_deg), nsets = nsets)
}

fit_proba_weighting_function = function(l, gene_length){
    # prepare a table with 1 for sign DEG and other 0 otherwise, for all comparisons, with the length of the genes
    genes = l$sign_fc_deg %>% pull(genes)
    gene_length_df = dplyr::data_frame(genes = names(gene_length), gene_length = gene_length)
    gene_vector = l$sign_fc_deg %>%
        select(-genes) %>%
        transmute_all(funs(1*!is.na(.))) %>%
        mutate(genes = genes) %>%
        full_join(gene_length_df, by='genes') %>%
        replace(., is.na(.), 0)
    gene_length_vec = gene_vector %>% select(c(genes, gene_length)) %>% deframe()
    # fit the probability weighting function
    comp = head(colnames(gene_vector),-2)
    l$pwf = lapply(
        comp,
        function(x){
            suppressMessages(nullp(
                gene_vector %>% select(c(genes, !!as.name(x))) %>% deframe(),
                'mm10',
                'geneSymbol',
                plot.fit=F,
                bias.data=gene_length_vec))
            }
        )
    names(l$pwf) = comp
    return(l)
}

extract_cat_de_genes = function(comp, interesting_cat, sign_fc_deg, full_go_genes, file_prefix){
    cat = interesting_cat %>%
        filter(!is.na(!!as.name(comp))) %>%
        pull(category)
    sign_deg = sign_fc_deg %>%
        filter(!is.na(!!as.name(comp))) %>%
        mutate(genes = toupper(genes)) %>%
        pull(genes)
    cat_de_genes = full_go_genes %>%
        filter(values %in% cat) %>%
        filter(ind %in% sign_deg)
    write.table(cat_de_genes, paste(file_prefix, gsub("[(),]", "", gsub(" ", "_", comp)), sep=""), sep="\t", quote=FALSE, row.names=FALSE)
    return(cat_de_genes)
}

get_significant_cat = function(mat, annot, pvalue_threshold){
    return(mat %>%
        mutate_at(vars(-category), funs(replace(., . > pvalue_threshold | . == 0, NA))) %>%
        filter_at(vars(-category), any_vars(!is.na(.))) %>%
        left_join(annot, by = 'category'))
}

get_interesting_cat = function(l, cat_type, pvalue_threshold = 0.05){
    # extract the adjusted p-values and ratios
    over_p_adjust = dplyr::data_frame(category = character())
    under_p_adjust = dplyr::data_frame(category = character())
    ratios = dplyr::data_frame(category = character())
    for(x in names(l$wall)) {
        over = dplyr::data_frame(!!x := p.adjust(l$wall[[x]][, "over_represented_pvalue"], method="BH"),
                                  category=l$wall[[x]]$category)
        over_p_adjust = over_p_adjust %>%
            full_join(over, by = 'category')
        under = dplyr::data_frame(!!x := p.adjust(l$wall[[x]][, "under_represented_pvalue"], method="BH"),
                                  category =l$ wall[[x]]$category)
        under_p_adjust = under_p_adjust %>%
            full_join(under, by = 'category')
        # compute ratio
        comp_ratio = l$wall[[x]] %>%
            mutate(!!x := numDEInCat/numInCat) %>%
            select(c(category, !!as.name(x)))
        ratios = ratios %>%
            full_join(comp_ratio, by="category")
    }
    l$ratios = ratios
    # extract the categories with a significant p-values (over or under represented pvalue) and add annotation
    annot = l$wall[[1]] %>% select(if(cat_type == "GO") c("category", "term", "ontology") else c("category"))
    l$over = get_significant_cat(over_p_adjust, annot, pvalue_threshold)
    l$under = get_significant_cat(under_p_adjust, annot, pvalue_threshold)
    return(l)
}

get_col_ramp_all_ont = function(mat, min_col, max_col){
    mat = mat %>%
        select(-c(category, term, ontology))
    values = unique(na.omit(unlist(mat)))
    if(length(values) > 1){
        cuts = cut(values,
            breaks=exp(log(10)*seq(log10(min(values)), log10(max(values)), len = 100)),
            include.lowest = TRUE)
        df = data.frame(values=values, color=colorRampPalette(c(min_col, max_col))(99)[cuts]) %>%
            rownames_to_column("comparison") %>%
            mutate(comparison = gsub(")[0-9]+", ")", comparison))
        return(df)
    }else if (length(values) > 0){
        df = data.frame(values=values, color=colorRampPalette(c(min_col, max_col))(99)[50]) %>%
            rownames_to_column("comparison") %>%
            mutate(comparison = gsub(")[0-9]+", ")", comparison))
    }else {
        return(data.frame(values=numeric(), color=character(), comparison=character()))
    }
}

extract_GO_terms = function(l, dir_path, go_term_to_exclude_fp){
    # GO analysis
    go_dir_path = dir_path#paste("../results/dge/", dir_path, "go/", sep="")
    dir.create(go_dir_path, showWarnings = FALSE)
    l$GO = list()
    # calculate the over and under expressed GO categories among the DE genes
    l$GO$wall = lapply(l$pwf, function(x) suppressMessages(goseq(x, 'mm10', 'geneSymbol')))
    # extract interesting pathways/categories and export them
    l$GO = get_interesting_cat(l$GO, "GO")
    write.table(l$GO$over, paste(go_dir_path, "full_over_represented_GO", sep=""), sep = "\t", quote = FALSE)
    write.table(l$GO$under, paste(go_dir_path, "full_under_represented_GO", sep=""), sep = "\t", quote = FALSE)
    # exclude GO terms
    go_term_to_exclude = read.table(go_term_to_exclude_fp, h = F)$V1
    l$GO$over = l$GO$over %>%
        filter(!category %in% go_term_to_exclude)
    write.table(l$over, paste(go_dir_path, "over_represented_GO", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
    l$GO$under = l$GO$under %>%
        filter(!category %in% go_term_to_exclude)
    write.table(l$GO$under, paste(go_dir_path, "under_represented_GO", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
    # extract list of genes involved in the over and under represented GO
    full_go_genes = stack(getgo(l$sign_fc_deg$genes, 'mm10', 'geneSymbol'))
    lapply(names(l$GO$wall),
        function(x) extract_cat_de_genes(x, l$GO$over, l$sign_fc_deg, full_go_genes, paste(go_dir_path, "over_repr_", sep = "")))
    lapply(names(l$GO$wall),
        function(x) extract_cat_de_genes(x, l$GO$under, l$sign_fc_deg, full_go_genes, paste(go_dir_path, "under_repr_", sep = "")))
    # get list of terms and description
    l$GO$terms = union(l$GO$over, l$GO$under) %>%
        select(c(category, term, ontology))
    l$GO$cat = l$GO$terms %>% pull(category)
    # get colors
    l$GO$over_repr_colors = get_col_ramp_all_ont(l$GO$over, "red", "lightpink")
    #l$GO$under_repr_colors = get_col_ramp_all_ont(l$GO$under, "blue", "lightblue")
    return(l)
}

get_number = function(bool_expr) return(sum(bool_expr*1))
get_repartition = function(dge_mat) return(cbind(get_number(dge_mat[,1]>0 & dge_mat[,2]>0), get_number(dge_mat[,1]>0 & dge_mat[,2]<0), get_number(dge_mat[,1]<0 & dge_mat[,2]>0), get_number(dge_mat[,1]<0 & dge_mat[,2]<0)))
get_repartition_3col = function(dge_mat) return(cbind(
    get_number(dge_mat[,1]>0 & dge_mat[,2]>0 & dge_mat[,3]>0),
    get_number(dge_mat[,1]>0 & dge_mat[,2]>0 & dge_mat[,3]<0),
    get_number(dge_mat[,1]>0 & dge_mat[,2]<0 & dge_mat[,3]>0),
    get_number(dge_mat[,1]>0 & dge_mat[,2]<0 & dge_mat[,3]<0),
    get_number(dge_mat[,1]<0 & dge_mat[,2]>0 & dge_mat[,3]>0),
    get_number(dge_mat[,1]<0 & dge_mat[,2]>0 & dge_mat[,3]<0),
    get_number(dge_mat[,1]<0 & dge_mat[,2]<0 & dge_mat[,3]>0),
    get_number(dge_mat[,1]<0 & dge_mat[,2]<0 & dge_mat[,3]<0)))

quantile_breaks = function(xs, n = 10) {
    breaks = quantile(unlist(xs), probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
}

plot_count_heatmap = function(genes, samples, annot){
    plot_heatmap(norm_counts, genes, samples, annot)
}

plot_heatmap = function(count, genes, samples, annot, show_rownames=FALSE){
    data = count[genes,samples]
    breaks = quantile_breaks(data, n = 11)
    pheatmap(data,
             cluster_rows=F,
             cluster_cols=F,
             show_rownames=show_rownames,
             show_colnames=F,
             annotation_col=annot,
             annotation_color=annot_colors,
             breaks=breaks,
             color=inferno(10))
}

plot_fc_heatmap = function(sign_fc_deg, annot){
    data_vector = unlist(sign_fc_deg)
    data_vector = data_vector[!is.na(data_vector)]
    breaks = quantile_breaks(data_vector, n = 11)
    data = sign_fc_deg
    data[is.na(data)] = 0
    pheatmap(data,
        cluster_rows=T,
        cluster_cols=F,
        show_rownames=F,
        show_colnames=F,
        annotation_col=annot,
        annotation_color=annot_colors,
        breaks=breaks,
        color=rev(brewer.pal(10, "RdBu")))
}

get_list = function(mapping){
    mapped_genes = mappedkeys(mapping)
    return(as.list(mapping[mapped_genes]))
}

capFirst = function(s) {
    paste(substring(s, 1, 1), tolower(substring(s, 2)), sep = "")
}

get_top_go = function(mat, full_mat, top_nb, comp){
    return(mat %>%
        top_n(top_nb, desc(!!as.name(comp))) %>%
        select(c(category, !!as.name(comp))) %>%
        full_join(full_mat, by="category"))
}

plot_top_go = function(l, ont, top_nb){
    l = l$GO
    # extract the top GO for the ontology
    over_mat = l$over %>% filter(ontology == ont)
    under_mat = l$under %>% filter(ontology == ont)
    over_top_go = dplyr::data_frame(category=character())
    under_top_go = dplyr::data_frame(category=character())
    for(comp in names(l$wall)){
        over_top_go = get_top_go(over_mat, over_top_go, top_nb, comp)
        under_top_go = get_top_go(under_mat, under_top_go, top_nb, comp)
    }
    over_top_go = over_top_go %>% mutate(type = "over")
    under_top_go = under_top_go %>% mutate(type = "under")
    conserved_go = c(over_top_go$category, under_top_go$category)
    # get the ratios
    ratios = melt(as.data.frame(l$ratios %>% filter(category %in% conserved_go))) %>%
        rename(ratio = value)
    # format the matrix
    mat = melt(as.data.frame(bind_rows(over_top_go, under_top_go))) %>%
        rename(p_value = value) %>%
        left_join(ratios, by = c("category", "variable")) %>%
        rename(comparison = variable) %>%
        left_join(l$terms, by = 'category') %>%
        mutate(term = factor(term)) %>%
        mutate(comparison = factor(comparison, levels=names(l$wall)))
    # plot
    size_scale_lim = c(min(mat$ratio), max(mat$ratio))
    col_scale_lim = c(min(mat$p_value, na.rm=T), max(mat$p_value, na.rm=T))
    plot1 = mat %>%
      filter(type == "over") %>%
      ggplot() +
      geom_point(aes(x=comparison, y=term, size=ratio, col=p_value)) +
      labs(x = "", y = "Over represented") +
      theme(axis.text.x = element_text(angle=90, hjust=1), axis.text.y = element_text(size = rel(.75))) +
      scale_colour_gradient(limits=col_scale_lim, low="red", high="blue", guide='none') +
      scale_size(limits=size_scale_lim, range=c(0.2,2))
    plot2 = mat %>%
      filter(type == "under") %>%
      ggplot() +
      geom_point(aes(x=comparison, y=term, size=ratio, col=p_value)) +
      labs(x = "", y = "Under represented") +
      theme(axis.text.x = element_blank(), axis.text.y = element_text(size = rel(.75))) +
      scale_colour_gradient(limits=col_scale_lim, low="red", high="blue") +
      scale_size(limits=size_scale_lim, range=c(0.2,2), guide='none')
    grid.newpage()
    grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "first"))
}

get_deg_colors = function(comp_deg, comp, connected_gene_colors, module_nb){
    deg_col = connected_gene_colors
    deg_col = connected_gene_colors
    sign_pos_FC = rownames(comp_deg$sign_fc_deg)[comp_deg$sign_fc_deg[,comp] > 0]
    sign_neg_FC = rownames(comp_deg$sign_fc_deg)[comp_deg$sign_fc_deg[,comp] < 0]
    deg_col[which(names(deg_col) %in% sign_pos_FC)] = module_nb + 1
    deg_col[which(names(deg_col) %in% sign_neg_FC)] = module_nb + 2
    return(deg_col)
}

get_col_ramp = function(mat, ont, min_col, max_col){
    mat = mat %>%
        filter(ontology==ont) %>%
        select(-c(category, term, ontology))
    values = na.omit(unlist(mat))
    if(length(values) > 0){
        cuts = cut(values,
            breaks=exp(log(10)*seq(log10(min(values)), log10(max(values)), len = 100)),
            include.lowest = TRUE)
        return(data.frame(values=sort(values), color=colorRampPalette(c(min_col, max_col))(99)[cuts]))
    }else{
        return(data.frame(values=numeric(), color=character()))
    }
}

extract_GO_ont = function(l, ont){
    l[[ont]] = list()
    l[[ont]]$go_terms = l$terms %>%
        filter(ontology==ont) %>%
        pull(category)
    l[[ont]]$over_repr_colors = get_col_ramp(l$over, ont, "red", "lightpink")
    l[[ont]]$under_repr_colors = get_col_ramp(l$under, ont, "blue", "lightblue")
    return(l)
}

create_GO_network = function(l, ont, ont_go){
    net = list()
    # Get GO term for ontology and that are at least once over-represented or under-represented
    net$go_terms = l$terms %>%
        filter(ontology==ont) %>%
        pull(category)
    # Get similarity between GO terms
    go_sim = mgoSim(net$go_terms, net$go_terms, semData=ont_go, measure="Wang", combine=NULL)
    # Extract adjency matrix by keeping only the distance > 0.5
    go_sim = as.data.frame(go_sim) %>%
        transmute_all(funs((. >0.5)*1))
    go_sim = as.matrix(go_sim)
    diag(go_sim) = 0
    # Keep only the GO with at least one connection
    #l$go_sim = l$go_sim[rowSums(l$go_sim)>0,rowSums(l$go_sim)>0]
    # Build the network
    net$go_net = graph_from_adjacency_matrix(go_sim, diag = FALSE, weighted = TRUE, mode="undirected")
    net$go_net_layout = layout_with_fr(net$go_net)
    # Extract the colors
    net$over_repr_colors = get_col_ramp(l$over, ont, "red", "lightpink")
    net$under_repr_colors = get_col_ramp(l$under, ont, "blue", "lightblue")
    l[[ont]] = net
    return(l)
}

get_GO_network_col = function(l, ont, comp){
    over = l$over %>%
        filter(ontology==ont) %>%
        filter(!is.na((!!as.name(comp)))) %>%
        select(c(category,!!as.name(comp))) %>%
        rename(values=!!as.name(comp)) %>%
        inner_join(l[[ont]]$over_repr_colors, by='values') %>%
        select(c(color, category))
    under = l$under %>%
        filter(ontology==ont) %>%
        filter(!is.na((!!as.name(comp)))) %>%
        select(c(category,!!as.name(comp))) %>%
        rename(values=!!as.name(comp)) %>%
        inner_join(l[[ont]]$under_repr_colors, by='values') %>%
        select(c(color, category))
    colors = over %>%
        union(under) %>%
        union(dplyr::data_frame(color = "white", category = l[[ont]]$go_terms)) %>%
        mutate(category=factor(category, levels=l[[ont]]$go_terms)) %>%
        arrange(color)
    return(colors)
}

plot_GO_network = function(l, ont, comp){
    # extract under and over represented GO and their value
    colors = get_GO_network_col(l, ont, comp)
    # plots
    plot(l[[ont]]$go_net,
         vertex.label=NA,
         vertex.size=4,
         vertex.color=colors$color,
         layout=l[[ont]]$go_net_layout)
}

get_GO_network_col_all_ont = function(l, comp){
    over = l$over %>%
        filter(!is.na((!!as.name(comp)))) %>%
        select(c(category,!!as.name(comp))) %>%
        rename(values=!!as.name(comp)) %>%
        left_join(l$over_repr_colors %>% filter(comparison == comp) %>% distinct(color, values), by='values') %>%
        select(c(color, category))
    under = l$under %>%
        filter(!is.na((!!as.name(comp)))) %>%
        select(c(category,!!as.name(comp))) %>%
        rename(values=!!as.name(comp)) %>%
        left_join(l$under_repr_colors %>% filter(comparison == comp) %>% distinct(color, values), by='values') %>%
        select(c(color, category))
    white = dplyr::data_frame(color = "white", category = l$cat) %>%
        filter(!category %in% over$category) %>%
        filter(!category %in% under$category)
    colors = over %>%
        union(under) %>%
        union(white) %>%
        mutate(category=factor(category, levels=l$cat)) %>%
        arrange(color)
    return(colors)
}

extract_KEGG_pathways = function(l, dir_path){
    # KEGG analysis
    kegg_dir_path = dir_path#paste("../results/dge/", dir_path, "kegg/", sep="")
    dir.create(kegg_dir_path, showWarnings = FALSE)
    l$KEGG = list()
    # calculate the over and under expressed KEGG pathways among the DE genes
    l$KEGG$wall = lapply(l$pwf, function(x) suppressMessages(goseq(x, 'mm10', 'geneSymbol', test.cats="KEGG")))
    # extract interesting pathways/categories and export them
    l$KEGG = get_interesting_cat(l$KEGG, "KEGG")
    write.table(l$KEGG$over, paste(kegg_dir_path, "over_represented_KEGG", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(l$KEGG$under, paste(kegg_dir_path, "under_represented_KEGG", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
    # extract list of genes involved in the over and under represented KEGG
    full_kegg_genes = stack(getgo(l$sign_fc_deg$genes, 'mm10', 'geneSymbol', fetch.cats="KEGG"))
    lapply(names(l$KEGG$wall),
        function(x) extract_cat_de_genes(x, l$KEGG$over, l$sign_fc_deg, full_kegg_genes, paste(kegg_dir_path, "over_repr_", sep = "")))
    lapply(names(l$KEGG$wall),
        function(x) extract_cat_de_genes(x, l$KEGG$under, l$sign_fc_deg, full_kegg_genes, paste(kegg_dir_path, "under_repr_", sep = "")))
    return(l)
}


plot_kegg_pathways = function(kegg_cats, fc_deg, dir_path){
    dir.create(dir_path, showWarnings = FALSE)
    tmp_fc_deg = fc_deg
    fc_deg = tmp_fc_deg %>% select(-genes) %>% as.matrix()
    rownames(fc_deg) = tmp_fc_deg %>% pull(genes)
    for(cat in kegg_cats){
        if(cat!="01100"){
            suppressMessages(pathview(gene.data=fc_deg, pathway.id=cat, species="Mus musculus", gene.idtype="SYMBOL", same.layer = F))
            if(dim(fc_deg)[2] > 1){
                file.rename(from=paste(dir_path,"mmu", cat, ".pathview.multi.png",sep=""),
                            to=paste(dir_path,"mmu", cat, ".pathview.multi.png",sep=""))
            }else{
                file.rename(from=paste(dir_path,"mmu", cat, ".pathview.png",sep=""),
                            to=paste(dir_path,"mmu", cat, ".pathview.png",sep=""))
            }
        }
    }
}

investigate_gene_set = function(mat){
    print(dim(mat)[1])
    print(cor.test(mat[,1],mat[,2]))
}

investigate_enrichement = function(set, all_genes){
    deg = 1*(all_genes %in% set)
    names(deg) = all_genes
    pwf = nullp(deg, 'mm10', 'geneSymbol', plot.fit=F)
    res = matrix(0,nrow=2,ncol=2, dimnames=list(c("over", "under"),c("GO","KEGG")))
    # GO
    GO_wall = goseq(pwf,'mm10', 'geneSymbol')
    over_represented_GO = GO_wall[GO_wall$over_represented_pvalue < 0.05,c("category","term","ontology")]
    under_represented_GO = GO_wall[GO_wall$under_represented_pvalue < 0.05,c("category","term","ontology")]
    res[1,1] = dim(over_represented_GO)[1]
    res[2,1] = dim(under_represented_GO)[1]
    # plot ontology barplot
    GO_ontology_counts = merge(count(over_represented_GO, var="ontology"), count(under_represented_GO, var="ontology"), by="ontology")
    rownames(GO_ontology_counts) = GO_ontology_counts[,1]
    GO_ontology_counts = GO_ontology_counts[,-1]
    colnames(GO_ontology_counts) = c("over_represented_GO", "under_represented_GO")
    GO_ontology_counts = as.matrix(GO_ontology_counts)
    barplot(t(GO_ontology_counts), beside = TRUE, col=c("green4", "red4"))
    legend("topright", c("over", "under"), fill=c("green4", "red4"))
    # KEGG pathways
    KEGG_wall = goseq(pwf,'mm10', 'geneSymbol', test.cats="KEGG")
    over_represented_KEGG = KEGG_wall[KEGG_wall$over_represented_pvalue < 0.05,]
    under_represented_KEGG = KEGG_wall[KEGG_wall$under_represented_pvalue < 0.05,]
    res[1,2] = dim(over_represented_KEGG)[1]
    res[2,2] = dim(under_represented_KEGG)[1]
    print(res)
}


order_by_min_na = function(data){
    na_nb = apply(is.na(data), 2, sum)
    return(data[order(data[,which.min(na_nb)]),])
}

plot_fc_heatmap_with_modules = function(data, fc_annot, connected_gene_colors){
    module_nb = length(unique(connected_gene_colors))
    # get genes per modules and order inside module by log2FC of the comparison with the less NA
    module_gene_fc_deg = matrix(0, ncol = 0, nrow = 0)
    for(i in 1:module_nb){
        genes_in_mod = names(connected_gene_colors)[connected_gene_colors == i & names(connected_gene_colors) %in% rownames(data)]
        module_gene_fc_deg = rbind(module_gene_fc_deg, order_by_min_na(data[genes_in_mod,]))
    }
    # add genes not in modules (also ordered)
    genes_in_mod_to_keep = rownames(module_gene_fc_deg)
    genes_not_in_mod = rownames(data)[!rownames(data) %in% genes_in_mod_to_keep]
    module_gene_fc_deg = rbind(module_gene_fc_deg, order_by_min_na(data[genes_not_in_mod,]))
    # add annotation
    annot_row = data.frame(module = as.factor(c(connected_gene_colors[genes_in_mod_to_keep], rep("No module", length(genes_not_in_mod)))))
    rownames(annot_row) = c(genes_in_mod_to_keep, genes_not_in_mod)
    mod_pal = c(head(pal2, -2),'white')
    names(mod_pal) = as.factor(c(1:module_nb,"No module"))
    annot_colors$module = mod_pal
    # get breaks for the colors
    data_vector = unlist(module_gene_fc_deg)
    data_vector = data_vector[!is.na(data_vector)]
    breaks = quantile_breaks(data_vector, n = 11)
    # plot heatmap
    pheatmap(module_gene_fc_deg,
            cluster_rows=F,
            cluster_cols=F,
            show_rownames=F,
            show_colnames=F,
            annotation_col=annot_col,
            annotation_row=annot_row,
            annotation_colors = annot_colors,
            breaks=breaks,
            color=rev(brewer.pal(10, "RdBu")))
}


order_rows = function(z){
    # perform hierarchical clustering inside the module to order the genes
    hc = hclust(dist(z), method = "complete")
    return(z[hc$order,])
}

get_genes_in_mod = function(mod_id, connected_gene_colors, gene_subset){
    # Get the genes in a given module and in the subset of genes
    return(names(connected_gene_colors)[connected_gene_colors == mod_id &
                                        names(connected_gene_colors) %in% gene_subset])
}


plot_z_score_heatmap_with_modules = function(z_scores, genes, col_order, annot_col, genes_in_modules, title, show_rownames=FALSE, fp=NULL, fontsize_row=10){
    data = z_scores
    genes = genes[genes %in% rownames(data)]
    # get genes ordered by modules
    genes_in_mod = c()
    for(x in names(genes_in_modules)){
        g_in_mod_id = intersect(genes_in_modules[[x]], genes)
        new_c = rep(x, length(g_in_mod_id))
        names(new_c) = g_in_mod_id
        genes_in_mod = c(genes_in_mod, new_c)
    }
    # order the z-score matrix by genes in genes_in_mod
    data = data[names(genes_in_mod), col_order]
    annot_colors$module=pal2
    # plot heatmap
    xx = pheatmap(data[names(genes_in_mod), col_order],
        cluster_rows=F,
        cluster_cols=F,
        show_rownames=show_rownames,
        show_colnames=F,
        annotation_col=annot_col,
        annotation_row=data.frame( module=genes_in_mod),
        annotation_colors = annot_colors,
        color=rev(brewer.pal(11, "RdBu")),
        breaks = seq(-3.5, 3.5, length=12),
        main = title,
        fontsize_row=fontsize_row)
    if(!is.null(fp)){
        pdf(fp)
        grid::grid.newpage()
        grid::grid.draw(xx$gtable)
        dev.off()
    }else{
        grid::grid.draw(xx$gtable)
    }
}

plot_z_score_heatmap = function(z_scores, de_genes, col_order, annot_col, title, col_for_clust, fp=NULL, show_rownames=FALSE, fontsize_row=10){
    # get z_score for the DE genes and correct column order
    de_genes = de_genes[de_genes %in% rownames(z_scores)]
    data = z_scores[de_genes,]
    # cluster rows
    hc = hclust(dist(data[,col_for_clust]), method = "complete")
    # plot
    xx = pheatmap(data[hc$order, col_order],
        cluster_rows=F,
        cluster_cols=F,
        show_rownames=show_rownames,
        show_colnames=F,
        annotation_col=annot_col,
        annotation_row=NULL,
        annotation_colors = annot_colors,
        color=rev(brewer.pal(11, "RdBu")),
        breaks = seq(-3.5, 3.5, length=12),
        fontsize_row=fontsize_row,
        main = title)
    if(!is.null(fp)){
        pdf(fp)
        grid::grid.newpage()
        grid::grid.draw(xx$gtable)
        dev.off()
    }else{
        grid::grid.draw(xx$gtable)
    }
}

plot_z_score_heatmap_2 = function(z_scores, de_genes, col_order, annot_col, title, col_for_clust, fp=NULL, show_rownames=FALSE, fontsize_row=10){
    # get z_score for the DE genes and correct column order
    de_genes = de_genes[de_genes %in% rownames(z_scores)]
    data = z_scores[de_genes,]
    # cluster rows
    hc = hclust(dist(data[,col_for_clust]), method = "complete")
    # plot
    xx = pheatmap(data[hc$order, col_order],
        cluster_rows=F,
        cluster_cols=F,
        show_rownames=show_rownames,
        show_colnames=F,
        annotation_col=annot_col,
        annotation_row=NULL,
        annotation_colors = annot_colors,
        color=inferno(100),
        border_color = "NA",
        breaks = seq(0, 5, length=85),
        fontsize_row=fontsize_row,
        main = title)
    if(!is.null(fp)){
        pdf(fp)
        grid::grid.newpage()
        grid::grid.draw(xx$gtable)
        dev.off()
    }else{
        grid::grid.draw(xx$gtable)
    }
}


plot_top_deg_in_modules = function(fc, comp, genes_in_modules){
    fc = fc %>%
        filter(!is.na((!!as.name(comp)))) %>%
        select(c(genes, !!as.name(comp)))
    mod_genes = lapply(genes_in_modules, function(mod){
        mod_fc = fc %>%
            filter(genes %in% mod)
        neg_genes = mod_fc %>%
            filter((!!as.name(comp)) < 0) %>%
            top_n(10, desc(!!as.name(comp))) %>%
            pull(genes)
        pos_genes = mod_fc %>%
            filter((!!as.name(comp)) > 0) %>%
            top_n(10, !!as.name(comp)) %>%
            pull(genes)
        return(list('pos_genes' = pos_genes, 'neg_genes' = neg_genes, 'max' = max(c(length(pos_genes), length(neg_genes)))))
    })
    # plot empty plot
    div = 2
    top_neg = 10*1/div
    top_pos = 2*top_neg + .5
    top = top_pos + 1
    w_offset = .5
    w = 2*length(genes_in_modules) + w_offset
    h = top
    prev_w = 7
    prev_h = 7
    options(repr.plot.width=w, repr.plot.height=h)
    par(mar = c(0,0,0,0))
    plot(c(0, w), c(0, h), type= "n", axes=FALSE, ann=FALSE)
    text(w/2, top, paste("Top DE genes in", comp), cex = 2.8, font = 2)
    text(0, top_neg + .5 + top_neg/2,"FC > 0", srt = 90, cex = 2)
    text(0, top_neg/2,"FC < 0", srt = 90, cex = 2)
    # parse modules
    t = sapply(1:length(mod_genes), function(i){
        mod = names(mod_genes)[i]
        mod_info = mod_genes[[mod]]
        # plot rect on the top
        rect((i-1)*2 + w_offset, top-.8, i*2 + w_offset, top - 1, col = pal2[mod], border = NA)
        text((i-1)*2 + 1 + w_offset, top-.7, mod, cex = 1.5, adj = c(0.5, 0))
        # add pos gene names to plot
        pos_genes = mod_info[['pos_genes']]
        if(length(pos_genes) > 0){
            t = sapply(1:length(pos_genes), function(x){
                text((i-1)*2 + w_offset, top_pos - x/div, pos_genes[x], cex = 1.2, adj = c(0, 0.5), offset = c(1,NA))
            })
        }
        # add neg gene names to plot
        neg_genes = mod_info[['neg_genes']]
        if(length(neg_genes) > 0){
            t = sapply(1:length(neg_genes), function(x){
                text((i-1)*2 + w_offset, top_neg - x/div, neg_genes[x], cex = 1.2, adj = c(0, 0.5), offset = c(1,NA))
            })
        }
    })
}