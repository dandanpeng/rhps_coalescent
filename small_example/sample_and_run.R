sample_indices_ls <- list()

ms_trees_list <- list()
rent_trees_list <- list()
relate_trees_list <- list()
tsinfer_trees_list <- list()
argweaver_samples_list <- list()
argweaver_trees_list <- list()
ms_haplotypes_list <- list()

for(i in 1:n.loci){
    write.traj(traj.fn, trajs[[i]])
    n_der <- n_ders[i]
    nsite <- nsites[i]
    theta <- theta_ls[i]
    
    counter = 0
    success = 1
    
    while (success == 1 & counter <= 100){
        ms.string <- paste(ms_dir, "mssel ", as.character(2000), " 1 ", as.character(2000 - n_der), " ", as.character(n_der), " ", traj.fn, " ", sel_site, " -r ", as.character(nsite), " ", as.character(len_hap), " -t ", as.character(theta), " -T -L > ", msout.fn, sep = "")
        system(ms.string)
    
        #cut down mssel output to just the part rent+ needs
        process.string <- paste("sed '/positions: /,$!d' ", msout.fn, " | sed 's/.*positions: //' > ", rent_in_fn, sep = "")
        system(process.string)
        
        # sample from mssel output (actually transferred mssel ouput which can be fed to RENT+) 
        sample_rent_input(rent_in_fn, rent_sample_in, sample_num, i)
        
        # read the haplotypes from the text file into a matrix
        hap_mat <- hap_to_mat(rent_sample_in)
        
        #################################### RENT+ #################################
        rent.string <- paste("java -jar ", rentplus_fn, " -t -l ", as.character(len_hap), " ", rent_sample_in, sep = "")
        success = system(rent.string)
        counter = counter + 1
    }
    
    #pull trees from ms
  	msoutlines <- readLines(msout.fn)
  	dists <- numeric(0)
  	for(m in 1:(length(msoutlines) - 4)){
  		suppressWarnings(dists[m] <- as.numeric(substring(strsplit(msoutlines[m + 4], "]")[[1]][1],2)))
  	}
  	dists <- dists[!is.na(dists) & dists != Inf]
  	posits.rs <- cumsum(dists)
  	ms_to_extract <- which.max(posits.rs[posits.rs < sel_site]) + 1
  	if(length(ms_to_extract)  == 0){ms_to_extract <- 1}
  	pre_tree <- read.tree(text = msoutlines[4 + ms_to_extract])
  	ms_trees_list[[i]] <- drop.tip(pre_tree, setdiff(as.character(c(1:2000)), as.character(sample_indices_ls[[i]])))
	
    print("start to run rent+")
    
   
    #Read in trees at selected site in rent+
    rent_trees_fn <- paste(rent_sample_in, ".trees", sep = "")
    rent_trees <- readLines(rent_trees_fn)
    
    #which tree is for the selected site?
    rent_dists <- numeric(0)
    for(r in 1:length(rent_trees)){
        rent_dists[r] <- strsplit(rent_trees[r], "\t")[[1]][1]
    }
    rent_extract_index <- which.min(abs(as.integer(rent_dists) - sel_site))
    rent_to_extract <- strsplit(rent_trees[rent_extract_index], "\t")[[1]][1]    

    rent_sel_ind <- grep(paste(rent_to_extract, "\t", sep = ""), rent_trees)
    
    if(length(rent_sel_ind) > 1){
        rent_sel_ind <- rent_sel_ind[1]
    }
    
    rent_sel_tree_newick <- sub(paste(rent_to_extract, "\t", sep = ""), "", rent_trees[rent_sel_ind])
    rent_trees_list[[i]] <- read.tree(text = paste(rent_sel_tree_newick, ";", sep = ""))
    
    
    snp_pos <- recover_pos(rent_sample_in, len_hap)
    close_sel_site <- snp_pos[which.min(abs(snp_pos) - sel_site)]
    derived_info <- hap_to_mat(rent_sample_in)
    
    ############################### RELATE #############################
    #transfer rent+ input to RELATE recognized format
    relate_input(derived_info, snp_pos)
    
    print("start to run RELATE")
   
    setwd(new_temp)
    relate.string <- paste("../relate/Relate --mode All -m", as.character(u),
                           "-N", as.character(N*2),
                           "--haps data.haps", 
                           "--sample data.sample",
                           "--map data.map",
                           "--output data")
    system(relate.string)
    
    #Read in trees at selected site in RELATE
    relate_sel_tree <- paste("../relate/RelateExtract --mode AncToNewick",
                          "--anc data.anc",
                          "--mut data.mut",
                          "--bp_of_interest", as.character(close_sel_site),
                          "--output data")
    print("execute relateextract")
    system(relate_sel_tree)
    relate_trees_list[[i]] <- read.tree(text = readLines("data.newick"))
    setwd("..")
    
    for(k in 1:length(relate_trees_list[[i]]$tip.label)){
      relate_trees_list[[i]]$tip.label[k] <- as.character(as.numeric(relate_trees_list[[i]]$tip.label[k]) + 1)
    }
    
    ################################ tsinfer #############################

    #run tsinfer
    print("start to run tsinfer")
    #source("TSinfer.R")
    nwtree = tsinfer_input(derived_info, snp_pos, close_sel_site)
    nwtree = divide_newick(nwtree,4*N) # 4N gives the resulting tree a height of approx. 2
    
    #Read in tree at selected site in tsinfer
    tsinfer_sel_tree_newick <- nwtree
    tsinfer_trees_list[[i]] <- read.tree(text = tsinfer_sel_tree_newick)
    
    ###############################  ARGWeaver #############################
    
    print("start to run ARGWeaver")
    
    #transfer rent+ input to argweaver recognized format
    aw_input("rent_sample_in.txt", snp_pos, len_hap, new_temp)   
    
    #run ARGweaver
   
    argweaver.string <- paste("arg-sample", "-s", argweaver_in_fn, "-N", N, "-r", re, "-m", u, "-n", as.character(aw_sample * 10), "--overwrite -o", argweaver_out_fn)
    system(argweaver.string)
    
    #Read in tree at selected site in ARGweaver
    for(sample in 0:aw_sample){
        aw.smc <- read.table(paste(argweaver_out_fn, ".", as.character(10*sample), ".smc.gz", sep = ""), skip = 2, header = FALSE, fill = TRUE)
        sel_row <- aw.smc[aw.smc$V1 == "TREE" & aw.smc$V2 <= sel_site & aw.smc$V3 >= sel_site, ]
        sel_tree <- read.tree(text = sel_row$V4)
        ori_label <- sel_tree$tip.label
        sel_tree$tip.label <- as.character(as.numeric(ori_label) + 1) # the original tip labels are zero-based numbering and we change them to one-based 
        argweaver_samples_list[[sample + 1]] <- sel_tree
    }
    argweaver_trees_list[[i]] <- argweaver_samples_list

    ms_haplotypes_list[[i]] <- readLines(rent_sample_in) #save the haplotypes produced by ms (and fed to rent+) in a list entry.
    #the result will be a character vector. The first entry in the vector is relative (0-1) positions, with the selected
    #site at 0.05. The remaining entries are string haplotypes of 0s and 1s.ts
    
	print(paste("simulation locus", as.character(i), "of trait", as.character(k), "complete"))
}

save.image("test.RData")
#Run through list and check for sites where 2+ variants had same coordinates as
#selected site. At those places, check which tree(s) are monophyletic for derived tips.
#if more than one, select randomly from among them.

check_samecor_sel_site(ms_trees_list)
check_samecor_sel_site(rent_trees_list)
check_samecor_sel_site(relate_trees_list)
check_samecor_sel_site(tsinfer_trees_list)
aw_check_samecor_sel_site(argweaver_trees_list)

##rescale branch length of relate trees
tmrca.relate <- numeric(0)

for(i in 1:n.loci){
  tmrca.relate[i] <- max(branching.times(relate_trees_list[[i]]))
}


for(i in 1:n.loci){
  relate_trees_list[[i]]$edge.length <- relate_trees_list[[i]]$edge.length*2/mean(tmrca.relate)
}

#rescale branxh length of tsinfer trees
for(i in 1:n.loci){
  tsinfer_trees_list[[i]]$edge.length <- tsinfer_trees_list[[i]]$edge.length*2
}

##rescale branch length of argweaver trees
##The original branch lengths are given in units of generations.
for(i in 1:n.loci){
    for(j in 1:length(argweaver_trees_list[[i]])){
        argweaver_trees_list[[i]][[j]]$edge
        argweaver_trees_list[[i]][[j]]$edge.length <- argweaver_trees_list[[i]][[j]]$edge.length/(2*N)
    }
}



#split into ancestral and derived trees and retrieve coalescence times.


split_tree("ms", ms_trees_list)
split_tree("rent", rent_trees_list)
split_tree("relate", relate_trees_list)
split_tree("tsinfer", tsinfer_trees_list)
aw_split_tree("argweaver", argweaver_trees_list)

#store coalesenct times of the whole tree, anc tree and der tree
coal_time_ls(ms_trees_list, anc_trees_ms, der_trees_ms, "ms")
coal_time_ls(rent_trees_list, anc_trees_rent, der_trees_rent, "rent")
coal_time_ls(relate_trees_list, anc_trees_relate, der_trees_relate, "relate")
coal_time_ls(tsinfer_trees_list, anc_trees_tsinfer, der_trees_tsinfer, "tsinfer")
aw_coal_time_ls(argweaver_trees_list, anc_trees_argweaver, der_trees_argweaver)

# store the number of lineages of ref allele (col1) and alt allel (col2) into a list of matrix
lins_ls(ms_trees_list, "ms")
lins_ls(rent_trees_list, "rent")
lins_ls(relate_trees_list, "relate")
lins_ls(tsinfer_trees_list, "tsinfer")
aw_lins_ls(argweaver_trees_list, times.c.argweaver)

fn_str <- paste("loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), sep = "")


save.image(paste(out_dir, "/", fn_str, ".RData", sep = ""))


	


    

