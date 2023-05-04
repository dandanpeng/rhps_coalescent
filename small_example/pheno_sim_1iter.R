

#Doc Edge, 6/13/18
#Goal: Simulate phenotypic scores (possibly under selection) that are additive functions
# of genotypes at many loci. 


#Plan: simulate n_loci allele frequency trajectories.
#to do this, select an effect size from normal dist with mean 0 and sd
#equal to heritability over n_loci (from Martin et al). locus undergoes selection 
#according to effect size and selection gradient of trait (which changes through time). 
#Generate an allele frequency trajectory obeying the selection coefficient.
#Do this for all loci and save "phenotype" trajectory.
#Run through ms and save trees to estimate phenotype trajectory using coalescent times.
#Run haplotypes through RENT+ to see whether estimated trees recover phenotype trajectory.


trajs <- list()
sel.parts <- list()
eff_sizes <- numeric(0)
ss <- numeric(0)
curr.freqs <- numeric(0)


#need a constant that reflects a harmonic series from 1 to 2N-1, but only the
#terms where i/(2N) \in [.01, .99]. for N large this will be very close to the 
#analogous integral, which is ln(.99) - ln(.01) ~= 4.595
ser <- 1:(2*N-1)
harm.const <- sum(1/ser[ser >= .01*2*N & ser <= .99*2*N])


#Simulate and write n.loci allele-frequency trajectories. None of these should have fixed,
#so reject and do it again if fixed. Reject if minor allele frequency is < 0.01.
shift.achieved <- 0
while(shift.achieved == 0){
	pre.sel.freqs <- numeric(0)
	post.sel.freqs	<- numeric(0)
	for(i in 1:n.loci){
		fix <- 1
		while(fix == 1){
			eff_size <- rnorm(1,0,sqrt(herit * sd.trait^2 * harm.const /n.loci))
			s_locus <- eff_size * sel.intens / sd.trait #Charlesworth B3.7.7
			ss[i] <- s_locus
			p0 <- gen.neut.sfs.condit(1, 2*N, .01) #select from neutral sfs conditional on common.
			pre.sel.freqs[i] <- p0 
			sel.part <- sel.traj(p0, s = s_locus, N = N, t = t-t.off)
			post.sel.freqs[i] <- sel.part[nrow(sel.part),2]
			if(t.off > 0){
				post.drift <- neut.traj.time(p0 = sel.part[nrow(sel.part),2], N, t.off)
				sel.part <- rbind(sel.part, cbind(post.drift[,1] + t - t.off, post.drift[,2] ))			
			}
			if(sel.part[nrow(sel.part),2] >= 0.01 & sel.part[nrow(sel.part),2] <= 0.99){fix = 0}
		}
		eff_sizes[i] <- eff_size 
		curr.freqs[i] <- sel.part[nrow(sel.part),2]
		sel.parts[[i]] <- sel.part
	}	
	#Check whether the achieved shift is close to the target.
	trait.sd <- sqrt(sum(2 * curr.freqs * (1 - curr.freqs) *eff_sizes^2))
	sel.shift <- sum(2 * eff_sizes * (post.sel.freqs - pre.sel.freqs)) / trait.sd
	target.shift <- (sel.intens * sd.trait) * (t - t.off) * 2 * N
	if(target.shift*.95 <= sel.shift & sel.shift <= target.shift*1.05){shift.achieved <- 1}
#	print("trait.attempted")
}


for(i in 1:n.loci){
	sel.part <- sel.parts[[i]]
	driftup <- neut.traj(pre.sel.freqs[i], N, loss = TRUE)
	traj.fwd <- rbind(cbind(driftup[,1], rev(driftup[,2])), cbind(sel.part[-1,1] + max(driftup[,1]), sel.part[-1,2] ))
	traj <- traj.fwd
	traj[,2] <- rev(traj.fwd[,2])
	trajs[[i]] <- traj
}

eff_sizes_unscaled <- eff_sizes
eff_sizes <- eff_sizes / trait.sd #The SD of the polygenic score is set
#to be 1 in the present


rm(ser)
rm(harm.const)
rm(eff_size)
rm(s_locus)
rm(p0)
rm(sel.part)
#rm(post.drift)
rm(post.sel.freqs)
rm(pre.sel.freqs)


mat.trajs <- matrixify.list.of.trajs(trajs) #matrix of derived allele freq at each locus

phen.traj <- as.numeric(2 * eff_sizes %*% t(mat.trajs))  #phenotypes vector
pt.time <- seq(0, by = 1/(2*N), length.out = max(sapply(trajs, length)/2) ) #desired time interval




n_ders <- numeric(n.loci)

ms_trees_list <- list()
rent_trees_list <- list()
relate_trees_list <- list()
tsinfer_trees_list <- list()
argweaver_samples_list <- list()
argweaver_trees_list <- list()
ms_haplotypes_list <- list()


#for each locus, write trajectory, run ms to simulate a sample and tree,
#and run rent+ to infer tree. Save both the ms and the rent+ trees, as well
#as the number of derived alleles for each locus.

for(i in 1:n.loci){
	write.traj(traj.fn, trajs[[i]])
	curr.freq.der <- curr.freqs[i]
	n_der <- rbinom(1, n_chroms, curr.freq.der) #number of chroms with derived allele.
	if(n_der == 0){n_der <- n_der + 1}
	if(n_der == n_chroms){n_der <- n_der - 1}
	n_ders[i] <- n_der
	
	
	counter = 0
	success = 1 # when RENT+ running fails
	
	while (success == 1 & counter <= 100){
    	#run mssel
    	ms.string <- paste(ms_dir, "mssel ", as.character(n_chroms), " 1 ", as.character(n_chroms - n_der), " ", as.character(n_der), " ", traj.fn, " ", sel_site, " -r ", as.character(4*N*(1 - dbinom(0,len_hap, re))), " ", as.character(len_hap), " -t ", as.character(4*N* (1 - dbinom(0, len_hap, u))), " -T -L > ", msout.fn, sep = "")
    	system(ms.string)
    
        ######################## format processing #############################
        
    	#cut down mssel output to just the part rent+ needs
    	process.string <- paste("sed '/positions: /,$!d' ", msout.fn, " | sed 's/.*positions: //' > ", rent_in_fn, sep = "")
    	system(process.string)
	
      #run rent+
      rent.string <- paste("java -jar ", rentplus_fn, " -t -l ", as.character(len_hap), " ", rent_in_fn, " ", msout.fn, sep = "")
      success = system(rent.string)
      counter = counter + 1
    }
	    
    #transfer rent+ input to RELATE recognized format
    source("relate_input.R")

	  #transfer rent+ input to argweaver recognized format
	  #use_python("")
    aw_input("rent_in.txt", len_hap, temp)   


    ms_haplotypes_list[[i]] <- readLines(rent_in_fn) #save the haplotypes produced by ms (and fed to rent+) in a list entry.
    #the result will be a character vector. The first entry in the vector is relative (0-1) positions, with the selected
    #site at 0.05. The remaining entries are string haplotypes of 0s and 1s.
    
    ######################## run softwares #################################
    #run RELATE
    setwd(temp)
    relate.string <- paste("../relate/Relate --mode All -m", as.character(u),
                           "-N", as.character(N*2),
                           "--haps data.haps", 
                           "--sample data.sample",
                           "--map data.map",
                           "--output data")
    system(relate.string)
    setwd("..")
    
    #run tsinfer
    source("TSinfer.R", verbose = TRUE)
    nwtree = divide_newick(nwtree,4*N) # 4N gives the resulting tree a height of approx. 2
     
    #run ARGweaver
    argweaver.string <- paste("arg-sample", "-s", argweaver_in_fn, "-N", N, "-r", re, "-m", u, "-n", as.character(aw_sample * 10), "--overwrite -o", argweaver_out_fn)
    system(argweaver.string)
    
    ############################# extract the local tree ##########################

    #Read in trees at selected site in rent+
    rent_trees_fn <- paste(rent_in_fn, ".trees", sep = "")
    rent_trees <- readLines(rent_trees_fn)
    #which tree is for the selected site?
    rent_sel_ind <- grep(paste(as.character(sel_site), "\t", sep = ""), rent_trees)
    
    if(length(rent_sel_ind) > 1){
        rent_sel_ind <- rent_sel_ind[1]
    }
    
    rent_sel_tree_newick <- sub(paste(as.character(sel_site), "\t", sep = ""), "", rent_trees[rent_sel_ind])
    rent_trees_list[[i]] <- read.tree(text = paste(rent_sel_tree_newick, ";", sep = ""))


    #Read in trees at selected site in RELATE
    setwd(temp)
    relate_sel_tree <- paste("../relate/RelateExtract --mode AncToNewick",
                          "--anc data.anc",
                          "--mut data.mut",
                          "--bp_of_interest 100000",
                          "--output data")
    print("execute relateextract")
    system(relate_sel_tree)
    relate_trees_list[[i]] <- read.tree(text = readLines("data.newick"))
    setwd("..")

    for(k in 1:length(relate_trees_list[[i]]$tip.label)){
      relate_trees_list[[i]]$tip.label[k] <- as.character(as.numeric(relate_trees_list[[i]]$tip.label[k]) + 1)
    }

    #Read in tree at selected site in tsinfer
    tsinfer_sel_tree_newick <- nwtree
    tsinfer_trees_list[[i]] <- read.tree(text = tsinfer_sel_tree_newick)

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
	ms_trees_list[[i]] <- read.tree(text = msoutlines[4 + ms_to_extract])
	
	print(paste("simulation locus", as.character(i), "of trait", as.character(k), "complete"))
}

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




