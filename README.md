# Clustering-and-HJ--Agg

Clustering script with can:

## Cluster specifically pi-stacked molecules if they are within the cutoff

cc.pi_stack_clustering(num_mol,pi_selection,cutoff,out_file_name)

num_mol = number of molecules 

pi_selection = selection containing each residue to loop over and cluster

cutoff = the cutoff distance to determine if 2 molecules are in the same cluster

out_file_name = output file name

## Cluster entire molcules with with multiple residues by clustering two molecules if any pair of residues is within the cutoff

cc.nres_cluster(num_mol,num_res_in_mol,cutoff,out_file_name)

num_mol = number of molecules

num_res_in_mol = number of residues in each molecule

cutoff = clustering cutoff distance

out_file_name = output file name

## Cluster using the above nres_cluster command but includes an added function which calculates number of J-agg,H-agg and random paris

cc.nres_cluster_w_hj(num_mol,num_res_in_mol,pi_selection,cutoff,out_file_name)

inputs same as above in cc.nres_cluster