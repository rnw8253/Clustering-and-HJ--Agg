import sys
import numpy as np
import MDAnalysis as mda
import math

class Cluster():

    def computePbcDist2(self,r1,r2,box):
        dist2 = 0

        for j in range(0,3):
                temp = r1[j]-r2[j]
                if temp < -box[j]/2.0:
                        temp += box[j]
                elif temp > box[j]/2.0:
                        temp -= box[j]
                dist2 += temp*temp

        return dist2;


    def pi_stacked_nneighbor_matrix(self, sel, cutoff, num_mol, min_dist=True):
        # Determine which pi stacking cores are within cutoff
        # if min_dist == True then cluster_matrix will be filled with their distances else its filled with 1
        
        # selection of all the molecules that will be analyzed 
	
        cluster_matrix = np.zeros((num_mol,num_mol),dtype=np.float)
        
        # loop over each pair of molecules
        for i,res1 in enumerate(sel.residues):
            for j,res2 in enumerate(sel.residues):
                if res2.resid > res1.resid:
                    # Calculate the distance between each pair of moleculse
                    dist = math.sqrt(self.computePbcDist2(res1.center_of_mass(),res2.center_of_mass(),u.dimensions[:3]))
                    # fill cluster_matrix marking which pairs are within the cutoff distance
                    if(min_dist):
                        if dist <= cutoff:
                            cluster_matrix[i,j] = dist
                            cluster_matrix[j,i] = dist
                    else:
                        if dist <= cutoff:
                            cluster_matrix[i,j] = 1
                            cluster_matrix[j,i] = 1
	                
	    
        return cluster_matrix
	
    def min_res_dist_nneighbor_matrix(self, sel, num_mol, cutoff, min_dist=False):
        # Determine which pi stacking cores are within cutoff
        # if min_dist == True then cluster_matrix will be filled with their distances else its filled with 1
	    
        cluster_matrix = np.zeros((num_mol,num_mol),dtype=np.int)
        
        # loop over each pair of molecules
        for i in range(num_mol-1):
            for j in range(i+1,num_mol):
                min_dist = 100
                # loop over each pair of residues within the molecule parir
                for res1 in sel[i].residues:
                    for res2 in sel[j].residues:
                        # Calculate the distance between each pair of residues
                        dist = math.sqrt(self.computePbcDist2(res1.center_of_mass(),res2.center_of_mass(),u.dimensions[:3]))
                        # fill cluster_matrix marking which pairs are within the cutoff distance
                        if(min_dist):
                            if dist < min_dist:
                                min_dist = dist
                                if min_dist <= cutoff:
                                    cluster_matrix[i,j] = min_dist
                                    cluster_matrix[j,i] = min_dist
                        else:
                            if dist <= cutoff:
                                cluster_matrix[i,j] = 1
                                cluster_matrix[j,i] = 1
	
        return cluster_matrix
	
    def cluster_alg(self, num_mol,cluster_matrix):
        num_clusts=0
        sel = np.zeros(num_mol,dtype=np.int)
        todo = np.zeros(num_mol,dtype=np.int)
        cluster = np.zeros((num_mol,num_mol),dtype=np.int)
        cluster_size = np.zeros(num_mol,dtype=np.int)
        for i in range(num_mol):
            if sel[i] == 0:
                sel[i] = 1
                count1=0
                count2=0
                cluster[num_clusts][0] = i
                todo[0]=i
                cluster_size[num_clusts] += 1
                while count1 <= count2:
                    j=todo[count1]
                    for k in range(num_mol):
                        if sel[k] == 0 and cluster_matrix[k,j] != 0:
                            sel[k]=1
                            count2 += 1
                            todo[count2]=k
                            cluster[num_clusts][count2] = k
                            cluster_size[num_clusts] += 1
                    count1 += 1
                num_clusts += 1
	
        return cluster,cluster_size,num_clusts
	
    def hj_cluster(self, cluster_matrix,cluster,cluster_size,num_clusts):
        # calculate the total number of H, and random aggregates
        # method adapted from:
	    
        pi_6=0.8660254038
	
        # loop over every cluster
        for clust in range(num_clusts):
            # loop over every pair of molcules within the cluster within the cutoff
            for mol1 in range(cluster_size[clust]-1):
                for mol2 in range(mol1+1,cluster_size[clust]):
                    if cluster_matrix[cluster[clust][mol1],cluster[clust][mol2]] != 0:
                        # define vector pointing from COM of one molecule to another                                                                           
                        COM = sel[cluster[clust][mol1]].center_of_mass() - sel[cluster[clust][mol2]].center_of_mass()
                        for l in range(3):
                            if COM[l] < -hbox[l]:
                                COM[l] += box[l]
                            elif COM[l] > hbox[l]:
                                COM[l] -= box[l]
                        # define vectors pointing from NA to NB for both molecules in pair                                                                     
                        mol1 = NA_sel[cluster[clust][mol1]].atoms[0].position - NB_sel[cluster[clust][mol1]].atoms[0].position
                        mol1_norm = mol1/np.linalg.norm(mol1)
                        mol2 = NA_sel[cluster[clust][mol2]].atoms[0].position - NB_sel[cluster[clust][mol2]].atoms[0].position
                        mol2_norm = mol2/np.linalg.norm(mol2)
	
                        # define vector perpendicular to ring system                                                                                           
                        ring1 = C5_sel[cluster[clust][mol1]].atoms[0].position - C1_sel[cluster[clust][mol1]].atoms[0].position
                        ring2 = C3_sel[cluster[clust][mol1]].atoms[0].position - C1_sel[cluster[clust][mol1]].atoms[0].position
                        Z = np.cross(ring2,ring1)
                        Z /= np.linalg.norm(Z)
                        # calculate delta-Z                                                                                                                    
                        angle = abs(np.dot(mol1_norm,mol2_norm))
                        if angle >= pi_6:
                            delta_Z = abs(np.dot(COM,Z))
                        # calculate delta-D                                                                                                                    
                            if delta_Z < 5:
                                delta_D = abs(np.dot(COM,mol1_norm))
                                if delta_D < 14:
                                    # calculate stacking angle                                                                                     
                                    phi = math.atan(delta_Z/delta_D)
                                    phi_trans = 0.5585053606
                                    # calculate number of h and j aggregates                                                                       
                                    if phi <= phi_trans:
                                        jagg += 1
                                    if phi > phi_trans:
                                        hagg += 1
                                else:
                                    random_agg += 1
                            else:
                                random_agg += 1
                        else:
                            random_agg += 1
	
        # calculate the number of pairs and                                                                                                                                    
        total_agg = jagg + hagg + random_agg
        percent_j = np.float(jagg/total_agg)
        percent_h = np.float(hagg/total_agg)
        percent_r = np.float(random_agg/total_agg)
	
        return jagg,hagg,random_agg,total_agg_percent_j,percent_h,percent_r
	
    def calc_num_clust_sizes(self,num_mol,num_clusts,cluster_size):
        # calculate the total number of each size cluster
        size = np.zeros(num_mol,dtype=np.int)
        for k in range(num_clusts):
            size[cluster_size[k]-1] += 1
        return size
	
    def write_clusters(self,size,num_mol,num_clust,out):
        # write out the number of clusters and the number of each sized clusters
        out.write('%d' %(num_clust))
        for l in range(num_mol-1):
            out.write(' %d' %(size[l]))
            l += 1
        out.write(' %d\n' %(size[num_mol-1]))
    def write_hj_agg_clust(self,out2,agg,hagg,random_agg,total_agg_percent_j,percent_h,percent_r):
        out2.write("%8d %8d %8d %10.5f %8d %10.5f %8d  %10.5f\n" %(m, total_agg, jagg, percent_j, hagg, percent_h, random_agg, percent_r))
        return l

    def nres_cluster_w_hj(self,num_mol,num_res_in_mol,pi_selection,cutoff,out_file_name):
        # Script to run nearest-residue clustering and hj analysis on clusters
        sel = []
        c_sel = []
        NA_sel = []
        NB_sel = []
        C1_sel = []
        C3_sel = []
        C5_sel = []

        # create a list of molecules to cluster using num_mol and rum_res_in_mol
        print "Creating selections for all molecules"
        for i in range(num_mol):
            resA = i*num_res_in_mol + 1
            resB = i*num_res_in_mol + num_res_in_mol
            sel.append(u.select_atoms("resid %s:%s" %(resA,resB))) 

        # create selections for hj analysis
        pi_stacked_sel = u.select_atoms(pi_selection)
        print "Creating selections for hj_agg"
        # create a list of atoms from which the vectors for determing hj aggregation will be caluculated from
        # A vector from NA to NB should point across the pi-stacking surface
        # C1, C3, and C5 should form a plane on the pi stacking surface
        for res in pi_stacked_sel.resids:
            print "Residue number: %s" %(res)
            c_sel.append(u.select_atoms("resid %s and name C*" %(res)))
            NA_sel.append(u.select_atoms("resid %s and name NA" %(res)))
            NB_sel.append(u.select_atoms("resid %s and name NB" %(res)))
            C1_sel.append(u.select_atoms("resid %s and name C13" %(res)))
            C3_sel.append(u.select_atoms("resid %s and name C19" %(res)))
            C5_sel.append(u.select_atoms("resid %s and name C23" %(res)))        
            

        out = open("%s.hj_clust.dat" %(out_file_name),"w")
        for ts in u.trajectory:
            box = u.dimensions[:3]
            hbox = u.dimensions[:3]/2.
            print ts.frame
            cluster_matrix = cc.min_res_dist_nneighbor_matrix(sel, num_mol, cutoff, min_dist=False)
            cluster,cluster_size,num_clusts = cc.cluster_alg(num_mol, cluster_matrix)
            size = cc.calc_num_clust_sizes(num_mol,num_clusts,cluster_size)
            agg,hagg,random_agg,total_agg_percent_j,percent_h,percent_r = cc.hj_cluster(cluster_matrix,cluster,cluster_size,num_clusts)
            cc.write_clusters(size, num_mol, num_clusts, out)
        
        out.close

    def nres_cluster(self,num_mol,num_res_in_mol,cutoff,out_file_name):
        sel = []
        print "Creating selections for all molecules"
        for i in range(num_mol):
            resA = i*num_res_in_mol + 1
            resB = i*num_res_in_mol + num_res_in_mol
            sel.append(u.select_atoms("resid %s:%s" %(resA,resB))) 

        # Script to run nearest-residue clustering 
        out = open("%s.nres_clust.dat" %(out_file_name),"w")
        for ts in u.trajectory:
            box = u.dimensions[:3]
            hbox = u.dimensions[:3]/2.
            print ts.frame
        
            cluster_matrix = cc.min_res_dist_nneighbor_matrix(sel, num_mol, cutoff)
            cluster,cluster_size,num_clusts = cc.cluster_alg(num_mol, cluster_matrix)
            size = cc.calc_num_clust_sizes(num_mol,num_clusts,cluster_size)
            cc.write_clusters(size, num_mol, num_clusts, out)
        
        out.close

    def pi_stack_clustering(self,num_mol,pi_selection,cutoff,out_file_name):
        print "Creating selections for pi-stacking molecules"
        pi_stacked_sel = u.select_atoms(pi_selection)   
        # Script to run pi-stacked clustering
        out = open("%s.pi_clust.dat" %(out_file_name),"w")
        for ts in u.trajectory:
            box = u.dimensions[:3]
            hbox = u.dimensions[:3]/2.
            print ts.frame
            cluster_matrix = cc.pi_stacked_nneighbor_matrix(pi_stacked_sel, 5.5, num_mol)
            cluster,cluster_size,num_clusts = cc.cluster_alg(num_mol, cluster_matrix) 
            size = cc.calc_num_clust_sizes(num_mol,num_clusts,cluster_size)
            cc.write_clusters(size, num_mol, num_clusts, out)

        out.close

###################################
########## Main Program ###########
###################################

top_file = sys.argv[1]
traj_file = sys.argv[2]
out_file_name = sys.argv[3]

u = mda.Universe(top_file,traj_file)

num_mol = 30
num_res_in_mol = 17
cutoff = 5.5
pi_selection = "resname PDI"



cc = Cluster()


#cc.pi_stack_clustering(num_mol,pi_selection,cutoff,out_file_name)
cc.nres_cluster(num_mol,num_res_in_mol,cutoff,out_file_name)
#cc.nres_cluster_w_hj(num_mol,num_res_in_mol,pi_selection,cutoff,out_file_name)
