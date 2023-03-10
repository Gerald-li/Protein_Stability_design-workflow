# Protein_Stability_design-workflow
This is a process of protein stability modification designed by exploring weak regions of a protein.

Required software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Python v3.7.7 (Anaconda installation recommended)
- Rosetta v3.10 or later
- Pymol v2.4
- openmpi v4.0.5
- numpy v1.18.1 or later
- CNA web server (Constraint Network Analysis, https://cpclab.uni-duesseldorf.de/cna/main.php)

Usage
^^^^^
Input informaiton:
	Input structure: A Pdb structure (prepared structure and standard residues are allowed) is needed for stability calculation. 
	Input CNA dat: A dat file (ends with 'local_indices.dat') downloaded from the result page of the CNA web server is needed for Step 3 in the following workflow.

The main steps for this workflow:
	
	Step 1: - run " python step1_findholes.py input.pdb cpunumbers(how many cpus used for mpi) radiuscutoff(default = 2.5)"
			#for example "python step1_findholes.py 1PS6_chainA.pdb 10 2.5"  
			
		   
	Step 2: - run "python step2_detectholetype.py input.pdb holes.pdb(last step generated) radius(default = 4.0)"
			#for example "python step2_detectholetype.py 1PS6_chainA.pdb holes.pdb 4.0"
			
	Step 3: - Creat a directory named 'cna_results', and copy the '*local_indices.dat' file to this directory.
			#for example "mv 1PS6_chainA_local_indices.dat cna_results"
			
	Step 4: - run "python step4_holes_join_to_cna.py local_indices.dat(downloaded from CNA web server)"
			#for example "python step4_holes_join_to_cna.py ./cna_results/1PS6_chainA_local_indices.dat"
			
	Step 5: - run "python step5_preparepdb_rosetta.py pdbfilename nstruct(regulates the number of outputs per input structure) cpunumber"
			## Notes: if cpunumber > 1, the Rosetta mpi version must be installed correctly, or set cpunumber = 1.
			## Notes: larger nstruct number will spend more time.
			#for example "python step5_preparepdb_rosetta.py 1PS6_chainA.pdb 50 10"
			
	Step 6: - run "python step6_get_new_mutlist_withWT.py cpunumber(how many CPUs were used for Rosetta mpirun)"
			#for example "python step6_get_new_mutlist_withWT.py 10"
			
	Step 7: - run "python step7_mut_and_relax pdbfilename nstruct relaxrange(residues that how far from mutant residue will be used for relaxation) cpunumber"
			## Notes: if cpunumber > 1, the Rosetta mpi version must be installed correctly, or set cpunumber = 1.
			## Notes: larger nstruct and relaxrange spend more time.
			#for example "python step7_mut_and_relax pdbfilename 50 12 10"
			
	Last results: At last a file named "bestscoreandpdb.csv" and best-scored pdb will be saved to the directory named "selected_bestscorepdb". The "selected_bestscorepdb" can be found in the new directory named by input pdb file name. 
			
Statement:
	This workflow relies on academic free software such as ROSETTA and should not be used for commercial purposes without permission.
	
