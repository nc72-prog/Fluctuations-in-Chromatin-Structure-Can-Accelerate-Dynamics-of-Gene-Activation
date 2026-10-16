# Fluctuations-in-Chromatin-Structure-Can-Accelerate-Dynamics-of-Gene-Activation
Monte Carlo Simulation and plotting code for PNAS publication. Code used to generate analytical results in figure 3 and for simulation results in figure 4. 

Authors: 
Natalie Collins 
Isabella Muradyan 
Anatoly Kolomeisky 





What the code does

The code computes the mean time a transcription factor (TF) needs to find a specific target site on DNA when the DNA is covered by chromatin that opens and closes dynamically, and compares it with the search time on chromatin-free DNA.

The DNA is modeled as a lattice of L sites with the target at site m. A TF in solution binds to accessible DNA, slides along it by one-dimensional diffusion (rate u), and dissociates back into solution (rate koff). Chromatin uncovers DNA segments of l sites with rate α/l and covers them again with rate β; a TF bound to DNA when chromatin closes is ejected into solution. Every search starts with the DNA fully covered and the TF in solution.

The repository covers three cases:

Intermediate opening (1 < l < L): chromatin uncovers a segment of l sites at a time (l = 10 in the simulations). Mean search times are obtained from exact stochastic (Gillespie) simulations.
Single-site opening (l = 1): chromatin uncovers one site at a time. Mean search times T_00 are calculated analytically from the backward master equation solutions derived in the paper.
Whole-DNA opening (l = L): chromatin uncovers the entire DNA at once. Mean search times T_00 are calculated analytically from the backward master equation solutions derived in the paper.

The sliding length is λ = √(u/k_off); the heatmap scripts set k_off = u/λ² for each panel.

In each case, the results are compared with the analytical search time on chromatin-free DNA, T_free (Shvets, Kochugaeva & Kolomeisky, Molecules 23, 2106 (2018)). The ratio T_free / T_00 shows where chromatin accelerates the search (ratio > 1, green in the heatmaps) or slows it down (ratio < 1, red)

Repository contents
File	Description	Output
figure3a.py : Analytical T_00 versus β for single-site opening (l = 1), for several α.	
figure3b.py :	Analytical T_00 versus closing rate β for whole-DNA opening (l = L), for several opening rates α.
L=1heatmap.py Heatmaps of T_free / T_00 over (α, β) for l = 1, at sliding lengths λ = 10, 100, 1000.
L=lheatmap.py	Heatmaps of T_free / T_00 over (α, β) for l = L, at λ = 10, 100, 1000.	
l=ngillespie.cpp	Gillespie simulation of the search with segment opening (l = 10) over a grid of (α, β).	
l=ngillespie.py	Plots T_free / T_chrom from the simulation output.	



Requirements & Operating system

The code was developed and run on Red Hat Enterprise Linux release 8.10 (Ootpa). It uses only standard Python packages and standard C++, so it also runs on macOS and Windows without modification. No GPU or special hardware is required.

Software and packages
Software	Version used	Purpose
Python	3.6.8	Analytical calculations and plotting
NumPy	1.19.5	Numerical calculations
Matplotlib	3.3.4	Figures
g++ (GCC) 8.5.0 20210514 (Red Hat 8.5.0-29)	Gillespie simulations


Running the analysis

Run all commands from the repository folder. Each Python script saves its figure as a PNG in the same folder and opens it in a window. Parameters are set at the top of each script and can be changed there.

Figure 3: analytical results
bash
python3 figure3a.py        # T_00 vs beta, single-site opening (l = 1)
python3 figure3b.py        # T_00 vs beta, whole-DNA opening (l = L)


python3 L=1heatmap.py      # T_free / T_00 heatmaps, l = 1 (will create the top panel of figure 4)
python3 L=lheatmap.py      # T_free / T_00 heatmaps, l = L (will create the bottom panel of figure 4) 


Gillespie simulations for intermediate segment opening. You will need to change the koff parameter in the .cpp file and in the .py file to create a new image for each middle panel. Likewise, you may want to change the name of the .dat file for each simulation. 

Compile and run the simulation:

bash
g++ -O2 -o l=ngillespie.cpp -o l=ngillespie
./l=ngillespie


Then plot the results:

bash
python3 l=ngillespie.py


Parameters are set at the top of each script and can be changed there.

Citation
If you use this code, please provide proper citation. 

[See MIT License file]

Contact
Anatoly B. Kolomeisky, Rice University — tolya@rice.edu
