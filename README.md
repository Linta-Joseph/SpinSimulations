# SpinSimulations
Matlab code for simulating spin systems. 

**DecouplingSequences_Magnusterms** - Numerically calculates Magnus expansion terms up to a set maximum order in the effective Hamiltonian engineered by various decoupling sequences. Main script to run: *MagnusTerms_main* 

**DecouplingSequences_Unifidelity** - Calculates fidelity of the experimental unitary for various dipolar decoupling sequences. Main script to run: *Calculate_Unifidelity*

**MBL_EE** - Sets and evolves a random product initial state of a spin Hamiltonian and calculates the bipartite entanglement entropy at each timestep. Main script to run: *ESS_and_EE*

**MBL_EES** - Diagonalizes a spin Hamiltonian and calculates the Schmidt decomposition, bipartitie entanglement entropy, and the entanglement spectrum statistics of mid-spectrum eigenstates for a grid of relevant Hamiltonian parameters. Main script to run: *EES_main*

**MBL_OTOC** - Time evolves a collective spin operator under a spin Hamiltonian and calculates the ZZ OTOC at each timestep. Main script to run: *Calculate_MQC_Widths*

**Parametres** - Calculates and plots the inter-pulse delays of the Hamiltonian engineering pulse sequence to achieve a certain set of Hamiltonian parameters experimentally. 
