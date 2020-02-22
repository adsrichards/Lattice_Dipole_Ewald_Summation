# Ewald Summation for Lattice Dipoles

The lattice sum for an infinite arrangement of dipoles in 3 dimensions is conditionally convergent and shouldn't be approximated by truncation. Ewald summation transforms this sum into two fast converging sums.

We may simulate an infinite lattice system by imposing unlimited periodic boundary conditions on a finite lattice (think of a room with mirrors for walls where your image is replicated infinitely). Since the coupling strength of the infinite sum has no dependence on the orientation of dipoles, the coupling strength of the finite system may be computed once and stored in memory for use in Monte Carlo simulations.

My code calculates the coupling matrix using the Ewald summation method and runs Monte Carlo simulations for classical (Heisenberg) point dipoles on the simple cubic lattice at various temperatures in parallel. Results from an older version of this code are featured in a publication by Bovo et. al. on special tempuratures of magnetic systems with competing ferromagnetic and antiferromagnetic order.

[Nature Communications volume 9, Article number: 1999 (2018)](https://doi.org/10.1038/s41467-018-04297-3)

![](EwaldSum_Summarized-1.png)
![](EwaldSum_Summarized-2.png)
