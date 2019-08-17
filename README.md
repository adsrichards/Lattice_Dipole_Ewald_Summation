# Ewald Summation for Lattice Dipoles

The lattice sum for an infinite arrangement of dipoles in 3 dimensions is conditionally convergent and can't be approximated by truncation. Ewald summation transforms this sum into two fast converging sums.

We may simulate an infinite lattice system by imposing unlimited periodic boundary conditions on a finite lattice (think of a room with mirrors for walls where your image is replicated infinitely). Since the coupling strength of the infinite sum has no dependence on the orientation of dipoles, the coupling strength of the finite system may be computed once and stored in memory for use in Monte Carlo simulations.
