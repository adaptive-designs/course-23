install.packages("librarian")
librarian::shelf(clinfun, curtailment)


# 1:
ph2single(pu = 0.05, pa = 0.3, ep1 = 0.05, ep2 = 1-0.8, nsoln = 5)


# 2:
ph2simon(pu = 0.05, pa = 0.3, ep1 = 0.05, ep2 = 1-0.8)

# Personal opinion.
# Single-stage design is the most simple to run.
# Minimax has the same N_max as the single-stage, but lower ESS(H0).
# Optimal has lowest ESS(H0) (7.9, vs 9.1 for minimax), but greater N-max (18 vs 14).
# Admissible has a low ESS(H0) (8.4) with only a slight increase in N_max (15 vs 14).


# 3:
# For each design, the values qLo and qHi denote the range of weights for N_max
# for which the trial is optimal. Allocating equal weighting to N_max and ESS(H0)
# means giving both optimality criteria a weight of 0.5.

# As such, the most appropriate design is the design that contains 0.5 in the 
# [qLo, qHi] range. In this case, this is the minimax design, which is the most
# appropriate design for N_max weights in the range [0.421, 1.000]. The parameters
# of this design are

# r1=0, n1=7, r=2, n=14.

# The probability that the trial will end early if the response rate is equal
# to 0.05 is PET(p0)=0.6983.


# 4:
response.lo <- ph2simon(pu = 0.05, pa = 0.2, ep1 = 0.05, ep2 = 1-0.8)
response.lo

# For p1=0.2, the most appropriate design for the weight (0.5) is the admissible design.
# This design has a greater N_max and ESS(H0) compared to earlier designs.

response.hi <- ph2simon(pu = 0.05, pa = 0.4, ep1 = 0.05, ep2 = 1-0.8)
response.hi

# For p1=0.4, the most appropriate design for the weight (0.5) is the minimax design.
# This design has a lower N_max and ESS(H0) compared to earlier designs.


# 5:
mander <- find2stageDesigns(nmin=6,
                            nmax=25,
                            p0=0.05,
                            p1=0.3,
                            alpha=0.05,
                            power=0.8,
                            benefit = T)
mander