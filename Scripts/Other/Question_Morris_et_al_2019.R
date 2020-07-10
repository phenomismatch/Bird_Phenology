#adjacency matrix - given in Morris et al. 2019 Section 2
A <- matrix(c(0,1,0,0,1,0,1,0,0,1,0,1,0,0,1,0), 
            nrow = 4, byrow = TRUE)

#row-standardized (scaled) A
B <- matrix(c(0,1,0,0,0.5,0,0.5,0,0,0.5,0,0.5,0,0,1,0), 
            nrow = 4, byrow = TRUE)

#diagonals are number of neighbors for each node
D <- matrix(c(1,0,0,0,0,2,0,0,0,0,2,0,0,0,0,1), 
            nrow = 4, byrow = FALSE)

#identity matrix
I <- diag(1, nrow = 4)

#set to 1 as with intrinsic CAR
alpha <- 1

#Given in section 2.1 of Morris et al. 2019
#Q = D(I - alpha * A)
#Should be B [scaled adjacency matrix]?
D %*% (I - (alpha * A))
#Morris Stan case study
#Q = D - alpha * A
D - alpha * A
#CAR Stan case study
#Q = D(I - alpha * B)
D %*% (I - (alpha * B))
#Q = D - alpha * A - same as Morris Stan case study
#Ver Hoef et al. 2018 Ecological Monographs Eq. 15
#Q = D - alpha * A - same as Morris Stan case study and CAR Stan case study


#Given in section 2.2 of Morris et al. 2019:
#D(I − alpha * A) simplifies to D − A
all.equal(D %*% (I - (alpha * A)), D - A)
#I think this should be B (scaled adjacency matrix) instead of A (adjacency matrix)
all.equal(D %*% (I - (alpha * B)), D - A)

#inverse of D multiplied by A = B
solve(D) %*% A