
# tests
A = ExtMatrix(2*ones(2,2))
@argcheck (A[1,1],A[0,0],A[-1,0],A[3,1])==(2,1,1,1)