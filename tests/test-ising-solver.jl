include("..//ising-solver.jl")

#test
bs = initialize(2,2,0.1)


# test display
display(bs)
displaysimp(bs)

# test on a lower triangle
bs = triangle(0.5)
exactsummation(bs)
fullsimplification!(bs)
exactsummation(bs)
bs.R

diagpropagation!(2,2,bs)
display(bs)
exactsummation(bs)
diagpropagation!(1,2,bs)
display(bs)
exactsummation(bs)
diagpropagation!(2,1,bs)
display(bs)
exactsummation(bs)
diagpropagation!(2,2,bs)
display(bs)
exactsummation(bs)
bs.R
bs = triangle(0.5)
exactsummation(bs)

bs = triangle(0.5)


# compare exact summation and bondpropagation
bs = initialize(2,2,0.5)
r1 = exactsummation(bs)
r2 = fullsimplification!(bs)
log(2,r2/r1)

bs = initialize(2,4,0.5)
r1 = exactsummation(bs)
r2 = fullsimplification!(bs)
log(2,r2/r1)

bs = initialize(3,4,0.4)
r1 = exactsummation(bs)
r2 = fullsimplification!(bs)
r3 = exactsummation(bs)
log(2,r2/r1)
log(2,r3/r1)

# check ising on a line
bs = initialize(1,2,0.5)
r1 = exactsummation(bs)
r2 = fullsimplification!(bs)

