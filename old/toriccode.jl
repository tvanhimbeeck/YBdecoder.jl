using Distributions

# toric codes (on a torus)
L = 10 
H = 10 
px = 0.2
pz = 0.2

Zerror = (rand(Bernoulli(pz),H,L),rand(Bernoulli(pz),H,L))
nZerros = sum( sum(Zerror[i]) for i in 1:2)

Xerror = (rand(Bernoulli(px),H,L),rand(Bernoulli(px),H,L))
nZerros = sum( sum(Xerror[i]) for i in 1:2)

mymod(x,L) = (x+L-1)%L+1
sz = [ xor(Xerror[1][x,y], Xerror[1][mymod(x+1,H),y], Xerror[2][x,y], Xerror[2][x,mymod(y-1,L)]) for x in 1:H, y in 1:L]
sx = [ xor(Zerror[1][x,y], Zerror[1][x,mymod(y+1,H)], Zerror[2][x,y], Zerror[2][mymod(x-1,L),y]) for x in 1:H, y in 1:L]

#the lattice is regrouped in 2x2 in the configuration [Q1 X; Z Q2]
# sz= syndrome measured by Z checks (ie. about the Xerrors) etc.


# visualization
disp = Matrix{Any}(undef,0,2L)
for i = 1:H
    row = hcat([ [(Xerror[1][i,j]|>Int,Zerror[1][i,j]|>Int) ("sx = $(sx[i,j]|>Int)");
                  ("sz = $(sz[i,j]|>Int)") (Xerror[2][i,j]|>Int,Zerror[2][i,j]|>Int) ] for j in 1:L]...)
    disp = vcat(disp,row)
end
disp

