using SparseArrays, ArgCheck
include("star-triangle.jl")
# code for computing the partition function of the Ising model on a H by L grid with given bonds
L = 4
H = 3

# grid
nV = H*L
nEh = (L-1)*H
nEv = L*(H-1)
nF = (L-1)*(H-1)
nbonds = nEh + nEv + nF

# problem variables
sv = zeros(Bool,nEv)
sh = zeros(Bool,nEh)
q = 0.5

# function for descibing a square grid
Ev(i) = (@argcheck 0<=i<=nEv-1; (i÷L, i%L))
Ev(x,y) = (@argcheck (0<=x<=H-2) && (0<=y<=L-1); x*L+y)
Eh(i) = (@argcheck 0<=i<=nEh-1; (i÷(L-1),i%(L-1)))
Eh(x,y) = (@argcheck (0<=x<=H-1) && (0<=y<=(L-2)); x*(L-1)+y)
Ed(i) = (@argcheck 0<=i<=nF-1; (i÷(L-1),i%(L-1)))
Ed(x,y) = (@argcheck ((0<=x<=H-2) && (0<=y<=L-2)); x*(L-1)+y)

# test 
@argcheck (i = rand(0:nEv-1); Ev(Ev(i)...)==i)
@argcheck ((x,y) = (rand(0:H-2),rand(0:L-1)); Ev(Ev(x,y))==(x,y))

function initialize(q)
    # initialize variables
    bondh = [q for i in 1:nEh]
    bondv = [q for i in 1:nEv]
    bondd = Float64[1 for i in 1:nF]
    return (bondh,bondv,bondd)
end


function integratecorner!(x,y,bonds,R)
    (bondh,bondv,bondd) = bonds
    @argcheck 0<=x<H-1 && 0<=y<L-1
    @argcheck y==0||bondh[Eh(x,y-1)+1]==1
    @argcheck x==0||bondv[Ev(x-1,y)+1]==1
    
    v = bondh[Eh(x,y)+1]
    w = bondv[Ev(x,y)+1]
    @info (v,w) 
    bondh[Eh(x,y)+1] = 1.0
    bondv[Ev(x,y)+1] = 1.0
    bondd[Ed(x,y)+1] *= (v+w)/(1+v*w)
    R *= (1+v*w)
end

function diagpropagation!(x,y,bonds,R)
    
    @argcheck 0<=x<H-2 && 0<=y<L-2
    (bondh,bondv,bondd) = bonds
    v0 = bondd[Ed(x,y)+1]
    (a,b) = (bondh[Eh(x+1,y)+1],bondv[Ev(x,y+1)+1])
    (c,d) = (bondv[Eh(x+1,y+1)+1],bondh[Eh(x+1,y+1)+1])
    (R1,(v1,B,A)) = DeltaY(v0,a,b)
    (R2,(v2,D,C)) = YDelta(v1,c,d)

    R *= R1*R2
    bondd[Ed(x,y)+1] = 1
    bondd[Ed(x+1,y+1)+1] = v2
    bondh[Eh(x+1,y)+1] = A
    bondv[Ev(x,y+1)+1] = B
    bondv[Ev(x+1,y+1)+1] = C
    bondh[Eh(x+1,y+1)+1] = D
end

function diagterminination(x,y,bonds)
    @argcheck 0<=x<H-1 && 0<=y<L-1
end

function represent(bonds)
    (bondh,bondv,bondd) = bonds
    M = spzeros(Number,2H-1,2L-1)
    for i = 0:H-1, j = 0:L-1
        M[2i+1,2j+1] = 0
    end
    for i = 0:H-1, j = 0:L-2
        M[2i+1,2j+2] = bondh[Eh(i,j)+1]
    end
    for i = 0:H-2, j = 0:L-1
        M[2i+2,2j+1] = bondv[Ev(i,j)+1]
    end
    for i = 0:H-2, j = 0:L-2
        M[2i+2,2j+2] = bondd[Ed(i,j)+1]
    end
    display(M)
end

bonds = initialize(0.5)
represent(bonds)
integratecorner!(0,0,bonds,R)
represent(bonds)
diagpropagation!(0,0,bonds,R)
represent(bonds)