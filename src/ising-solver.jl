

## initalizations of ising models
function initialize(H::Int,L::Int,q::T) where {T<:Number}
    # initialize variables
    bondh = q*ones(T,H,L-1)
    bondv = q*ones(T,H-1,L)
    bondd = ones(T,H-1,L-1)
    return Bondstate(H,L,bondh,bondv,bondd)
end

function triangle(q)
    H=2
    L=2
    h = [1;q;;]
    v = [1 q]
    d = [q;;]
    return Bondstate(H,L,h,v,d)
end



## displays
function display(bs::Bondstate{T}) where {T}
    
    M = spzeros(T,2bs.H-1,2bs.L-1)
    for i = 1:bs.H, j = 1:bs.L
        M[2i-1,2j-1] = 0
    end
    for i = 1:bs.H, j = 1:bs.L-1
        M[2i-1,2j] = bs.h[i,j]
    end
    for i = 1:bs.H-1, j = 1:bs.L
        M[2i,2j-1] = bs.v[i,j]
    end
    for i = 1:bs.H-1, j = 1:bs.L-1
        M[2i,2j] = bs.d[i,j]
    end
    Base.display(M)
end

function displaysimp(bs::Bondstate)
    
    M = ["  " for i=1:2bs.H-1, j in 1:2bs.L-1]
    for i = 1:bs.H, j = 1:bs.L
        M[2i-1,2j-1] = "o"
    end
    for i = 1:bs.H, j = 1:bs.L-1
        if bs.h[i,j] == 1
            M[2i-1,2j] = " "
        elseif bs.h[i,j] == 0
            M[2i-1,2j] = "=="
        else
            M[2i-1,2j] = "--"
        end
    end
    for i = 1:bs.H-1, j = 1:bs.L
        if bs.v[i,j] == 1 
            M[2i,2j-1] = " "
        elseif bs.v[i,j] == 0
            M[2i,2j-1] = "||"
        else
            M[2i,2j-1] = "|" 
        end
    end
    for i = 1:bs.H-1, j = 1:bs.L-1
        M[2i,2j] = bs.d[i,j]!=1 ? "/" : " "
    end
    Base.display(M)
end



## bondpropagation algorithm
function diagpropagation!(x::Int,y::Int,bs::Bondstate{T}) where {T<:Number}

    @argcheck (0 <= x <= bs.H) && (0<= y <= bs.L)
    @argcheck bs.d[x-1,y]==1 && bs.d[x,y-1]==1

    
    @info "process vertex $((x,y))"

    #R1 = exactsummation_fixedcorners(bs,0,0)

    v0 = bs.d[x-1,y-1]
    (a,b) = (bs.h[x,y-1],bs.v[x-1,y])
    (c,d) = (bs.v[x,y],bs.h[x,y])
    
    @argcheck map(typeof,(v0,a,b,c,d))==(T,T,T,T,T)

    input = (v0,a,b,c,d)
    (R,ans) = transfer_exceptions(input...)
    (v2,A,B,C,D) = ans
    @argcheck !any(isnan.((v2,A,B,C,D))) == true "NaN error when applying diagtransfer $((v0,a,b,c,d)) to $((v2,A,B,C,D)) at vertex $((x,y))"

    bs.R *= R
    bs.d[x-1,y-1] = 1
    bs.d[x,y] *= v2
    bs.h[x,y-1] = A
    bs.v[x-1,y] = B
    bs.v[x,y] = C
    bs.h[x,y] = D
    
    #R2 = exactsummation_fixedcorners(bs,0,0)
#=     if abs(R1/R2-1)>0.001
        @info "non exact bond propagation R1 = $(R1) -> R2 = $(R2) at vertex $((x,y))"
        @info "with parameters $((v0,a,b,c,d)) to $((v2,A,B,C,D))"
    end =#
    return bs.R
end

function diagonalsimplification!(bs::Bondstate,delta1::Int,delta2::Int)
    
    x = Int((delta2 + delta1)/2)
    y = Int((delta2 - delta1)/2)
    @argcheck (1<=x<=bs.H) && (1<=y<=bs.L) "$((x,y,delta1,delta2))"

    while x<=bs.H && y<=bs.L
        #@info "process diagonal $(x+y)"
        diagpropagation!(x,y,bs)
        #display(bs)
        x+=1
        y+=1
    end
end

function fullsimplification!(bs)
    for delta2 = 2:bs.L+bs.H
        for delta1 = max(2-delta2,delta2-2*bs.L):2:min(2*bs.H-delta2,delta2-2)
            diagonalsimplification!(bs,delta1,delta2)
        end
    end
    #@argcheck all(bs.h[1:H-1,:] .==1) && all(bs.v .==1) && all(bs.d .==1) && all((bs.h[H,:])[1:L-2].==0)
    #@argcheck isapprox(bs.R,real(bs.R);rtol = 1e-8 )
    #@argcheck isapprox(bs.h[bs.H,bs.L-1],real(bs.h[bs.H,bs.L-1]); rtol=1e-8)
    return (real(bs.R),real(bs.h[bs.H,bs.L-1]))
end
        
function exactsummation(bs)
    sum = 0.0
    for i in 0:2^(bs.H*bs.L)-1
        I = digits(i,base=2,pad=bs.H*bs.L)
        I = reshape(I,bs.H,bs.L)
        product = 1
        product *= bs.H==1 ? 1 : prod( bs.v[x,y]^xor(I[x,y],I[x+1,y]) for x in 1:bs.H-1, y in 1:bs.L)
        product *= bs.L==1 ? 1 : prod( bs.h[x,y]^xor(I[x,y],I[x,y+1]) for x in 1:bs.H, y in 1:bs.L-1)
        product *= bs.H==1||bs.L==1 ? 1 : prod( bs.d[x,y]^xor(I[x+1,y],I[x,y+1]) for x in 1:bs.H-1, y in 1:bs.L-1)
        #@info (i,I,product)
        sum += product
    end
    @argcheck isapprox(bs.R*sum, real(bs.R*sum);rtol=1e-8)
    return real(bs.R*sum)
end

function exactsummation_fixedcorners(bs,u,v)
    sum = 0.0
    for i in 0:2^(bs.H*bs.L)-1
        I = digits(i,base=2,pad=bs.H*bs.L)
        I = reshape(I,bs.H,bs.L)
        if I[bs.H,1]!=u || I[bs.H,bs.L]!=v
            continue
        end
        product = 1
        product *= bs.H==1 ? 1 : prod( bs.v[x,y]^xor(I[x,y],I[x+1,y]) for x in 1:bs.H-1, y in 1:bs.L)
        product *= bs.L==1 ? 1 : prod( bs.h[x,y]^xor(I[x,y],I[x,y+1]) for x in 1:bs.H, y in 1:bs.L-1)
        product *= bs.H==1||bs.L==1 ? 1 : prod( bs.d[x,y]^xor(I[x+1,y],I[x,y+1]) for x in 1:bs.H-1, y in 1:bs.L-1)
        #@info (i,I,product)
        sum += product
    end    
    @argcheck isapprox(bs.R*sum, real(bs.R*sum);rtol=1e-8)
    return real(bs.R*sum)
end

function exactsummation_with_BC(bs,u,v) # exact summation with boundary conditions where left and right boundaries are ferromagnetic shorts and no diagonal coupling
      
    sum = 0.0
    for i in 0:2^(bs.H*(bs.L-2))-1
        I = digits(i,base=2,pad=bs.H*(bs.L-2))
        I = reshape(I,bs.H,bs.L-2)
        product = 1
        product *= prod( bs.v[x,y+1]^xor(I[x,y],I[x+1,y]) for x in 1:bs.H-1, y in 1:bs.L-2) # vertical edges first and last collumns have a trivial product

        product *= prod( bs.h[x,y+1]^xor(I[x,y],I[x,y+1]) for x in 1:bs.H, y in 1:bs.L-3)   # horizontal edges middle section
        product *= prod( bs.h[x,1]^xor(u,I[x,1]) for x in 1:bs.H)                           # horizontal edges first collumn
        product *= prod( bs.h[x,bs.L-1]^xor(v,I[x,bs.L-2]) for x in 1:bs.H)                 # horizontal edges last collumn

#        product *= prod( bs.d[x,y]^xor(I[x+1,y],I[x,y+1]) for x in 1:bs.H-1, y in 1:bs.L-1)
        #@info (i,I,product)
        sum += product
    end
    @argcheck isapprox(bs.R*sum, real(bs.R*sum);rtol=1e-8)
    return real(bs.R*sum)
end

#= function integratecorner!(x,y,bonds,R)
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
end =#