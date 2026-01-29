# provides star-triangle relationships in the xor representation

using ArgCheck

Yval(J,ind) = sum( prod( J[i]^xor(s,ind[i]) for i = 1:3) for s = 0:1)
Dval(K,ind) = prod( K[(i)+1]^xor(ind[(i+1)%3+1],ind[(i+2)%3+1]) for i = 0:2)


function YDelta(j1::Number,j2::Number,j3::Number) 
    @argcheck j1>=0 && j2>=0 && j3>=0 "negative/Nan numbers in YDelta : $((j1,j2,j3))"
    
    z0 = 1 + j1*j2*j3
    z1 = j1 + j2*j3
    z2 = j2 + j1*j3
    z3 = j3 + j1*j2
    k1 = sqrt(z2*z3/(z1*z0))
    k2 = sqrt(z1*z3/(z2*z0))
    k3 = sqrt(z1*z2/(z3*z0))
    return (z0,(k1,k2,k3))
end

dual(x::Number) = (1-x)/(1+x)

function DeltaY(k1::Number,k2::Number,k3::Number)
    if isequal((k1,k2,k3),(1,1,1))
        return (1,(1,1,1))
    end
    if isequal((k1,k3),(1,1))
        return (1,(k2,1,1))
    end
    if isequal((k1,k2),(1,1))
        return (1,(k2,1,1))
    end
    K = (k1,k2,k3)
    Kprim = map(dual,K)
    Jprim = YDelta( Kprim...)[2]
    J = map(x->dual(x),Jprim)
    ind = (1,1,1)
    
    R = Dval(K,ind)/Yval(J,ind)
    return (R,J)
end

#=    o               o
    / |               |
  v0  b               B
 /    |               |
o--a--o--d-- => o--A--o--D--o
      |               |    /
      c               C  v2
      |               | /
                      0
=#
function transfer(v0,a,b,c,d)
    
    @argcheck (v0==1 ? true : (a!=1 || b!=1) )h    

    (R1,(v1,B,A)) = DeltaY(v0,a,b)
    (R2,(v2,D,C)) = YDelta(v1,c,d)
    R = R1*R2
    return (R,(v2,A,B,C,D))
end


## checks

# test transfer: particular case of integrating out a corner and propagating empty diagonal
(a,b,c,d) = [rand() for _ in 1:4]
@argcheck collect(transfer(1,1,1,a,b)[2]) ≈ [(a+b)/(1+a*b),1,1,1,1]
@argcheck collect(transfer(1,a,b,c,d)[2]) ≈ [1,a,b,c,d]


# test YDelta 
J = [rand() for _ in 1:3]
(R,K) = YDelta(J...)
ind = digits(rand(0:7),base=2,pad=3)
@argcheck Yval(J,ind) ≈ R*Dval(K,ind)

# test DeltaY
K = [rand() for _ in 1:3]
(R,J) = DeltaY(K...)
ind = digits(rand(0:7),base=2,pad=3)
@argcheck R*Yval(J,ind) ≈ Dval(K,ind)
