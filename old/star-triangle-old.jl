# provides star-triangle relationships

function YDelta(j1::Number,j2::Number,j3::Number) 
    z0 = 1 + j1*j2*j3
    z1 = j1 + j2*j3
    z2 = j2 + j1*j3
    z3 = j3 + j1*j2
    k0 = (z0*z1*z2*z3)^(1/4)/(j1*j2*j3)^(1/2)
    k1 = (z0*z1/(z2*z3))^(1/2);
    k2 = (z0*z2/(z1*z3))^(1/2);
    k3 = (z0*z3/(z1*z2))^(1/2);
    return (k0,(k1,k2,k3))
end

dual(x::Number) = (1-x)/(1+x)

function DeltaY(k1::Number,k2::Number,k3::Number)
    K = (k1,k2,k3)
    Kprim = map(dual,K)
    @info Kprim
    Jprim = YDelta( Kprim...)[2]
    @info Jprim
    J = map(x->-dual(x),Jprim)
    ind = (1,1,1)
    R = Dval(K,ind)/Yval(J,ind)
    return (R,J)
end

# consistency checks
yYval(J,ind) = sum( prod( J[i]^(s*ind[i]/2) for i = 1:3) for s = -1:2:1)
Dval(K,ind) = prod( K[(i)+1]^(ind[(i+1)%3+1]*ind[(i+2)%3+1]/2) for i = 0:2)

# test YDelta 
J = [rand() for _ in 1:3]
(R,K) = YDelta(J...)
ind = (-1).^digits(rand(0:7),base=2,pad=3)
@argcheck Yval(J,ind) ≈ R*Dval(K,ind)

# test DeltaY
K = [rand() for _ in 1:3]
(R,J) = DeltaY(K...)
ind = (-1).^digits(rand(0:7),base=2,pad=3)
ind = [1,1,1]
@argcheck R*Yval(J,ind) ≈ Dval(K,ind)
