# provides star-triangle relationships in the xor representation


Yval(J,ind) = sum( prod( J[i]^xor(s,ind[i]) for i = 1:3) for s = 0:1)
Dval(K,ind) = prod( K[(i)+1]^xor(ind[(i+1)%3+1],ind[(i+2)%3+1]) for i = 0:2)
function YDelta(j1::Number,j2::Number,j3::Number) 
    @argcheck !isnan(j1) && !isnan(j2) && !isnan(j3) "Nan numbers in YDelta : $((j1,j2,j3))"
    


    if j1==0                        # j1 is special
        return (1,(1,j3,j2))
    elseif abs(j1)==Base.Inf
        return (1,(1,1/j3,1/j2))
    elseif j2 == 0                  # j2 is special 
        return (1,(j3,1,j2))
    elseif abs(j2) == Base.Inf
        return (1,(1/j3,1,1/j2))
    elseif j3 == 0                  # j3 is special 
        return (1,(j2,j1,1))
    elseif abs(j3) == Base.Inf
        return (1,(1/j2,1/j1,1)) 
    end
    z0 = 1 + j1*j2*j3
    z1 = j1 + j2*j3 
    z2 = j2 + j1*j3 
    z3 = j3 + j1*j2 

    k = sqrt(z1*z2*z3/z0)
    k1 = k/z1
    k2 = k/z2
    k3 = k/z3
    return (1,(k1,k2,k3))
end

#dual(x::Number) = (1-x)/(1+x)
function dual(x::T)::T where {T<:Number}
    if abs(x) == Base.Inf
        return -1
    elseif x == -1
        return Base.Inf
    end
    return (1-x)/(1+x)
end

function DeltaY(k1::Number,k2::Number,k3::Number)
    
    K = (k1,k2,k3)
    Kprim = map(dual,K)
    
    (_,Jprim) = YDelta( Kprim...)
    J = map(dual,Jprim)
    
    #ind = (0,0,0)
    #R = Dval(K,ind)/Yval(J,ind)
  
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
    (R1,(v1,B,A)) = DeltaY(v0,a,b)
    (R2,(v2,D,C)) = YDelta(v1,c,d)
    R = R1*R2
    @argcheck !any(isnan.((v2,A,B,C,D))) "NaN error when applying diagtransfer from $((v0,a,b,c,d)) to $((v2,A,B,C,D))"
    return (R,(v2,A,B,C,D))
end
function transfer_exceptions(v0::T,a::T,b::T,c::T,d::T)::Tuple{T,NTuple{5,T}} where {T<:Number}
        
    #exceptions
    if (v0,a,b,c)==(one(T),one(T),one(T),one(T))
        # this is the left bottom corner        -> leave as is 
        return (one(T),(one(T),one(T),one(T),one(T),d))
    elseif (v0,b,c,d)==(1,1,1,1)    
        # bottom right corner without diagonal  -> leave as is
        return (one(T),(one(T),a,one(T),one(T),one(T)))
    elseif (v0,a,b)==(1,1,1)        
        # corner that is integrated out
        return ((one(T)+c*d)/2, ((c+d)/(one(T)+c*d),one(T),one(T),one(T),one(T)))   ###### divied by 2 here #############
    elseif (b,c,d)==(0,0,1) || (b,c,d)==(0,1,1)      
        # diagonal arriving on right border (or bottom right corner) -> fuse the diagonal and horizontal edges
        return (1,(1,v0*a,b,c,d))
    elseif (v0,b,c)==(1,1,1)            
        # --o-- to ==o-- on the bottom border
        return (1+a*d,(one(T),zero(T),one(T),one(T),(a+d)/(one(T)+a*d)))
    elseif c==1
        #diagonal arriving on the bottom border -> choose solution where the bottom edge is positive and equal to 1
        (R,(v2,A,B,C,D)) = transfer(v0,a,b,c,d)
        if C≈-1
            return (R,(-v2,A,B,-C,-D)) 
        else 
            return (R,(v2,A,B,C,D))
        end
    else
        return transfer(v0,a,b,c,d)
    end
end


