## checks
transfer_exceptions(1,1,1,2,2)
transfer_exceptions(0.9,0.9,0,0,1)
transfer_exceptions(0.5,0.5,0.5,1,0.5)
transfer_exceptions(1,1,1,1,0.5)
# test transfer: particular case of integrating out a corner and propagating empty diagonal
(a,b,c,d) = [rand() for _ in 1:4]
@argcheck collect(transfer_exceptions(1,1,1,a,b)[2]) ≈ [(a+b)/(1+a*b),1,1,1,1]
@argcheck collect(transfer_exceptions(1,a,b,c,d)[2]) ≈ [1,a,b,c,d]








# test YDelta 
J =Complex[2*rand() for _ in 1:3]
(R,K) = YDelta(J...)
ind = digits(rand(0:6),base=2,pad=3)
@argcheck Yval(J,ind) ≈ R*Dval(K,ind)

J = [2,2,2]
(R,K) = YDelta(J...)
ind = digits(rand(0:6),base=2,pad=3)
@argcheck Yval(J,ind) ≈ R*Dval(K,ind)


# test DeltaY
K = Complex[2*rand() for _ in 1:3]
(R,J) = DeltaY(K...)
ind = digits(rand(0:6),base=2,pad=3)
@argcheck R*Yval(J,ind) ≈ Dval(K,ind)

# test with number larger than 1

K = Complex[1.1,0.9,0.9]
(R,J) = DeltaY(K...)
ind = digits(rand(0:6),base=2,pad=3)
@argcheck R*Yval(J,ind) ≈ Dval(K,ind)

for i = 0:7
    ind = digits(i,base=2,pad=3)
    Base.display((R*Yval(J,ind),Dval(K,ind)))
end