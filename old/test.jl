include(pwd()*"\\LDPCcodes.jl")

# Surface code parameters qubits on the corners
H = 5
L = 7
nqubits = (H*L-1)/2

# generate error y
ex = [rand(px) in 1:nqubits]


syndromez = [rand(px) for i in 1:nqubits] 

#coordinate of points on a grid assuming H,L are odd and coordinates starting at 0
coord00(H,L,i) = (ptperline=(L+1)÷2; (2*(i÷ptperline),2*(i%ptperline)))
coord01(H,L,i) = (ptperline=(L-1)÷2; (1+2*(i÷ptperline),2*(i%ptperline)))
coord10(H,L,i) = (ptperline=(L+1)÷2; 1+(2*(i÷ptperline),2*(i%ptperline)))
coord01(H,L,i) = (ptperline=(L-1)÷2; 1+(2*(i÷ptperline),1+2*(i%ptperline)))

#surface code with qubits on the border
Hrep(n) = [(i==j||i+1==j) ? 1 : 0 for i in 1:n-1, j in 1:n]


(QZ,QX) = (3,4)

HZ = [kron(Hrep(QZ),Matrix(I,QX,QX)) kron(Matrix(I,QZ-1,QZ-1),Hrep(QX)-1)]|> sparse

rand(0.1)

surfacecode = HGproduct(Hrep(n),Hrep(n))
