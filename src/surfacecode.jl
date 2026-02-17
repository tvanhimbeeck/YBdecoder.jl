#= unrotated surcace code with top corner of the form
Q - X - (etc) 
|   |
Z - Q - (etc)
|   |
(etc) 
=#

function SCsyndromez(Xerror::Tuple{Matrix{Bool},Matrix{Bool}})::Matrix{Bool}
    (H,L) = size(Xerror[1])
    # compute syndrome of Z checks and return a H-1 x L matrix
    sz = [( 
        if y == 1
            xor(Xerror[1][x,y], Xerror[1][x+1,y], Xerror[2][x,y]) 
        elseif y == L
            xor(Xerror[1][x,y], Xerror[1][x+1,y], Xerror[2][x,y-1]) 
        else
            xor(Xerror[1][x,y], Xerror[1][x+1,y], Xerror[2][x,y], Xerror[2][x,y-1]) 
        end
        ) for x in 1:H-1, y in 1:L]
    return sz
end


function cosetrep(sz,x0)
    # coset representative in horizontal gauge with bottom right qubit in state x
    (H,L) = size(sz).+[1,0] 
    
    Xrep = ( zeros(Bool,H,L), Matrix{Bool}(undef,H-1,L-1) )
    #Xrep[1] = zeros(Int,H,L)
    #Xrep[2] = zeros(Int,H-1,L-1)
    for y = 1:L-1
        if y == 1
            Xrep[2][:,y] = sz[:,y]
        else
            Xrep[2][:,y] = map(xor, Xrep[2][:,y-1], sz[:,y])
        end
    end

    for x = H:-1:1
        if x == H
            Xrep[1][x,L] = x0
        else
            Xrep[1][x,L] = xor(Xrep[1][x+1,L], Xrep[2][x,L-1], sz[x,L] )
        end
    end
    return Xrep
end

function cosetrep_from_error(Xerror::Tuple{Matrix{Bool},Matrix{Bool}})::Tuple{Matrix{Bool},Matrix{Bool}}
    Xrep = Xerror
    (H,L) = size(Xerror[1])
    for y = 1:L-1
        for x = 1:H
            x0 = Xrep[1][x,y]
            Xrep[1][x,y] = 0

            Xrep[1][x,y+1] = xor(Xrep[1][x,y+1], x0)
            if x!= 1
                Xrep[2][x-1,y] = xor(Xrep[2][x-1,y], x0)
            end
            if x!=H
                Xrep[2][x,y] = xor(Xrep[2][x,y], x0)
            end
        end
    end
    return Xrep
end



@argcheck (L = 15; # number of qubits in the top row
H = 15; # number of qubits in the left collumn
px = 0.05;
q = px/(1-px);
Xerror = (rand(Bernoulli(px),H,L),rand(Bernoulli(px),H-1,L-1)) ;
nZerros = sum( sum(Xerror[i]) for i in 1:2);
sz = SCsyndromez(Xerror);
SCsyndromez(cosetrep(sz,0))==sz)

#mapping to bound state


function surfacecode_simulateerrorX(H::Int,L::Int,px::Number) 
    # generate H x L surface code with qubits on the corners and Z check at (2,1)
    # returns 
    #   H-1 x L array of Z syndromes
    #   (H x L) and (H-1 x L-1) arrays of qubits X errors

    Xerror = (rand(Bernoulli(px),H,L),rand(Bernoulli(px),H-1,L-1)) 
    #nXerrors = sum( sum(Xerror[i]) for i in 1:2)
    sz = SCsyndromez(Xerror)
    #@info "$H x $L surface code with $nXerrors X errors violating $(sum(sz)) Z checks"
    return (sz,Xerror)
end

function surfacecode_bondstateX(sz,px)
    # given a H-1 x L array of Z syndromes
    # generate an H x L+1 ising model 
    #       where the bottom corners are free 
    #       and the left and right borders are shorts
    (H,L) = size(sz).+[1,0]

    (Hbs,Lbs) = (H,L+1)

    
    q = ComplexF64(px/(1-px))
    Xrep = cosetrep(sz,0)
    bondh = map( x-> (x==false ? q : 1/q), Xrep[1])
    
    
    bondv = ComplexF64[ (j==1 || j==Lbs) ? zero(ComplexF64) : (Xrep[2][i,j-1] == true ? 1/q : q)
                for i in 1:Hbs-1, j in 1:Lbs
            ]   
    #= bondv = Matrix{ComplexF64}(undef,Hbs-1,Lbs)
    bondv[:,1] .= 0
    bondv[:,2:Lbs-1] = map(x->ComplexF64(q^(1 -2*x)),Xrep[2])#q.^(1 .- 2*Xrep[2])
    bondv[:,Lbs] .= 0 =#
    #bondv = [zeros(Hbs-1,1) q.^(1 .- 2*Xrep[2]) zeros(Hbs-1,1)]
    bondd = ones(ComplexF64,Hbs-1,Lbs-1)
    return Bondstate(Hbs,Lbs,bondh,bondv,bondd)
end

function surfacecode_decode(sz,px)
    bs = surfacecode_bondstateX(sz,px)
    (R,a) = fullsimplification2!(bs)
    coset = Int(abs(a)<=1)
    loglikelyhood = log(10,abs(a))
    return (coset,loglikelyhood)
end

function surfacecode_checksolution(Xerror,sz,coset)
    Xerrorrep = cosetrep_from_error(Xerror)
    Xrepsol = cosetrep(sz,1-coset)
    success = all( Xerrorrep[i]==Xrepsol[i] for i in 1:2)
    return success
end 
   

function surfacecode_simulate(H,L,px)
    (sz,Xerror) =  surfacecode_simulateerrorX(H,L,px)
    (coset,loglikelyhood)  = surfacecode_decode(sz,px)
    success = surfacecode_checksolution(Xerror,sz,coset)
    #@info "decoding successful : $success, in t=$(t) seconds, with loglikelyhood = $loglikelyhood"
    return success
end 

function surfacecode_simulate_many(H::Int,L::Int,px;maxtrials=1000,maxfailure=100)
    nfailure = 0
    ntrials = 0
    nabort = 0
    decodingsuccess = 0

#=     prog = ProgressUnknown(desc="Decoding: ntrials=")
    naborts = 0
    for n in 1:100
        naborts += (rand()>=0.5)
        next!(prog; showvalues = [("naborts",naborts)])
        if n==100
            finish!(prog)
            break
        end
        sleep(0.1)
    end =#


    prog = ProgressUnknown(desc="Simulating $(H)x$(L) surface code with p=$px",showspeed=true,spinner=true)
    while true
        ntrials += 1
        decodingsuccess = 0
        try 
            decodingsuccess = surfacecode_simulate(H,L,px)
            nfailure += (1-decodingsuccess)
        catch 
            nabort += 1
        end
        next!(prog; showvalues = [("n failures(/$maxfailure)",nfailure),("failure rate(/$maxfailure)",nfailure/ntrials),("ntrials(/$maxtrials)",ntrials)])
        if nfailure + nabort == maxfailure || ntrials == maxtrials
            finish!(prog)
            break
        end
    end
    failurerate = nfailure/ntrials
    abortrate = nabort/ntrials
    @info "L=$L, H=$H, p=$px, failurerate=$failurerate, abortrate=$abortrate, ntrials=$ntrials"
    return (failurerate,abortrate)
end 

function surfacecode_errorcurve(H::Int,L::Int,prange;maxtrials=1000,maxfailure=100)
    errorrate = zeros(length(prange))
    abortrate = zeros(length(prange))
    for i in eachindex(prange)
        #@info "p = $(prange[i])"
        (er,ar)  = surfacecode_simulate_many(H,L,prange[i];maxtrials=maxtrials,maxfailure=maxfailure)
        errorrate[i] = er
        abortrate[i] = ar
    end
    return (errorrate,abortrate)
end 

# proceeds in reverse order to avoid high suppression region
function surfacecode_errorcurve_fast(H::Int,L::Int,prange;maxtrials=1000,maxfailure=100)
    errorrate = zeros(length(prange))
    abortrate = zeros(length(prange))
    for i in length(prange):-1:1
        #@info "p = $(prange[i])"
        (er,ar)  = surfacecode_simulate_many(H,L,prange[i];maxtrials=maxtrials,maxfailure=maxfailure)
        errorrate[i] = er
        abortrate[i] = ar
        if er+ar<=sqrt(maxfailure)/maxtrials
            errorrate[1:i-1] .= NaN # nan to indicate unmeasured value
            abortrate[1:i-1] .= NaN           
            break
        end
    end
    return (errorrate,abortrate)
end 

function surfacecode_thresholdcurve(Lrange,prange;maxtrials=1000,maxfailure=100,reverse=false)
    errorrates = zeros(length(Lrange),length(prange))
    abortrates = zeros(length(Lrange),length(prange))
    for l in 1:length(Lrange)
        @info "L = $(Lrange[l])"
        if reverse
            (er,ar) = surfacecode_errorcurve_fast(Lrange[l],Lrange[l],prange,maxtrials=maxtrials,maxfailure=maxfailure)
        else
            (er,ar) = surfacecode_errorcurve(Lrange[l],Lrange[l],prange,maxtrials=maxtrials,maxfailure=maxfailure)
        end
        errorrates[l,:] = er
        abortrates[l,:] = ar
    end
    return (errorrates,abortrates)
end



# test 
