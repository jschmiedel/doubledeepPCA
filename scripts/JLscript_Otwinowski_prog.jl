## Julia script to run Otwinowski deltaG method programmatically

using ArgParse
s = ArgParseSettings()
@add_arg_table s begin
    "--dataset", "-d"
        help = "dataset file"
    "--name", "-n"
        help = "dataset name"
    # "--flag1"
        # help = "an option without argument, i.e. a flag"
        # action = :store_true
    # "arg1"
        # help = "a positional argument"
        # required = true
end

parsed_args = parse_args(ARGS, s)

using CSV
using DataFrames
using SparseArrays
using DelimitedFiles
using NLopt
import StatsBase.predict
using Distributions
using KahanSummation
using LinearAlgebra

## define functions
function prepdata(df, seqname, kind, wt, c0n, c1n; delim = '-')
    n = nrow(df)
    I = collect(1:n)
    J = ones(Int64, n)
    j = 2
    aa = Char[' ']
    pos = Int64[0]
    code = Dict{String, Int64}()
    if kind == :listofmuts
        muts = df[!, seqname]
        for i = 1:n
            if muts[i] != wt
                for k in split(muts[i], delim)
                    if !haskey(code, k)
                        code[k] = j
                        push!(aa, k[end])
                        # push!(pos, parse(k[2:end-1])) ##issue here?
                        push!(pos, parse(Int64,k[2:end-1])) ## this has to be parsed to Int64
                        j += 1
                    end
                    push!(J, code[k])
                    push!(I, i)
                end
            end
        end
    elseif kind == :sequence
        seq = df[seqname]
        ns = length(wt)
        for i = 1:n
            for p = 1:ns
				a = seq[i][p]
	            if a != wt[p]
                    k = string(wt[p], p, a)
                    if !haskey(code, k)
                        code[k] = j
                        push!(aa, a)
                        push!(pos, p)
                        j += 1
                    end
                    push!(J, code[k])
                    push!(I, i)
                end
            end
        end
    end
    x = sparse(I, J, 1.0)
    g = trues(size(x, 2))
    g[1] = false
    # ham = ceil.(Int64, vec(sum(x,2))-1)
    ham = ceil.(Int64, vec(sum(x, dims = 2)) .- 1)
	
	c0 = collect(Missings.skipmissing(df[!, c0n]))
    c1 = collect(Missings.skipmissing(df[!, c1n]))
	# wti = find(df[seqname] .== wt)[1]
	wti = findall(df[!, seqname] .== wt)[1]
    f = log.((c1 .+ .5) ./ (c0 .+ .5) .* (c0[wti] .+ 0.5) ./ (c1[wti] .+ .5))
    fv = (c0 .+c1 .+1) ./ ((c0 .+ .5) .* (c1 .+ .5))
    fv[df[!, seqname] .!= wt] = fv[df[!, seqname] .!= wt] .+ fv[wti]
    Dict(:x => x, :c0 => c0, :c1 => c1, :code => code, :fv => fv, :wt => wt,
                :ham=> ham, :f =>f, :pos => pos, :aa => aa)
end

function energymodel(data, kind, maxe = Inf; kwargs...)
    l = size(data[:x], 2)
   
    energies = DataFrame(efi = zeros(l))
    pb0 = 0.01
    if kind == :fold_and_bind
        energies[!, :ebi] = zeros(l)
    end
    r = 1.0
    mi = Dict(:r => r, :pb0 => pb0, :kind => kind, :data => data,
        :energies => energies, :maxe => maxe)
    energymodel(mi; kwargs...)
end

function energymodel(mi; 
		alg = :LD_MMA, 
		maxit = 1000000, 
		tol = 1e-12,
		fixparameter = nothing,
		fixparametervalue = 0.0,
		fixparameterindex = 0,
		normalize = true
	)
    data = mi[:data]
    kind = mi[:kind]
    x = data[:x]
    c0 = data[:c0]
    c1 = data[:c1]
    n, l = size(x)

    efi = (mi[:energies][!, :efi])
    r = mi[:r]
    pb0 = mi[:pb0]

    efr = 2:(l+1)
    mem1 = zeros(n)
    mem2 = zeros(n)
    mem3 = zeros(n)
    mem4 = zeros(n)
    mll = zeros(n)
    best = [-Inf]
    prog = nothing
    if kind == :fold_and_bind
        ebr = (l+2):(2l+1)
        ebi = (mi[:energies][!, :ebi])
        pi = [log(r); efi; ebi; log(pb0)]
        mem6 = zeros(n)
        mem7 = zeros(n)
        opt = Opt(alg, length(pi))
        max_objective!(opt, (p, g) -> fold_and_bind_obj(p, g, 
                n, efr, ebr, x, c0, c1, mem1, mem2, mem3, mem4, mem7,
                mem6, mll, prog, best, normalize))
    elseif kind == :fold
        pi = [log(r); efi; log(pb0)]
        opt = Opt(alg, length(pi))
        max_objective!(opt, (p, g) -> fold_obj(p, g, n, efr, x, c0, c1,
                mem1, mem2, mem3, mem4, mll, prog, best, normalize))
    end    
    
    maxeval!(opt, maxit)
#     xtol_abs!(opt, default_tol)
    ftol_rel!(opt, tol)

    maxe = mi[:maxe]
    # lbounds = -maxe * ones(pi)
    lbounds = -maxe * ones(size(pi))
    lbounds[1] = -Inf
    lbounds[end] = -Inf
    # ubounds = maxe*ones(pi)
    ubounds = maxe * ones(size(pi))
    ubounds[end] = Inf
    ubounds[1] = Inf
    if fixparameter == :r
        ifixpar = 1
    elseif fixparameter == :pb0
        ifixpar = length(pi)
    elseif fixparameter == :efi
        ifixpar = efr[fixparameterindex]
    elseif fixparameter == :ebi
        ifixpar = ebr[fixparameterindex]
    end
    if fixparameter != nothing
        lbounds[ifixpar] = fixparametervalue
        ubounds[ifixpar] = fixparametervalue
        pi[ifixpar] = fixparametervalue
    end
    if any(lbounds != -Inf)
        lower_bounds!(opt, lbounds)
    end
    if any(ubounds != Inf)
        upper_bounds!(opt, ubounds)
    end
    (ll, p, ret) =  optimize!(opt, pi)
    if normalize
        ll = ll*n
    end
    if ret != :FTOL_REACHED
        @warn("FTOL not reached, $ret")
    end
    
    dout = Dict(:kind => kind, :ll => ll, :pb0 => exp(p[end]), :r => exp(p[1]), 
        :return => ret, :data => data, :maxe=>maxe)
    
    dout[:energies] = DataFrame(pos = data[:pos], aa = data[:aa], efi = p[efr])
    #energies[:efi] = p[efr]
    if kind == :fold_and_bind
        dout[:energies][!, :ebi] = p[ebr]
    end
    dout[:pred] = predict(dout, x) 
    
    return dout
end

function predict(m, x::SparseMatrixCSC{Float64,Int64})
    prediction = DataFrame(ef = x*m[:energies][!, :efi])
    p0 = m[:pb0]
    ex = exp.(prediction[!, :ef])
    pwt = exp(m[:energies][!, :efi][1])
    if m[:kind] == :fold_and_bind
        prediction[!, :eb] = x*m[:energies][!, :ebi]
        ex = exp.(prediction[!, :eb]) .* (1 .+ ex)
        pwt = exp(m[:energies][!, :ebi][1]) .* (1 .+ pwt)
    end
    p = (1 .+ p0 .* ex) ./ (1 .+ ex)
    pwt = (1 .+ p0 .* pwt) ./ (1 .+ pwt)
    prediction[!, :p] = p
    prediction[!, :f] = log.(p/pwt)
    return prediction
end

function boot(m, kind=:naive; kwargs...)
    mboot = copy(m)
    mboot[:data] = copy(m[:data])
    rp = m[:r]*m[:pred][!, :p]
    if kind == :thermo
        mboot[:data][:c0] = rand.(Poisson.((m[:data][:c0].+m[:data][:c1])./(1+rp)))
        mboot[:data][:c1] = rand.(Poisson.((m[:data][:c0].+m[:data][:c1])./(1+rp).*rp))
    elseif kind == :naive
        mboot[:data][:c0] = rand.(Poisson.(m[:data][:c0]))
        mboot[:data][:c1] = rand.(Poisson.(m[:data][:c1]))
    end
    energymodel(mboot; kwargs...)
end

function bootsearch(mi, i, tol=1e-4, cb = nothing; kwargs...)
    mdls = []
    mbest = copy(mi)
    mbest[:data] = mi[:data]
    while length(mdls) < i
        mboot = boot(mbest, :naive; kwargs...)
        mboot[:data] = mi[:data]
        m = energymodel(mboot; kwargs...)
        dll = m[:ll]-mbest[:ll]
        if dll > tol
            mbest = m
            mdls = []
        else
            delete!(mboot, :pred)
            push!(mdls, mboot)
        end
        # display("delta $dll, i $(length(mdls))")
        print("delta $dll, i $(length(mdls))")
        if cb != nothing
            cb(mbest)
        end
    end
    mbest[:boot] = mdls
    mbest[:energies_boot] = boot_stat(mbest)
    # display("total delta LL $(mbest[:ll]-mi[:ll])")
    print("total delta LL $(mbest[:ll]-mi[:ll])")
    return mbest
end


function boot_stat(m)
    pboot = deepcopy(m[:boot][1][:energies])
    for mb in m[:boot][2:end]
        append!(pboot, mb[:energies])
    end
    ci = (d,s) -> DataFrame( 
                            med = median(d[s]), 
                            upper = quantile(d[s], .975), 
                            lower = quantile(d[s], .025))
    #return pboot
    bbs = by(pboot, [:pos, :aa], d -> ci(d, :efi))
    rename!(bbs, :med => :f_med)
    rename!(bbs, :upper => :f_upper)
    rename!(bbs, :lower => :f_lower)
    if m[:kind] == :fold_and_bind
        bbsb = by(pboot, [:pos, :aa], d -> ci(d, :ebi))
        rename!(bbsb, :med => :b_med)
        rename!(bbsb, :upper => :b_upper)
        rename!(bbsb, :lower => :b_lower)
        bbs = join(bbs, bbsb, on = [:pos, :aa])
    end
    bbs = join(m[:energies], bbs, on = [:pos,:aa])
    #bbs = vcat(bbs, stack(m[:energies], [:efi, :ebi]))
    return bbs
end

function fold_and_bind_obj(p, g, n, efr, ebr, x, c0, c1, ef, eb, 
        dldef, dldeb, dr, dp0, ll, prog, best, normalize)
    efi = view(p, efr)
    ebi = view(p, ebr)
    r = exp(p[1])
    pb0 = exp(p[end])
    # A_mul_B!(ef, x, efi)
    mul!(ef, x, efi)
    # A_mul_B!(eb, x, ebi)
    mul!(eb, x, ebi)
    Threads.@threads for i = 1:n
        exfb = exp(ef[i]+eb[i])
        exfbb = exfb + exp(eb[i])
        r1px = r*(1 + pb0*exfbb)
        nt = c0[i] + c1[i]
        ll[i] = -nt*log(1+exfbb+r1px) + c0[i]*log1p(exfbb) + c1[i]*log(r1px)

        rdldr = -nt*r1px/(1+exfbb+r1px) + c1[i]
        dr[i] = rdldr
        dp0[i] = rdldr*exfbb/(1+pb0*exfbb)*pb0
        temp = rdldr*(pb0-1)/(1+pb0*exfbb)/(1+exfbb)
        dldef[i] = exfb*temp
        dldeb[i] = exfbb*temp
    end
    g[1] = sum_kbn(dr)
    g[end] = sum_kbn(dp0)
    # gef = view(g, efr)
    # At_mul_B!(gef, x, dldef)
    # gef = transpose(x) * dldef
    g[efr] = transpose(x) * dldef
    g[efr[1]] = sum_kbn(dldef)
    # geb = view(g, ebr)
    # At_mul_B!(geb, x, dldeb)
    # geb = transpose(x) * dldeb
    g[ebr] = transpose(x) * dldeb
    g[ebr[1]] = sum_kbn(dldeb)
    if normalize
        g .= g./n
        dll = mean(ll)
    else
        dll = sum_kbn(ll)
    end
    return dll
end

#load data
datafr = CSV.read(parsed_args["dataset"], delim=" ")

#prepare data
data = prepdata(datafr, :Mut, :listofmuts, "0X", :DNA, :SelAll)

#calculate model
m = energymodel(data, :fold_and_bind, 15)

m2 = bootsearch(m, 100, 1e-4)

#save data
CSV.write(string(parsed_args["name"],"_3state_energies.csv"),
	m2[:energies]
)
CSV.write(string(parsed_args["name"],"_3state_globalpar.csv"),
	DataFrame(
		b_scale = m2[:r],
		b_bgr = m2[:pb0],
		s_dGwt = m2[:energies][1,:efi],
		b_dGwt = m2[:energies][1,:ebi]
	)
)
CSV.write(string(parsed_args["name"],"_3state_energies_boot.csv"),
    m2[:energies_boot])

