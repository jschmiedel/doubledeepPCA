## Otwinowski three state model from binding measurements

export JULIA_NUM_THREADS=14
julia

# first, implement code from jupyter notebook locally

cd("Dropbox (Personal)/Science/CRG/DMS2struct/datasets/Otwinoski2018_ProteinGthermo-master")
cd("doubledeepPCA/Otwinoski2018_ProteinGthermo-master")

######### cell 1
# using DataFrames
# using JLD
# @load "precomputedmodels.jld"

######### #cell2
# using CSV
# using DataFrames
# smut = CSV.read("SMutList.txt", delim="\t")
# data0 = smut[1,:]
# data0[:Mut] = "WT"
# data0[:ham] = 0
# ns = names(smut)
# for (n1,n2) in zip(ns[2:7], ns[8:13])
#     data0[n1] = data0[n2]
# end
# smut[:ham] = 1
# smut = vcat(data0, smut)
# delete!(smut, [:Pos; ns[8:13]])
# dmut = CSV.read("DMutList.txt", delim="\t")
# delete!(dmut, ns[8:13])
# dmut[:ham] = 2
# datafr = vcat(smut, dmut)
# note that positions in :Mut are off by 1; you should +1 these positions

#use assembled data frame (from Otwinowski2018_prep_data.R)
using CSV
using DataFrames
using SparseArrays
using DelimitedFiles
using NLopt
import StatsBase.predict
using Distributions
using KahanSummation
using LinearAlgebra

##load data
datafr = CSV.read("GB1_all_data.txt", delim=" ")
#downsampled dataframe
datafr = CSV.read("GB1_singles_all_doubles_5k.txt", delim=" ")


######### #cell3
using RCall

R"""
library(ggplot2)
library(cowplot)
library(dplyr)
theme_set(theme_cowplot(11))
options(device = function(filename=getOption('rcalljl_filename'),...) png(filename, width=800, height=800, ...))
""";

######### #cell4
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

######### #cell5
# data = prepdata(datafr, :Mut, :listofmuts, "WT", :DNA, :SelAll)
data = prepdata(datafr, :Mut, :listofmuts, "0X", :DNA, :SelAll)
# data = prepdata(df = datafr, seqname = :Mut, kind = :listofmuts, wt = "0X", c0n = :DNA, c1n = :SelAll, delim = '-')



######### #cell 6

# energymodel(data; kind = :fold; maxe = 15; maxit = 10)
# energymodel(data, :fold, 15, maxit = 10000)

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


function fold_obj(p, g, n, efr, x, c0, c1, ef, dLde,
        dr, dp0, ll, prog, best, normalize)
    r = exp(p[1])
    pb0 = exp(p[end])
    # A_mul_B!(ef, x, view(p, efr))
    mul!(ef, x, view(p, efr))

    Threads.@threads for i = 1:n
        exf = exp(ef[i])
        ep = pb0*exf
        r1px = r*(1 + ep)
        nt = c0[i] + c1[i]
        ll[i] = -nt*log(1+exf+r1px) + c0[i]*log1p(exf) + c1[i]*log(r1px)
        temp3 = -nt*r1px/(1+exf+r1px) + c1[i]
        dr[i] = temp3
        temp = temp3*exf/(1+ep)
        dp0[i] = temp*pb0
        dLde[i] = -temp*(1-pb0)/(1+exf)
    end

    ## variable g not defined?!
    g[1] = sum_kbn(dr)
    g[end] = sum_kbn(dp0)
    # At_mul_B!(view(g, efr), x, dLde)
    g[efr] = transpose(x) * dLde
    g[2] = sum_kbn(dLde)

    dll = sum_kbn(ll)
    if normalize
        g .= g./n
        dll = dll/n
    end
    
    return dll
end


function predict(m, x::SparseMatrixCSC{Float64,Int64})
    prediction = DataFrame(ef = x*m[:energies][:efi])
    p0 = m[:pb0]
    ex = exp.(prediction[:ef])
    pwt = exp(m[:energies][!,:efi][1])
    if m[:kind] == :fold_and_bind
        prediction[:eb] = x*m[:energies][:ebi]
        ex = exp.(prediction[:eb]) .* (1 .+ ex)
        pwt = exp(m[:energies][:ebi][1]) .* (1 .+ pwt)
    end
    p = (1 .+ p0 .* ex) ./ (1 .+ ex)
    pwt = (1 .+ p0 .* pwt) ./ (1 .+ pwt)
    prediction[:p] = p
    prediction[:f] = log.(p/pwt)
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
        display("delta $dll, i $(length(mdls))")
        if cb != nothing
            cb(mbest)
        end
    end
    mbest[:boot] = mdls
    mbest[:energies_boot] = boot_stat(mbest)
    display("total delta LL $(mbest[:ll]-mi[:ll])")
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

# warning can take a long time
@time m1 = energymodel(data, :fold, 15, maxit = 100000)
# this can take even longer
m1 = bootsearch(m1, 100, 1e-4)


######### #cell 7
R"""
plot.yyhat1 = qplot(x=$(m1[:pred][:f]), y=$(m1[:data][:f]), geom="bin2d", bins=50) +
geom_abline(alpha=0.5) + coord_fixed() +
scale_fill_distiller(palette=7, direction=1) +
xlab("inferred fitness") + ylab("fitness")
ggsave(filename = 'twostate_logbindingprediction.pdf',width=6)
"""
######### #cell 8
R"""
qplot(data=$(m1[:pred]), x=ef, y=$(m1[:data][:f]), geom="bin2d", bins=200) +
geom_line(aes(y=f), alpha=0.3, linetype=2) +
scale_fill_distiller(palette=7, direction=1) +
xlab("prediction") + ylab("measurement") +
stat_summary_bin(geom="line", fun.y=median, alpha=0.3)
ggsave(filename = 'twostate_foldingenergyVSlogbinding.pdf')
"""

######### #cell 9
# stab = CSV.read("stabilities.csv")

# stab[!, :pos] = [parse(Int64, match(r"\d+", bn).match) for bn in stab[!, :mutation]]
# stab[!, :aa] = [match(r"(?<=\d)[A-Z]|\*", bn).match for bn in stab[!, :mutation]] #doesn't work because this creates a SubString
# estab1 = join(stab, m1[:energies], on = [:pos, :aa], kind = :inner)
# for s in [:efi]#, :f_lower, :f_upper]
#     estab1[!, s] = -estab1[!, s] * 0.001987 * (273+24)
# end
# R"print(cor.test($estab1$efi, $estab1$ddg))"

# R"""
# plot.ddg1 = qplot(data=$estab1, x=efi, y=ddg, alpha=I(.5)) +
# xlab("predicted energy kcal/mol") +
# ylab("measured energy kcal/mol") +
# coord_fixed() + geom_abline() #+ geom_errorbarh(aes(xmin=f_lower, xmax=f_upper))
# """

######### #cell 10
# ef = mem1
# eb = mem2
# dldef = mem3
# dldeb = mem4
# dr = mem7
# dp0 = mem6
# ll = mll
# p = pi
# g = p
# fold_and_bind_obj(p, g, 
#                 n, efr, ebr, x, c0, c1, mem1, mem2, mem3, mem4, mem7,
#                 mem6, mll, prog, best, normalize))
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
@time m3 = energymodel(data, :fold_and_bind, 15)

CSV.write("Ot_GB1_full_3state_energies.csv",m3[:energies])
CSV.write("Ot_GB1_full_3state_globalpar.csv",DataFrame(r = m3[:r],pb0 = m3[:pb0],s_dGwt = m3[:energies][1,:efi],b_dGwt = m3[:energies][1,:ebi]))
	
@time m4 = bootsearch(m3, 100, 1e-2)
CSV.write("Ot_GB1_ds1_boot_3state_energies.csv",m4[:energies])
CSV.write("Ot_GB1_ds1_boot_3state_energies_boot.csv",m4[:energies_boot])
CSV.write("Ot_GB1_ds1_boot_3state_globalpar.csv",DataFrame(r = m4[:r],pb0 = m4[:pb0],s_dGwt = m4[:energies][1,:efi],b_dGwt = m4[:energies][1,:ebi]))
