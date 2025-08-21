const mut_choice_map = reshape(1:12 |> collect,(3,4))

const powerset_topologies = [[bit == '1' ? -1 : bit == '2' ? 0 : 1 for bit in string(i;base = 3,pad = 10)] for i in 0:3^10-1];

mutable struct Trajectory

    sim_id :: Int64

    geno_traj :: Vector{Vector{Float64}}
    topologies :: Vector{Vector{Int64}}

    n_accepted_mutants :: Int64
    acceptance_ratio :: Float64

    mutation_number :: Vector{Int64}
    stripe_indicator :: Vector{Bool}
    H0 :: Int64
    wait_times :: Vector{Int64}

    fitness_traj :: Vector{Float64}
    fitness_traj_tuple :: Any

    top_edits :: Vector{Int64}
    weight_edits :: Any
    masked_hamming_distance_H0 :: Any
    masked_hamming_distance :: Any

    initial_fitness :: Float64
    final_fitness :: Float64
    mutant_info :: Any

    inc_metagraph_vertices :: Any
    inc_metagraph_parents :: Any 
    minimal_stripe_subgraphs :: Any

    parent_inclusion_indicator :: Any

    tt_label_probabilities :: Any
    tt_label_predictions :: Any
    tt_label_entropies :: Any
    tt_prediction_error :: Any

    gt_label_probabilities :: Any
    gt_label_predictions :: Any
    gt_label_entropies :: Any
    gt_prediction_error :: Any

    mss_probabilities :: Any
    mss_predictions :: Any
    mss_entropies :: Any
    mss_prediction_error :: Any

    train_test_indicator :: Any
    train_test_indicator_mss :: Any
    embeddings :: Any

    epistasis :: Any

    tt_kl_div :: Any
    gt_kl_div :: Any
    mss_kl_div :: Any

    gt_shap :: Any
    tt_shap :: Any
    mss_shap :: Any

    other :: Any
    
end

function Trajectory(sim_id::Int64,geno_traj_m::Matrix{Float64},fitness_traj_tuple::Vector{Tuple{Float64, Float64}},wait_times,mut_choices::Any,mut_types::Any,mut_sizes::Any,weight_names)

    geno_traj = [collect(i) for i in eachcol(geno_traj_m)]

    ######## process trajectory data

    topologies = map(w->Int.(sign.(w)),geno_traj)

    wait_times_v = vcat([1.],wait_times)

    n_accepted_mutants = length(fitness_traj_tuple)-1
    n_generated_mutants = length(fitness_traj_tuple)-1
    acceptance_ratio = n_accepted_mutants / n_generated_mutants

    mutation_number = [i for i in 0:n_accepted_mutants]
    stripe_indicator =  map(ft->ft[1] == 0,fitness_traj_tuple)
    H0 = minimum(findall(stripe_indicator))

    fitness_traj_add = map(ft->add_fitness(ft),fitness_traj_tuple)
    top_edits = compute_cumulative_edits(reduce(hcat,topologies))
    weight_edits = nothing 
    masked_hamming_distance_H0 = nothing
    masked_hamming_distance = nothing

    initial_fitness = fitness_traj_add[1]
    final_fitness = fitness_traj_add[end]

    ######## MST

    minimal_stripe_subgraphs = nothing 
    metagraph_vertices = nothing
    metagraph_parents = nothing

    ######## generate mutation data

    weight_id = map(mc-> [mut_choice_map[i] for i in mc],mut_choices)

    fitness_delta = get_fitness_delta(fitness_traj_add)

    n_weight_changes = map(m->length(m),weight_id)

    weight_id_label = map(v->join(map(x->weight_names[x],v),"|"),weight_id)

    start_fitness = fitness_traj_add[1:end-1]
    start_fitness_tuple = fitness_traj_tuple[1:end-1]
    mutant_fitness = fitness_traj_add[2:end]

    start_network = [collect(v) for v in eachcol(geno_traj_m[:,1:end-1])]

    is_new_int = [[sn[wi] == 0 for wi in mw_id] for (sn,mw_id) in zip(start_network,weight_id)]

    mut_types_combi = [[t for t in zip(new_int,mt)] for (new_int,mt) in zip(is_new_int,mut_types)]

    mutant_info = [(weight_id = weight_id[n],weight_id_label = weight_id_label[n],mut_type = mut_types_combi[n], new_interaction = is_new_int[n], mut_size = mut_sizes[n],start_fitness = start_fitness[n],start_fitness_tuple = start_fitness_tuple[n],mutant_fitness = mutant_fitness[n],fitness_delta = fitness_delta[n],start_network = start_network[n]) for n in 1:length(weight_id)]

    ####### predictions

    parent_inclusion_indicator = nothing

    tt_label_probabilities = nothing
    tt_label_predictions = nothing
    tt_label_entropies = nothing
    tt_prediction_error = nothing

    gt_label_probabilities = nothing
    gt_label_predictions = nothing
    gt_label_entropies = nothing
    gt_prediction_error = nothing

    mss_probabilities = nothing
    mss_predictions = nothing
    mss_entropies = nothing
    mss_prediction_error = nothing

    train_test_indicator = nothing
    train_test_indicator_mss = nothing
    embeddings = nothing

    epistasis = nothing

    tt_kl_div = nothing
    gt_kl_div = nothing
    mss_kl_div = nothing

    tt_shap = nothing
    gt_shap = nothing
    mss_shap = nothing

    other = nothing

    ####### instantiate 

    Trajectory(sim_id,geno_traj,topologies,n_accepted_mutants,acceptance_ratio,mutation_number,stripe_indicator,H0,wait_times_v,fitness_traj_add,fitness_traj_tuple,top_edits,weight_edits,masked_hamming_distance_H0,masked_hamming_distance,initial_fitness,final_fitness,mutant_info,metagraph_vertices,metagraph_parents,minimal_stripe_subgraphs,parent_inclusion_indicator,
                                                                                                            tt_label_probabilities,tt_label_predictions,tt_label_entropies,tt_prediction_error,gt_label_probabilities,gt_label_predictions,gt_label_entropies,gt_prediction_error,
                                                                                                            mss_probabilities,mss_predictions,mss_entropies,mss_prediction_error,train_test_indicator,train_test_indicator_mss,embeddings,epistasis,tt_kl_div,gt_kl_div,mss_kl_div,tt_shap,gt_shap,mss_shap,other)
end

function add_fitness(tuple_f)
    return tuple_f[1] + ((tuple_f[2]+1)/2) + 1
end

function create_full_fitness_traj(fitness_traj,wait_times)

    all_ff = []
    
    for (fitness,wt) in zip(fitness_traj[1:end-1],wait_times[2:end])

        full_fitness  = fill(fitness,wt)

        push!(all_ff,full_fitness)

    end

    all_ff_v = reduce(vcat,all_ff)

    all_ff_vf = vcat(all_ff_v,[fitness_traj[end]])

    return all_ff_vf
end

function return_fitness_delta(fitness_orig,fitness_mutant)
    Δf1 = fitness_orig[1] - fitness_mutant[1]
    Δf2 = fitness_orig[2] - fitness_mutant[2]

    return Δf1 != 0. ? Δf1 : Δf2
end

function compute_cumulative_edits(gt)

    dham = pairwise(Hamming(),gt)

    total_ham = 0

    cumulative_ham = []

    push!(cumulative_ham,total_ham)

    for i in 1:size(gt,2)-1

        total_ham += dham[i,i+1]

        push!(cumulative_ham,total_ham)
    end

    return Int64.(cumulative_ham)

end

function masked_hamming_distance(topology,target_topology)

    new_topology = copy(topology)

    z0 = findall(x->x == 0,target_topology)

    new_topology[z0] .= 0.

    Distances.evaluate(Hamming(),new_topology,target_topology)

end

function assign_minimal_subgraphs!(tr::Trajectory,fs,ls)
    tr.minimal_stripe_subgraphs = fill([],length(tr.topologies))
    tr.minimal_stripe_subgraphs[tr.H0] = fs
    tr.minimal_stripe_subgraphs[end] = ls

    tr.masked_hamming_distance_H0 = [[masked_hamming_distance(top,Int.(sign.(fs))) for top in tr.topologies],[masked_hamming_distance(tr.topologies[1],top) for top in tr.topologies]]
    tr.masked_hamming_distance = [[masked_hamming_distance(top,Int.(sign.(ls))) for top in tr.topologies],[masked_hamming_distance(tr.topologies[1],top) for top in tr.topologies]]
end

function create_inclusion_metagraph(trajectories::Vector{Trajectory})

    min_stripe_top = unique(reduce(vcat,[tr.minimal_stripe_subgraphs[[tr.H0,end]] for tr in trajectories]))

    n_min_stripe_top = length(min_stripe_top)

    min_stripe_top_complexity = map(top->sum(abs.(top)),min_stripe_top)

    min_stripe_top_ordered = min_stripe_top[sortperm(min_stripe_top_complexity)]

    vertex_top_map = Dict(n=>top for (n,top) in enumerate(min_stripe_top_ordered));
    vertex_complexity_map = Dict(n=>sum(abs.(top)) for (n,top) in enumerate(min_stripe_top_ordered));
    top_vertex_map = Dict(top=>n for (n,top) in enumerate(min_stripe_top_ordered));

    inclusion_matrix = zeros(Int64,(n_min_stripe_top,n_min_stripe_top))

    for i in 1:n_min_stripe_top
        for j in 1:n_min_stripe_top
            if i != j
                inclusion_matrix[i,j] = test_inclusion(min_stripe_top_ordered[j],min_stripe_top_ordered[i])
            end
        end
    end

    inc_metagraph = SimpleDiGraph(inclusion_matrix)

    return inc_metagraph, vertex_top_map, top_vertex_map, vertex_complexity_map,inclusion_matrix
end

function assign_inc_vertex_ids!(tr::Trajectory,top_vertex_map)
    tr.inc_metagraph_vertices = [n ∈ [tr.H0,length(tr.topologies)] ? top_vertex_map[tr.minimal_stripe_subgraphs[n]] : -1 for n in 1:length(tr.topologies)]
end

function assign_inc_parents!(tr::Trajectory,inclusion_matrix,vertex_complexity_map,minimal_motif_id)

    tr.inc_metagraph_parents = []

    for n in 1:length(tr.topologies)

        if n ∈ [tr.H0,length(tr.topologies)]

            if tr.inc_metagraph_vertices[n] ∈ minimal_motif_id

                push!(tr.inc_metagraph_parents,tr.inc_metagraph_vertices[n])
            else
                options = minimal_motif_id[findall(inclusion_matrix[minimal_motif_id,tr.inc_metagraph_vertices[n]] .== 1)]
                choice = argmax([vertex_complexity_map[v] for v in options])

                push!(tr.inc_metagraph_parents,options[choice])
            end
        else
            push!(tr.inc_metagraph_parents,-1)
        end 
    end
end

function assign_weight_edits!(tr)

    we = vcat([0],map(mi->length(mi[:weight_id]),tr.mutant_info))

    c_we = [sum(we[1:i+1]) for i in 0:length(we)-1]

    tr.weight_edits = c_we
end

function select_minimal_topologies(list_mss)
    if length(list_mss) == 1
        return sign.(list_mss[1][3])
    else
        id = argmax(map(x->x[2],list_mss))
        return sign.(list_mss[id][3])
    end
end

function assign_predictions!(tr::Trajectory,model,prediction_type)

    if prediction_type == :tt
        tt = reduce(hcat,tr.topologies)
        tt_dtrain = xgboost.DMatrix(tt[1:10,:] |> transpose |> collect, feature_types = c_types, feature_names = weight_names)
        tr.tt_label_probabilities = model.predict(tt_dtrain)
        tr.tt_label_predictions = mapslices(p->argmax(p),tr.tt_label_probabilities,dims = 2)
        tr.tt_label_entropies = mapslices(p->entropy(p),tr.tt_label_probabilities,dims = 2);

    elseif prediction_type == :gt
        gt = reduce(hcat,tr.geno_traj)
        gt_dtrain = xgboost.DMatrix(gt[1:10,:] |> transpose |> collect, feature_names = weight_names)
        tr.gt_label_probabilities = model.predict(gt_dtrain)
        tr.gt_label_predictions = mapslices(p->argmax(p),tr.gt_label_probabilities,dims = 2)
        tr.gt_label_entropies = mapslices(p->entropy(p),tr.gt_label_probabilities,dims = 2);
    else
        gt = reduce(hcat,tr.geno_traj)
        gt_dtrain = xgboost.DMatrix(gt[1:10,:] |> transpose |> collect, feature_names = weight_names)

        prediction_prob = []
        prediction_labels = []

        for m in model
            edge_prob = m.predict(gt_dtrain)
            edge_label = mapslices(x->argmax(x)-2 ,edge_prob,dims = 2)

            push!(prediction_prob,edge_prob)
            push!(prediction_labels,edge_label)
        end

        tr.mss_probabilities = reduce((x,y) -> cat(x,y,dims = 3),[reshape(mss_p,(size(mss_p)...,1)) for mss_p in prediction_prob])
        tr.mss_predictions = [r |> collect for r in eachrow(reduce(hcat,prediction_labels))]
        
        mss_entropies = []

        for n in 1:length(tr.topologies)
            ps = @view tr.mss_probabilities[n,:,:]
            e = entropy([calculate_probability(ps,t) for t in powerset_topologies])
            push!(mss_entropies,e)
        end

        tr.mss_entropies = mss_entropies
    end

end

function assign_predictions!(tr::Trajectory,model,prediction_type,predict_label_to_vertex)

    if prediction_type == :tt
        tt = reduce(hcat,tr.topologies)
        tt_dtrain = xgboost.DMatrix(tt[1:10,:] |> transpose |> collect, feature_names = weight_names)
        tr.tt_label_probabilities = model.predict(tt_dtrain)
        tr.tt_label_predictions = mapslices(p->predict_label_to_vertex[argmax(p)],tr.tt_label_probabilities,dims = 2)
        tr.tt_label_entropies = mapslices(p->entropy(p),tr.tt_label_probabilities,dims = 2);

    elseif prediction_type == :gt
        gt = reduce(hcat,tr.geno_traj)
        gt_dtrain = xgboost.DMatrix(gt[1:10,:] |> transpose |> collect, feature_names = weight_names)
        tr.gt_label_probabilities = model.predict(gt_dtrain)
        tr.gt_label_predictions = mapslices(p->predict_label_to_vertex[argmax(p)],tr.gt_label_probabilities,dims = 2)
        tr.gt_label_entropies = mapslices(p->entropy(p),tr.gt_label_probabilities,dims = 2);
    else
        gt = reduce(hcat,tr.geno_traj)
        gt_dtrain = xgboost.DMatrix(gt[1:10,:] |> transpose |> collect, feature_names = weight_names)

        prediction_prob = []
        prediction_labels = []

        for m in model
            edge_prob = m.predict(gt_dtrain)
            edge_label = mapslices(x->argmax(x)-2 ,edge_prob,dims = 2)

            push!(prediction_prob,edge_prob)
            push!(prediction_labels,edge_label)
        end

        tr.mss_probabilities = reduce((x,y) -> cat(x,y,dims = 3),[reshape(mss_p,(size(mss_p)...,1)) for mss_p in prediction_prob])
        tr.mss_predictions = [r |> collect for r in eachrow(reduce(hcat,prediction_labels))]
        
        mss_entropies = []

        for n in 1:length(tr.topologies)
            ps = @view tr.mss_probabilities[n,:,:]
            e = entropy([calculate_probability(ps,t) for t in powerset_topologies])
            push!(mss_entropies,e)
        end

        tr.mss_entropies = mss_entropies
    end

end

function assign_tt_other_prediction_errors!(tr::Trajectory,label,predict_id)

    tr.tt_prediction_error = [label == -1 ? !(pred ∈ predict_id) :  pred == label for pred in tr.tt_label_predictions]

end

function assign_gt_other_prediction_errors!(tr::Trajectory,label,predict_id)

    tr.gt_prediction_error = [label == -1 ? !(pred ∈ predict_id) :  pred == label for pred in tr.gt_label_predictions]

end

function return_order_by_count(v)

    v_un = unique(v)
    counts_v = [count(x->x==value,v) for value in v_un]

    order_v = sortperm(counts_v,rev = true)

    return v_un[order_v],counts_v[order_v]
end

function get_fitness_delta(f_traj)
    f_traj_diff = f_traj[2:end] .- f_traj[1:end-1]
    return f_traj_diff
end

##### Assignment tools

function test_inclusion(net_v,top_v)

    n1 = sign.(net_v)
    incl = 1

    for (n,w) in enumerate(top_v)

        if w != 0

            if n1[n] != w
                incl = 0
            end
        end
    end

    return incl
end

function assign_class(row)
    if all(row.==0)
        return "no assignment"
    else
        return fundamental_topologies[findall(x->x==1,row)]
    end
end

function determine_class(en_top,dyn_top)

    r = zeros(Int,size(en_top,1))

    for i in 1:size(en_top,1)
        id =  findall(x->x==1,en_top[i,:])
        if length(id) == 0
            r[i] = 0
        elseif length(id) == 1
            r[i] = id[1]
        else
            comb = en_top[i,:] .& dyn_top[i,:]

            if all(comb.==0)
                r[i] = 0
            else
                r[i] = findall(x->x==1,comb)[1]
            end
        end
    end

    return r
end

function entropy(v)
    length(v) != 1 ? -sum(v[v .!= 0] .* log.(v[v .!= 0])) / log(length(v)) : -sum(v[v .!= 0] .* log.(v[v .!= 0]))
end

function cross_entropy(p,q)

    total_entropy = 0

    for (pi,qi) in zip(p,q)
        if (pi == 0) || (qi == 0)
            nothing
        else
            total_entropy += pi*log(qi)
        end
    end
    
    return -total_entropy
end

#### epistasis

function type_epi(df_ab,df_aB,df_Ab,df_AB,fitness_eps)

    if abs(df_AB) < fitness_eps
        return :Neutral
    
    else

        # aA1 = Ab - ab # df_Ab
        # aA2 = AB - aB # (f_ab + df_AB) - (f_ab + df_aB) = df_AB - df_aB

        # bB1 = aB - ab # df_aB
        # bB2 = AB - Ab # (f_ab + df_AB) - (f_ab + df_Ab) = df_AB - df_Ab

        aA1 = df_Ab # df_Ab
        aA2 = df_AB - df_aB # (f_ab + df_AB) - (f_ab + df_aB) = df_AB - df_aB

        bB1 = df_aB # df_aB
        bB2 = df_AB - df_Ab # (f_ab + df_AB) - (f_ab + df_Ab) = df_AB - df_Ab

        aA1 = abs(aA1) < fitness_eps ? fitness_eps : aA1
        aA2 = abs(aA2) < fitness_eps ? fitness_eps : aA2
        bB1 = abs(bB1) < fitness_eps ? fitness_eps : bB1
        bB2 = abs(bB2) < fitness_eps ? fitness_eps : bB2

        if (bB1 * bB2 > 0) & (aA1 * aA2 > 0)

            if abs((df_Ab + df_aB) - df_AB) < fitness_eps
                return :ne
            else
                return :me
            end

        elseif (bB1 * bB2 < 0) & (aA1 * aA2 < 0)
            return :rse

        else
            return :se
        end
    end
end

function characterise_weight_int(weights,mst_weight)

    a = [i in mst_weight for i in weights]

    if a==[true,false]
        return [false,true]
    else
        a
    end
end

function characterise_weight_corr(weights,weight_correlations)
    abs(weight_correlations[Tuple(weights)...])
end

function calculate_mut_type_proportion_list(mut_type,types)
    
    total = length(mut_type)

    return [count(x->x ∈ ec,mut_type)/total for ec in types]
end

function get_mut_n(tr,range_l,range_u)

    d = [length(mi[:weight_id]) for mi in tr.mutant_info[range_l:range_u]]

    if length(d) >= 1
        return reduce(vcat,d)
    else
        return []
    end

end

function evaluate_epistasis_class_full_deltas(mut_tuple,grn_parameters,development,fitness_function,mut_op::MutationOperatorDual) # main

    n_mut = length(mut_tuple[:weight_id])

    if n_mut > 1

        mut_combi = [[bit == '1' ? true : false for bit in string(i;base = 2,pad = n_mut)] for i in 1:2^n_mut-1]

        prob_mutant = []

        for mut_id in mut_combi[1:end]

            new_network = noise_specified(mut_tuple[:start_network],mut_tuple[:weight_id][mut_id],mut_tuple[:mut_size][mut_id],mut_tuple[:mut_type][mut_id],mut_op)

            mutant = Individual(reshape(new_network,(3,4)),grn_parameters,development)

            mutant_fitness = fitness_function(mutant.phenotype)

            Δf1 = mutant_fitness[1] - mut_tuple[:start_fitness_tuple][1]
            Δf2 = mutant_fitness[2] - mut_tuple[:start_fitness_tuple][2]

            Δf = Δf1 != 0. ? Δf1 : Δf2

            push!(prob_mutant,Δf)

        end

        return vcat([0.],prob_mutant)

    else
        return [0.,1.]
    end
end

function characterise_epistasis_deltas(combi_result,fitness_eps) # main

    if length(combi_result) > 2
        ratio_new_mutant = combi_result ./ (combi_result[end] - fitness_eps)

        accept_new_mutant = ratio_new_mutant[2:end-1] .>= 1

        if !any(accept_new_mutant)
            rtype = :rse
        elseif all(accept_new_mutant)
            rtype = :ne
        else
            rtype = :se
        end

        return rtype,combi_result
    else
        return :sm, combi_result
    end
end

function evaluate_epistasis_types_full_deltas!(tr,grn_parameters,development,fitness_function,mut_op::MutationOperatorDual)
    all_class_epi = map(mi->evaluate_epistasis_class_full_deltas(mi,grn_parameters,development,fitness_function,mut_op),tr.mutant_info);
    tr.epistasis = all_class_epi
end

function return_fitness_delta(fitness_orig,fitness_mutant)
    Δf1 = fitness_orig[1] - fitness_mutant[1]
    Δf2 = fitness_orig[2] - fitness_mutant[2]

    return Δf1 != 0. ? Δf1 : Δf2
end

function extract_minimal_weights_with_fitness(epi_vector,fitness_eps)
    if length(epi_vector) > 1
        n_mut = Int(log2(length(epi_vector)))
        mut_combi = [[bit == '1' ? true : false for bit in string(i;base = 2,pad = n_mut)] for i in 1:2^n_mut-1]
        minimum(map(x->sum(x),mut_combi[epi_vector[2:end] .>= epi_vector[end] - fitness_eps])),n_mut
    else
        return NaN
    end
end


######## PredictionAnalytics

function find_restricted_mutant_id(we,restr)
    id = 0
    for w in we
        if w <= restr
            id +=1
        end
    end
    id
end 

function restricted_pred_scores_bins(trajectories,cohort_n,h_we,labels,null=false,model =:gt)

    traj_cohorts = map(tr->map(f->StatsBase.binindex(h_we, f),tr.weight_edits[1:tr.H0-1] ./ tr.weight_edits[tr.H0]),trajectories)

    # we bin the networks within each trajectory into bins representing the % of weight edits relative to S0. For each bin we collect the cohort of networks across trajectopries. 
    # If a trajectory has more than one network in that bin range we take the latest. We then calculate the prersited accuracy, f1 etc for each bin cohort

    in_cohort = findall(c->cohort_n ∈ c, traj_cohorts)

    labels_c = labels[in_cohort]
    
    if null
        restr_pred = [-1 for tr in trajectories[in_cohort]]
        restr_persist =  [true for tr in trajectories[in_cohort]]
        av_we = NaN
        restr_mut = NaN
    else
        restr_mut_id = map(tc->find_restricted_mutant_id(tc,cohort_n),traj_cohorts[in_cohort])

        av_we = mean([tr.weight_edits[id] ./ tr.weight_edits[tr.H0] for (tr,id) in zip(trajectories[in_cohort],restr_mut_id)])

        if model == :gt
            restr_pred = [tr.gt_label_predictions[id] for (tr,id) in zip(trajectories[in_cohort],restr_mut_id)]
            restr_persist =  [all(tr.gt_label_predictions[id:tr.H0-1] .== p) for (tr,(p,id)) in zip(trajectories[in_cohort],zip(restr_pred,restr_mut_id))]
        else
            restr_pred = [tr.tt_label_predictions[id] for (tr,id) in zip(trajectories[in_cohort],restr_mut_id)]
            restr_persist =  [all(tr.tt_label_predictions[id:tr.H0-1] .== p) for (tr,(p,id)) in zip(trajectories[in_cohort],zip(restr_pred,restr_mut_id))]
        end

        restr_mut = [tr.mutant_info[id:tr.H0-1] for (tr,id) in zip(trajectories[in_cohort],restr_mut_id)]
        
    end

    # per class recall

    class_recall = [] # = class accuracy

    all_tp_1 = []
    all_fp = []

    all_pred_class = []
    all_persist_class = []

    for n in 1:5

        pred_for_class_n = restr_pred[labels_c .== n] 
        correct_for_class_n = pred_for_class_n .== predict_label_to_vertex[n]
        persist_for_class_n = restr_persist[labels_c .== n] 

        TP = sum(correct_for_class_n .& persist_for_class_n)
        FP = length(pred_for_class_n) - TP

        recall = TP/(TP+FP)

        push!(class_recall,recall)
        push!(all_tp_1,TP)
        push!(all_fp,FP)
        push!(all_pred_class,pred_for_class_n)
        push!(all_persist_class,persist_for_class_n)
    end

    class_recall = [isnan(p) ? 0. : p for p in class_recall]

    # per class precision

    class_precision = []

    all_tp_2 = []
    all_fn = []

    all_class_pred = []
    all_persist_pred = []
    all_mutant_pred = []

    for n in 1:5

        class_for_pred_n = labels_c[(restr_pred .== predict_label_to_vertex[n]) .& restr_persist]
        correct_for_pred_n = class_for_pred_n .== n

        TP = sum(correct_for_pred_n)
        FN = length(class_for_pred_n) - TP

        precision = TP/(TP+FN)

        push!(class_precision,precision)
        push!(all_tp_2,TP)
        push!(all_fn,FN)
        push!(all_class_pred,labels_c[(restr_pred .== predict_label_to_vertex[n])])
        push!(all_persist_pred,restr_persist[restr_pred .== predict_label_to_vertex[n]])
        if !null
            push!(all_mutant_pred,restr_mut[restr_pred .== predict_label_to_vertex[n]])
        end
    end

    class_precision = [isnan(p) ? 0. : p for p in class_precision]

    @assert all(all_tp_1 .== all_tp_2)

    recall_micro = sum(all_tp_1)/(sum(all_tp_1) + sum(all_fp)) # equiv. to total accuracy
    recall_macro = mean(class_recall)

    precision_micro = sum(all_tp_1)/(sum(all_tp_1) + sum(all_fn)) 
    precision_macro = mean(class_precision)

    f1_class = 2 .* (class_precision .* class_recall) ./ (class_precision .+ class_recall) 
    f1_micro = 2 * (recall_micro * precision_micro) / (recall_micro + precision_micro)
    f1_macro = 2 * (recall_macro * precision_macro) / (recall_macro + precision_macro)

    f1_class = [isnan(p) ? 0. : p for p in f1_class]

    return recall_micro,f1_micro,f1_macro,f1_class,class_recall,av_we,all_pred_class, all_persist_class,all_class_pred,all_persist_pred,all_mutant_pred
end    

function create_label_H0_parent(tr,top_vertex_map,predict_id,predict_mst,predict_mst_complexity)
    if tr.minimal_stripe_subgraphs[tr.H0] ∈ predict_mst
        top_vertex_map[tr.minimal_stripe_subgraphs[tr.H0]]
    else
        incl = [test_inclusion(tr.minimal_stripe_subgraphs[tr.H0],pmst) for pmst in predict_mst]

        if sum(incl) != 0
            inc_id = findall(x->x==1, incl)
            choice = inc_id[argmax(predict_mst_complexity[inc_id])]
            predict_id[choice]
        else
            -1
        end

    end
end

function create_label_H0(tr,top_vertex_map,predict_id,predict_mst)
    if tr.minimal_stripe_subgraphs[tr.H0] ∈ predict_mst
        top_vertex_map[tr.minimal_stripe_subgraphs[tr.H0]]
    else
        -1
    end
end

####### Phenotype Evaluations 

function new_phenotype(ind::Individual,new_w::Matrix{Float64},development)
    Individual(remake(ind.genotype, p = (new_w,ind.genotype.p[2:end]...)),development).phenotype.u[end][3,:]
end

function pheno_characterise_v1(conc,lower_bound,upper_bound,min_stripe_width)

    low = findall(conc .< lower_bound)
    high = findall(conc .> upper_bound)

    left_transition = nothing
    right_transition = nothing

    max_height = maximum(conc)
    min_height = minimum(conc)

    if (length(low) != 0) & (length(high) != 0)

        if low[1] < high[1] # first low cell is before first high cell --> left boundary
            if high[1] - low[1] >= min_stripe_width
                left_transition = high[1] # true
            end

            if low[end] > high[end]
                if (high[end] - high[1] + 1 >= min_stripe_width) & (low[end] - high[end] >= min_stripe_width)
                    cts_high_region = !(any([(id > high[1]) & (id < high[end]) for id in low]))
                    if cts_high_region
                        right_transition = high[end] # true
                    end
                end
            end
        else # first high cell is before first low cell --> right boundary, could be left transition but would be inverse
            if (low[1] - high[1] >= min_stripe_width)
                right_transition = low[1]
            end
        end

    end

    return (left_transition,right_transition,max_height,min_height)

end

function characterise_mutation_v1(ph_profile_1,ph_profile_2,neutral_ind)

    if neutral_ind
        mutant_profile = [:neutral]
    else
        mutant_profile = [:unclassified]

        left_boundary_was_present  = ph_profile_1[1] != nothing
        right_boundary_was_present = ph_profile_1[2] != nothing

        left_boundary_now_present = ph_profile_2[1] != nothing
        right_boundary_now_present = ph_profile_2[2] != nothing


        if left_boundary_now_present & (! left_boundary_was_present)

            if right_boundary_was_present 
                push!(mutant_profile, :clb_wrb)
            else
                push!(mutant_profile, :clb)
            end

        elseif left_boundary_now_present & left_boundary_was_present

            move_boundary = ph_profile_1[1] != ph_profile_2[1]

            if move_boundary
                if right_boundary_now_present
                    push!(mutant_profile, :mlb_wbb)
                else
                    push!(mutant_profile, :mlb_wlb)
                end
            end

        elseif (! left_boundary_now_present) & left_boundary_was_present

            push!(mutant_profile, :dlb)
        else
            nothing
        end

        #############

        if right_boundary_now_present & (! right_boundary_was_present)

            if left_boundary_was_present 
                push!(mutant_profile, :crb_wlb)
            else
                push!(mutant_profile, :crb)
            end

        elseif right_boundary_now_present & right_boundary_was_present

            move_boundary = ph_profile_1[2] != ph_profile_2[2]

            if move_boundary
                if left_boundary_now_present
                    push!(mutant_profile, :mrb_wbb)
                else
                    push!(mutant_profile, :mrb_wrb)
                end
            end

        elseif (! right_boundary_now_present) & right_boundary_was_present

            push!(mutant_profile, :drb)
        else
            nothing
        end

        if length(mutant_profile) == 1
            if (ph_profile_2[3] > ph_profile_1[3]) | (ph_profile_2[4] < ph_profile_1[4])
                if right_boundary_was_present & left_boundary_was_present
                    push!(mutant_profile,:grow_wbb)
                elseif right_boundary_was_present
                    push!(mutant_profile,:grow_wrb)
                elseif left_boundary_was_present
                    push!(mutant_profile,:grow_wlb)
                else
                    push!(mutant_profile,:grow_wnb)
                end
            end
        end


        if (:mlb ∈ mutant_profile) & (:mrb ∈ mutant_profile)
            mutant_profile = [:nothing,:mbb_wbb]
        end

        if (:clb ∈ mutant_profile) & (:crb ∈ mutant_profile)
            mutant_profile = [:nothing,:cbb]
        end

        if (:dlb ∈ mutant_profile) & (:crb_wlb ∈ mutant_profile)
            mutant_profile = [:nothing,:switch_lr]
        end

        if (:drb ∈ mutant_profile) & (:clb_wrb ∈ mutant_profile)
            mutant_profile = [:nothing,:switch_rl]
        end
    end

    if length(mutant_profile) > 1
        return mutant_profile[2:end]
    else
        return [mutant_profile[1]]
    end
end

function get_first_boundary(tr)
    lb = findall(x->x==[:clb],tr)
    rb = findall(x->x==[:crb],tr)
    bb = findall(x->x==[:cbb],tr)

    v_list = []
    for v in [lb,rb,bb]
        if length(v) > 0
            push!(v_list,minimum(v))
        else
            push!(v_list,Inf)
        end
    end

    return argmin(v_list)
end

function replace_nothing_phenotypes(ph_list)

    new_ph_list = []
    
    last_p = nothing

    for p in ph_list 
        if !isnothing(p)
            last_p = p
            push!(new_ph_list,p)
        else
            push!(new_ph_list,last_p)
        end
    end

    return new_ph_list
end

####### Convenience

function cond_save(dir,fig,cond)
    if cond
        CairoMakie.save(dir,fig,pt_per_unit = 1)
    end
end

