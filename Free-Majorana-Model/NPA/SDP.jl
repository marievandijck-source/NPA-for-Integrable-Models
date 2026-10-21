function build_moment_sdp(
    basis::Vector{Word},
    H::Operator,
    optimizer,
)

    # --------------------------------------------------
    # Moment Matrix structure and content
    # --------------------------------------------------

    structure = moment_matrix_structure(basis)

    moments = moment_words(structure)

    word_to_index = Dict(
        word => index
        for (index, word) in enumerate(moments)
    )

    # --------------------------------------------------
    # Every Hamiltonian word must occur in the relaxation
    # --------------------------------------------------

    missing_words = [
        word
        for word in keys(H.coeffs)
        if !haskey(word_to_index, word)
    ]

    if !isempty(missing_words)
        error(
            "The following Hamiltonian words are absent " *
            "from the moment relaxation: $missing_words. " *
            "Increase the basis degree.",
        )
    end

    # --------------------------------------------------
    # JuMP model
    # --------------------------------------------------

    model = Model(optimizer)

    number_of_moments = length(moments)
    matrix_dimension = length(basis)

    # --------------------------------------------------
    # Complex variable y_W = ⟨W⟩ per canonical word
    # --------------------------------------------------
    @variable(
        model,
        y[1:number_of_moments] in ComplexPlane(),
    )

    # --------------------------------------------------
    # Hermitian PSD moment matrix
    # --------------------------------------------------
    @variable(
        model,
        Γ[1:matrix_dimension, 1:matrix_dimension]
        in HermitianPSDCone(),
    )

    # --------------------------------------------------
    # Connect moment matrix to the moment variables
    # --------------------------------------------------
    for i in 1:matrix_dimension
        for j in 1:matrix_dimension
            entry = structure[i, j]
            moment_index = word_to_index[entry.word]

            @constraint(
                model,
                Γ[i, j] ==
                entry.coefficient * y[moment_index],
            )
        end
    end

    # --------------------------------------------------
    # State normalization: ⟨I⟩ = 1
    # --------------------------------------------------
    identity_word = Word()
    identity_index = word_to_index[identity_word]

    @constraint(
        model,
        y[identity_index] == 1.0 + 0.0im,
    )

    # --------------------------------------------------
    # Adjoint consistency : ⟨W†⟩ = conj(⟨W⟩)
    # --------------------------------------------------
    for word in moments
        adjoint_sign, adjoint_canonical_word =
            adjoint(word)

        word_index = word_to_index[word]

        if !haskey(
            word_to_index,
            adjoint_canonical_word,
        )
            error(
                "Adjoint word $adjoint_canonical_word " *
                "is absent from the moment list.",
            )
        end

        adjoint_index =
            word_to_index[adjoint_canonical_word]

        @constraint(
            model,
            adjoint_sign * y[adjoint_index] ==
            conj(y[word_index]),
        )
    end

    # --------------------------------------------------
    # Hamiltonian expectation value
    # --------------------------------------------------

    energy = sum(
        coefficient * y[word_to_index[word]]
        for (word, coefficient) in H.coeffs
    )

   
    # --------------------------------------------------
    # Hermitian Hamiltonian has a real expectation value
    # --------------------------------------------------
    @constraint(model, imag(energy) == 0.0)

    @objective(
        model,
        Min,
        real(energy),
    )

    return (
        model=model,
        moments=moments,
        word_to_index=word_to_index,
        y=y,
        Γ=Γ,
        structure=structure,
        energy=energy,
    )
end


function solve_moment_sdp!(
    problem;
    verbose::Bool=true,
)
    model = problem.model

    if !verbose
        set_silent(model)
    end

    optimize!(model)

    status = termination_status(model)
    primal_status_value = primal_status(model)
    dual_status_value = dual_status(model)

    successful_statuses = (
        MOI.OPTIMAL,
        MOI.ALMOST_OPTIMAL,
        MOI.LOCALLY_SOLVED,
    )
    
    if !(status in successful_statuses)
        return (
            success=false,
            termination_status=status,
            primal_status=primal_status_value,
            dual_status=dual_status_value,
            objective=nothing,
            moments=nothing,
        )
    end

    moment_values = Dict(
        word => value(problem.y[index])
        for (index, word) in
            enumerate(problem.moments)
    )

    return (
        success=true,
        termination_status=status,
        primal_status=primal_status_value,
        dual_status=dual_status_value,
        objective=objective_value(model),
        moments=moment_values,
    )
end