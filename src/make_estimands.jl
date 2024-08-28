
using TMLE
using YAML

# This was obtained by running: scripts/variants_ukb_info.jl
const VARIANTS_INFO = Dict(
    "rs1805005"   => (["GG", "GT", "TT"], 0.119632),
    "rs974766"    => (["AA", "AC", "CC"], 0.082168),
    "rs9940128"   => (["GG", "GA", "AA"], 0.421814),
    "rs456998"    => (["GG", "GT", "TT"], 0.48168),
    "rs6059655"   => (["AA", "AG", "GG"], 0.0985599),
    "rs10132320"  => (["AA", "AG", "GG"], 0.00658415),
    "rs1129038"   => (["CC", "CT", "TT"], 0.257224),
    "rs9926016"   => (["GG", "GA", "AA"], 0.136166),
    "rs1805008"   => (["CC", "CT", "TT"], 0.0836136),
    "rs72926466"  => (["GG", "GC", "CC"], 0.124016),
    "rs4506565"   => (["AA", "AT", "TT"], 0.313938),
    "rs11868112"  => (["CC", "CT", "TT"], 0.409067),
    "rs1732170"   => (["TT", "TC", "CC"], 0.3964),
    "rs148515035" => (["TT", "TA", "AA"], 0.0188177),
    "rs6088372"   => (["CC", "CT", "TT"], 0.155658),
    "rs6456121"   => (["CC", "CT", "TT"], 0.295634),
    "rs62107261"  => (["TT", "TC", "CC"], 0.0463914),
    "rs117737810" => (["TT", "TC", "CC"], 0.0307375),
    "rs59103106"  => (["AA", "AC", "CC"], 0.034453),
    "rs3859191"   => (["GG", "GA", "AA"], 0.46476),
    "rs10419224"  => (["CC", "CT", "TT"], 0.0112278),
    "rs3129889"   => (["GG", "GA", "AA"], 0.139662),
    "rs62295911"  => (["GG", "GA", "AA"], 0.0564486),
    "rs502771"    => (["TT", "TC", "CC"], 0.272973),
    "rs9268219"   => (["TT", "TG", "GG"], 0.118571),
    "rs3129716"   => (["TT", "TC", "CC"], 0.143209),
    "rs8111699"   => (["CC", "CG", "GG"], 0.478006),
    "rs1805007"   => (["CC", "CT", "TT"], 0.0977936),
    "rs356219"    => (["GG", "GA", "AA"], 0.371947)
)

const CONFOUNDERS = [Symbol(:PC, i) for i in 1:6]

const EXTRA_COVARIATES = ["Genetic-Sex", "Age-Assessment"]

function make_treatment_levels(variants)
    variant_genotypes = [VARIANTS_INFO[string(v)][1] for v in variants]
    return NamedTuple{variants}(variant_genotypes)
end

MakeATE(treatment_levels, outcome) = factorialEstimand(ATE, treatment_levels, outcome;
    confounders=CONFOUNDERS, 
    outcome_extra_covariates=EXTRA_COVARIATES
    )

MakeAIE(treatment_levels, outcome) = factorialEstimand(AIE, treatment_levels, outcome;
    confounders=CONFOUNDERS, 
    outcome_extra_covariates=EXTRA_COVARIATES
    )

function geneatlas_traits_to_variants()
    return Dict(
        # White blood cell counts: Count
        "White blood cell (leukocyte) count" => (:rs3859191, :rs9268219),
        # Sarcoidosis : rare ~ 500 cases
        "sarcoidosis" => (:rs502771, :rs148515035), 
        # MS: ~ 1400 cases 
        "G35 Multiple sclerosis" => (:rs3129889, :rs62295911),
        # Other digestive diseases: ~ 13000 cases
        "K90-K93 Other diseases of the digestive system" => (:rs3129716, :rs72926466), 
        # BMI: continuous
        "Body mass index (BMI)" => (:rs9940128, :rs62107261),
        # T2D ~ 2800 cases
        "type 2 diabetes" => (:rs4506565, :rs117737810)
    )
end

function epistatic_traits_to_variants()
    return Dict(
        # Skin colour: https://www.nature.com/articles/s41467-018-07691-z
        "Skin colour" => [(:rs1805007, :rs6088372), (:rs1805005, :rs6059655), (:rs1805008, :rs1129038)],
        # Parkison: https://pubmed.ncbi.nlm.nih.gov/31234232/
        "G20 Parkinson's disease" => [(:rs1732170, :rs456998, :rs356219, :rs8111699), (:rs11868112, :rs6456121, :rs356219)],
        # MS: https://www.sciencedirect.com/science/article/pii/S0002929723000915
        "G35 Multiple sclerosis" => [(:rs10419224, :rs59103106)],
        # Psoriasis: https://www.sciencedirect.com/science/article/pii/S0002929723000915
        "psoriasis" => [(:rs974766, :rs10132320)],
    )
end

function MakeEstimands()
    estimands = []
    # All single variant effects identified trough the geneATLAS
    for (outcome, variants) ∈ geneatlas_traits_to_variants()
        # Make Single variant effects reported by geneATLAS
        for variant in variants
            treatment_levels = make_treatment_levels((variant,))
            push!(estimands, MakeATE(treatment_levels, outcome))
        end
        # Pairwise Interactions, likely to be non present
        treatment_levels = make_treatment_levels(variants)
        push!(estimands, MakeAIE(treatment_levels, outcome))
    end

    # All Epistatic Interactions identified through various sources
    for (outcome, variants_tuples) ∈ epistatic_traits_to_variants()
        for variants in variants_tuples
            treatment_levels = make_treatment_levels(variants)
            push!(estimands, MakeAIE(treatment_levels, outcome))
        end
    end

    return estimands
end

function SaveEstimands(;output="assets/estimands.yaml")
    estimands = MakeEstimands()
    TMLE.write_yaml(output, TMLE.Configuration(estimands=estimands))
end

SaveEstimands()