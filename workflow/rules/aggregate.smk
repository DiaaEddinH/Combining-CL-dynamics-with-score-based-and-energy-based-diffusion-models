# Define a generic aggregation rule
def aggregate_npz(experiment, obs):
	exp_samples = [s for s in SAMPLES if s.startswith(f"{experiment}")]
	return expand(
		"data/processed/{s}_{obs}.npz", s=exp_samples, obs=obs
	)

rule random_effects:
	input:
		lambda wc: aggregate_npz(wc.experiment, wc.obs)
	output:
		"data/processed/{experiment}_{obs}_random_effects.npz"
	threads: 2
	resources: mem_mb=2000
	conda: "../envs/analysis.yml"
	script: "../../libs/DiffusionModels/scripts/compute_random_effects.py"
