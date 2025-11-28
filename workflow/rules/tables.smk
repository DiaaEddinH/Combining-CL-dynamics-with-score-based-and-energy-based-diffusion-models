OBS = ["cumulants", "moments"]

HEADERS = {
    "QMODEL_cumulants": [
        *[f"Re $\\kappa_{2*(i+1)}$" for i in range(4)],
        *[f"Im $\\kappa_{2*(i+1)}$" for i in range(4)],
    ],
    "QMODEL_moments": [
        *[f"Re $\mu_{2*(i+1)}$" for i in range(4)],
        *[f"Im $\mu_{2*(i+1)}$" for i in range(4)],
    ],
	"GMODEL_cumulants": [
        *[f"$<x^{2*(i+1)}>_c$" for i in range(4)],
        *[f"$<y^{2*(i+1)}>_c$" for i in range(4)],
    ],
    "GMODEL_moments": [
        *[f"$<x^{2*(i+1)}>$" for i in range(4)],
        *[f"$<y^{2*(i+1)}>$" for i in range(4)],
    ],
}

rule build_tables:
	input:
		# npz="data/processed/{experiment}_EBM_{obs}_random_effects.npz",
		npz=expand("data/processed/{experiment}_{model}_{obs}_random_effects.npz", experiment=EXPERIMENT, obs=OBS, model=["SBM", "EBM", "EBMMH"])
	output:
		tex="assets/tables/{experiment}_{obs}.tex"
	params:
		experiment=lambda wc: wc.experiment,
		headers=lambda wc: HEADERS[f"{wc.experiment}_{wc.obs}"]
	conda: "../envs/analysis.yml"
	script: "../scripts/run_tables.py"
