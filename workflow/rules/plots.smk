# plots.smk — Drift field and potential plotting pipeline

from pathlib import Path

# User must specify NETWORK via --config NETWORK=EvenResNet | EvenLinear
NETWORK = config.get("NETWORK", "EvenResNet")   # default

PLOTS = [
    "potential",
    "distribution",
    "streamplot",
    "streamplot_w_thimbles",
]

def limited_plots(experiment, model):
	out = ["streamplot"]
	if model == "EBM":
		out.extend(["potential", "distribution"])
	if experiment == "QMODEL":
		out.append("streamplot_w_thimbles")
	return out

ALL_PLOTS = [
    f"assets/plots/{EXPERIMENT}_{model}_{kind}.pdf"
    for model in MODEL
    for kind in limited_plots(EXPERIMENT, model)
]

print(ALL_PLOTS)


# ---------------------------------------------------------------------
# Rule: Analytical CL drift (depends only on given parameters)
# ---------------------------------------------------------------------

rule plot_analytical_drift:
	output:
		pdf_cl="assets/plots/complex_langevin_drift.pdf",
		pdf_gaussian="assets/plots/gaussian_double_peak_drift.pdf"
	params:
		mass="1+1j",
		coupling=1,
		mu=(1, -1),
		sigma=0.25,
		x_lims=(-2, 2),
		y_lims=(-0.5, 0.5),
	conda: "../envs/analysis.yml"
	script: "../scripts/plot_analytical_drift.py"


# ---------------------------------------------------------------------
# Rule: Save drift NPZ for trained model
# ---------------------------------------------------------------------

rule save_model_drift_npz:
	input:
		config="configs/{experiment}_{model}_config.yaml"
	output:
		npz="data/processed/{experiment}_{model}_drift.npz"
	wildcard_constraints:
		model="EBM|SBM"
	params:
		network=NETWORK,
		experiment=lambda wc: wc.experiment,
		model=lambda wc: wc.model
	conda: "../envs/analysis.yml"
	script: "../scripts/save_model_drift_data.py"


# ---------------------------------------------------------------------
# Rule: Plot NPZ drift (potential, distribution, streamplot, thimbles)
# ---------------------------------------------------------------------

rule plot_model_drift:
	input: lambda wc: f"data/processed/{wc.experiment}_{wc.model}_drift.npz"
	output: "assets/plots/{experiment}_{model}_{kind}.pdf"
	wildcard_constraints:
		model="EBM|SBM"
	params: prefix=lambda wc: f"assets/plots/{wc.experiment}_{wc.model}"
	conda: "../envs/analysis.yml"
	script: "../scripts/plot_drift_data.py"
