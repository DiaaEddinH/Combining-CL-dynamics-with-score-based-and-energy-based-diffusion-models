from glob import glob
from pathlib import Path

# Discover all sample stems
SAMPLES = [
    Path(f).stem.replace("_samples", "")
    for f in glob(f"data/samples/{EXPERIMENT}_*_samples.npy")
]


rule compute_moments:
    input: lambda wc: f"data/samples/{wc.sample}_samples.npy"
    output: "data/processed/{sample}_moments.npz"
    threads: 2
    resources: mem_mb=2000
    conda: "../envs/analysis.yml"
    script: "../../libs/DiffusionModels/scripts/compute_moments.py"

rule compute_other_moments:
	input: lambda wc: f"data/samples/{wc.sample}_samples.npy"
	output: "data/processed/{sample}_other_moments.npz"
	threads: 2
	resources: mem_mb=2000
	conda: "../envs/analysis.yml"
	script: "../../libs/DiffusionModels/scripts/compute_other_moments.py"

rule compute_cumulants:
    input: lambda wc: f"data/samples/{wc.sample}_samples.npy"
    output: "data/processed/{sample}_cumulants.npz"
    threads: 2
    resources: mem_mb=2000
    conda: "../envs/analysis.yml"
    script: "../../libs/DiffusionModels/scripts/compute_cumulants.py"
