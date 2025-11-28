import pandas as pd
import numpy as np
from format_multiple_errors import format_column_errors


def build_table(experiment: str, obs: str):
    labels = ["SBM", "EBM", "EBMMH"]

    out_dict = {}

    for label in labels:
        npz = np.load(f"data/processed/{experiment}_{label}_{obs}_random_effects.npz")
        raw = dict({item: npz[item] for item in npz.files}, orient="index")

        df = pd.concat(
            {k: pd.DataFrame(v) for k, v in raw.items() if isinstance(v, np.ndarray)},
            axis=1,
        )

        df0 = df.xs(0, level=1, axis=1).drop("sigma_tot", axis=1)
        df1 = df.xs(1, level=1, axis=1).drop("sigma_tot", axis=1)

        df0["x"] = format_column_errors(df0, abbreviate=True)
        df1["y"] = format_column_errors(df1, abbreviate=True)

        x = df0["x"].reset_index(drop=True)[1::2]
        y = df1["y"].reset_index(drop=True)[1::2]

        x.index = [f"f0_{2*(i+1)}" for i in range(len(x))]
        y.index = [f"f1_{2*(i+1)}" for i in range(len(y))]

        flat_df = pd.DataFrame([pd.concat([x, y])]).to_dict(orient="records")[0]

        out_dict[label] = flat_df

    return pd.DataFrame.from_dict(out_dict, orient="index")


if __name__ == "__main__":
    obs = snakemake.wildcards.obs
    experiment = snakemake.params.experiment
    headers = snakemake.params.headers

    df = build_table(experiment, obs)
    df.to_latex(header=headers, buf=snakemake.output.tex)
