import torch
import numpy as np
import importlib
from model_setup import setup_model
from diffusion_models.utils import detach_to_numpy
from diffusion_models.datasets.datasets import QuarticCL, DoublePeak

EXPERIMENT_DATASET = {
    "QMODEL": QuarticCL(),
    "GMODEL": DoublePeak(mu=np.array([1, -1]), sigma=0.25, size=1_000_000),
}

EXPERIMENT_NETWORK = {"QMODEL": "EvenResNet", "GMODEL": "EvenLinear"}


def get_network_cls(name: str):
    module = importlib.import_module("diffusion_models.networks.networks")
    return getattr(module, name)


def save_model_drift_data(
    model,
    filename,
    time_point=1e-3,
    no_of_points=101,
    x_lims=(-3, 3),
    y_lims=(-3, 3),
    rescale=(1, 1),
):
    device = model.device
    a_x, a_y = rescale

    y = np.linspace(*x_lims, no_of_points)
    x = np.linspace(*y_lims, no_of_points)
    gridx, gridy = np.meshgrid(x, y)
    X = np.stack([gridx, gridy], axis=-1)

    requires_grad = "EnergyBasedModel" in str(type(model))
    t0 = torch.tensor(time_point, device=device)
    tensor_X = torch.tensor(
        X, dtype=torch.float32, device=device, requires_grad=requires_grad
    )

    v = model.network(tensor_X, t0)
    logP = 0.5 * torch.sum(v**2, dim=-1) / model.schedule.stddev(t0)
    logP = detach_to_numpy(logP)
    logP -= np.min(logP)

    P = np.exp(-logP)
    norm = np.sum(P) * (y[1] - y[0]) * (x[1] - x[0]) + 1e-12
    P /= norm

    if requires_grad:
        score = detach_to_numpy(model(tensor_X, t0)).squeeze()
    else:
        score = detach_to_numpy(v) / model.schedule.stddev(t0)

    score[..., 0] /= a_y
    score[..., 1] /= a_x

    np.savez(
        file=filename,
        gridx=a_x * gridx,
        gridy=a_y * gridy,
        score=score,
        logP=logP,
        pdf=P / (a_x * a_y),
    )


if __name__ == "__main__":
    config = snakemake.input.config
    output = snakemake.output[0]
    network_name = (
        EXPERIMENT_NETWORK[snakemake.params.experiment]
        if snakemake.params.model == "EBM"
        else "OddLinear"
    )
    if snakemake.params.experiment == "QMODEL" and snakemake.params.model == "SBM":
        network_name = "LinearNet"
    network_cls = get_network_cls(network_name)
    print(snakemake.params.experiment, snakemake.params.model, network_cls)
    model = setup_model(config_filepath=config, network_cls=network_cls, device="cpu")
    sig = EXPERIMENT_DATASET[snakemake.params.experiment].data.std(0)
    # np.loadtxt("data/raw/cl_K111_ccc.dat", delimiter=",", dtype=np.float32).std(0)
    save_model_drift_data(model, filename=output, rescale=tuple(sig))
