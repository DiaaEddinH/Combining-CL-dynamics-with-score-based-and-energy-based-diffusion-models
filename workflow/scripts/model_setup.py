import torch
import numpy as np
from pathlib import Path
from diffusion_models.config.config_loader import load_config
from diffusion_models.noise.noise_scheduler import GeometricSchedule
from diffusion_models.models.models import ScoreModel, EnergyBasedModel
from diffusion_models.utils import set_device, get_activation_func


def setup_model(config_filepath: str, network_cls, device=None):
    weight_dir = Path("data/weights")

    config = load_config(config_path=config_filepath)

    model_cls = EnergyBasedModel if "EBM" in config["file"] else ScoreModel
    schedule = GeometricSchedule(
        sigma_min=config["sigma_min"], sigma_max=config["sigma_max"]
    )

    weight_file = weight_dir / f"{config['file']}_weights.pt"

    device = device or config["device"]
    config["device"] = set_device(device)
    config["channels"] = config["hidden_channels"]
    config["activation"] = get_activation_func(config["activation"])()

    network = network_cls(**config)
    model = model_cls(network, schedule, config["device"])

    model._load_weights(weight_file)
    model.eval()
    model.ema.apply_shadow()

    return model
