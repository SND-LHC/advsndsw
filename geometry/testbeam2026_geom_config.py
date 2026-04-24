"""Configuration parameters for AdvSND geometry."""

import shipunit as u
from ShipGeoConfig import AttrDict, ConfigRegistry

if "tb_2026_mc" in globals():
    tb_2026_mc = globals()["tb_2026_mc"]
else:
    tb_2026_mc = False
if "n_target_layers" in globals():
    n_target_layers = globals()["n_target_layers"]
else:
    n_target_layers = 20

with ConfigRegistry.register_config("basic") as c:
    # cave parameters
    c.cave = AttrDict(z=0 * u.cm)

    c.Floor = AttrDict(z=48000.*u.cm)

    # AdvSND Target & Tracker structure
    c.AdvTarget = AttrDict(z=0 * u.cm)
    c.AdvTarget.testbeam2026 = tb_2026_mc
    # Target plates: altrenating 2-mm Fe, 3x1-mm W, 2-mm Fe
    # Tungsten plates
    c.AdvTarget.WPlateX = 19.2 * u.cm
    c.AdvTarget.WPlateY = c.AdvTarget.WPlateX
    c.AdvTarget.nWPlates = 3
    c.AdvTarget.WPlateZ = 1 * u.mm
    # Iron plates
    c.AdvTarget.FePlateX = 25 * u.cm
    c.AdvTarget.FePlateY = c.AdvTarget.FePlateX
    c.AdvTarget.FePlateZ = 2 * u.mm
    # Support frame of the plates
    c.AdvTarget.PlateFrameX = 59 * u.cm
    c.AdvTarget.PlateFrameY = 59 * u.cm
    c.AdvTarget.PlateFrameZ = 3 * u.mm
    c.AdvTarget.FrameBorderX = 7.25 * u.cm
    c.AdvTarget.FrameBorderY = c.AdvTarget.FrameBorderX

    # Target Tracking layers
    c.AdvTarget.TTZ = (
        4 * u.mm
    )  # due to module overlaps, the thickness of a layer is actually 2x4mm
    c.AdvTarget.nTT = n_target_layers

    # Full target volume including inactive area
    c.AdvTarget.X = 59 * u.cm
    c.AdvTarget.Y = c.AdvTarget.X
    c.AdvTarget.Z = c.AdvTarget.nTT * (2 * c.AdvTarget.TTZ
                                       + c.AdvTarget.FePlateZ
                                       + c.AdvTarget.nWPlates * c.AdvTarget.WPlateZ
                                       + c.AdvTarget.FePlateZ
                                       + c.AdvTarget.PlateFrameZ)
