# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import numpy as np

rgb2hex = lambda r, g, b: "#%02x%02x%02x" % (int(r * 256), int(g * 256), int(b * 256))

OFE_COLORS = (
    (49 / 256, 57 / 256, 77 / 256),  # Badass Blue
    (184 / 256, 87 / 256, 65 / 256),  # Feeling spicy
    (0, 147 / 256, 132 / 256),  # Feeling sick
    (217 / 256, 196 / 256, 177 / 256),  # Beastly grey
    (217 / 256, 196 / 256, 177 / 256),  # Sandy Sergio
    (238 / 256, 192 / 256, 68 / 256),  # Gold
    (0 / 256, 47 / 256, 74 / 256),  # otherBlue
)


def color_gradient(
    c1: tuple[float, float, float] = OFE_COLORS[1],
    c2: tuple[float, float, float] = OFE_COLORS[2],
    c3: tuple[float, float, float] = OFE_COLORS[1],
    mix: float = 0,
    hex: bool = True,
):
    c1 = np.array(c1)
    c2 = np.array(c2)
    c3 = np.array(c3)
    mix = np.array(mix, ndmin=1)

    if mix > 0.5:
        m = mix - 0.5
        c = (0.5 - m) * c2 + m * c3
    else:
        m = mix
        c = (0.5 - m) * c1 + m * c2

    if hex:
        return rgb2hex(*c)
    else:
        return c
