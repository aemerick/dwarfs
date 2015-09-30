from __future__ import division

import numpy as np


def generate_particle_distribution(DF, N):
    """
    Given a distribution function (df) object, uses the accept-reject method to populate the halo
    with particles
    """