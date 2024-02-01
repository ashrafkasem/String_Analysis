import pandas as pd
import numpy as np
from typing import Type
from core.policies.membrane_policy import Membrane
class String:
    def __init__(self, membrane : Type[Membrane],
                 radius: int = 7.5,
                 distance: int = 25,
                 gap: int = 35):
        self.membrane = membrane
        self.psii_strings = self.find_psii_strings()
        self.radius = radius
        self.distance = distance 
        self.gap = gap
    def find_psii_strings(self):
        # Implement your PSII string analysis logic here
        # You can access membrane.psii_positions to get the PSII positions
        # Apply your algorithm to identify strings and return the result
        # For example, you might want to store the strings as a list of lists of PSII positions
        psii_strings = []  # Placeholder, modify based on your logic
        return psii_strings