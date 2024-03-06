import numpy as np

# Store atom name and 3d coordinates
class Atom:
    def __init__(self, name = "", coord = np.array(3)):
        self.name = name
        self.coord = coord
    def __lt__(self, other):
        return self.name < other.name
    def __repr__(self):
        return f"Atom({self.name}, {self.coord})"