from dataclasses import dataclass
from sage.all import Graph

@dataclass(slots=True)
class Grafo2Terminales:
    grafo: Graph
    s: int
    t: int

    def __post_init__(self):
        if self.s not in self.grafo:
            raise ValueError(f"s={self.s} no es un vértice del grafo")
        if self.t not in self.grafo:
            raise ValueError(f"t={self.t} no es un vértice del grafo")

    def __eq__(self, other):
        if not isinstance(other, Grafo2Terminales):
            return NotImplemented
        if (self.s, self.t) != (other.s, other.t):
            return False
        return self.grafo.is_isomorphic(other.grafo)

    def show(self, **kwargs):
        vertex_colors = {
            "lime": [self.s, self.t],
            "gold": [v for v in self.grafo.vertices() if v not in (self.s, self.t)],
        }
        self.grafo.show(vertex_colors=vertex_colors, **kwargs)