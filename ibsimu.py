import math
from _ibsimu import *
import ezdxf


class Environment:
    def __init__(
        self,
        path,
        mesh=1e-3,
        size_x=150.0e-3,
        size_y=150.0e-3,
        size_z=2500.0e-3,
        z0=-40.0e-3,
        exclude=None,
    ):
        self.path = path
        self.solids = MyDXFFile(path)
        self.solids.set_warning_level(2)
        e = self.solids.get_entities()
        e.scale(e.selection_all(), self.solids, 1.0e-3)

        self.layers = self.get_layer_names(path, exclude=exclude)

        self.layer_ids = {layer: i for i, layer in enumerate(self.layers, 7)}

        origin = Vec3D(-size_x / 2, -size_y / 2, z0)
        size = Int3D(
            math.floor(size_x / mesh) + 1,
            math.floor(size_y / mesh) + 1,
            math.floor(size_z - z0 / mesh) + 1,
        )
        mode = GeometryMode.MODE_3D
        self.geometry = Geometry(mode, size, origin, mesh)

        self.set_defaults()

    def set_defaults(self):
        self.layer_ids = {}
        for layer_id, layer in enumerate(self.layers, 7):
            solid = DXFSolid(self.solids, layer)
            solid.cylindric()

            self.layer_ids[layer] = layer_id
            self.geometry.set_solid(layer_id, solid)
            potential = 0.0
            self.set_potential(layer, potential, BoundaryTypes.DIRICHLET)

    def set_potential(self, layer, potential, boundary_type):
        if isinstance(layer, str):
            layer = self.layer_ids[layer]

        self.geometry.set_boundary(layer, Bound(boundary_type, potential))

    def get_layer_names(self, path, exclude=None):
        if exclude is None:
            exclude = ()

        if isinstance(exclude, str):
            exclude = (exclude,)

        exclude = exclude + ("0", "Defpoints")

        doc = ezdxf.readfile(path)
        names = [layer.dxf.name for layer in doc.layers]
        names = [name for name in names if name not in exclude]
        return names


if __name__ == "__main__":
    exclude = "Einzel1_left_ground_shield"
    env = Environment("PIA_geom_V01.dxf", exclude=exclude)

