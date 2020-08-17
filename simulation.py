from sys import argv
from ibsimu import *

set_message()

mesh = 1e-3
lengthX = 150.0e-3
lengthY = 150.0e-3
lengthZ = 25.0e-3
startZ = -40.0e-3

# beamIntensity_mA = float(argv[1]) * 1e-3
V_Einzel1 = float(argv[1]) * 1.0e3
V_Einzel2 = float(argv[2]) * 1.0e3
V_Einzel3 = float(argv[3]) * 1.0e3

origin = Vec3D(-lengthX / 2, -lengthY / 2, startZ)
size = Int3D(
    math.floor(lengthX / mesh) + 1,
    math.floor(lengthY / mesh) + 1,
    math.floor((lengthZ - startZ) / mesh) + 1,
)

mode = GeometryMode.MODE_3D
geom = Geometry(mode, size, origin, mesh)

# // Define solids
geometry_file = "PIA_geom_V01.dxf"
dxfPIAlinear = MyDXFFile(geometry_file)
dxfPIAlinear.set_warning_level(2)
e = dxfPIAlinear.get_entities()
sel = e.selection_all()
e.scale(sel, dxfPIAlinear, 1.0e-3)

plasma = DXFSolid(dxfPIAlinear, "Plasma")
plasma.define_2x3_mapping(DXFSolid.rotz)  # alternatively: plasma.cylindric() also works
geom.set_solid(7, plasma)


Einzel1_left_ground = DXFSolid(dxfPIAlinear, "Einzel1_left_ground")
Einzel1_left_ground.define_2x3_mapping(DXFSolid.rotz)
geom.set_solid(8, Einzel1_left_ground)


Einzel1 = DXFSolid(dxfPIAlinear, "Einzel1")
Einzel1.define_2x3_mapping(DXFSolid.rotz)
geom.set_solid(9, Einzel1)

Einzel1_right_ground = DXFSolid(dxfPIAlinear, "Einzel1_right_ground")
Einzel1_right_ground.define_2x3_mapping(DXFSolid.rotz)
geom.set_solid(10, Einzel1_right_ground)

Einzel2_left_ground = DXFSolid(dxfPIAlinear, "Einzel2_left_ground")
Einzel2_left_ground.define_2x3_mapping(DXFSolid.rotz)
geom.set_solid(11, Einzel2_left_ground)

Einzel2 = DXFSolid(dxfPIAlinear, "Einzel2")
Einzel2.define_2x3_mapping(DXFSolid.rotz)
geom.set_solid(12, Einzel2)

Einzel2_right_ground = DXFSolid(dxfPIAlinear, "Einzel2_right_ground")
Einzel2_right_ground.define_2x3_mapping(DXFSolid.rotz)
geom.set_solid(13, Einzel2_right_ground)

Einzel3_left_ground = DXFSolid(dxfPIAlinear, "Einzel3_left_ground")
Einzel3_left_ground.define_2x3_mapping(DXFSolid.rotz)
geom.set_solid(14, Einzel3_left_ground)

Einzel3 = DXFSolid(dxfPIAlinear, "Einzel3")
Einzel3.define_2x3_mapping(DXFSolid.rotz)
geom.set_solid(15, Einzel3)

Einzel3_right_ground = DXFSolid(dxfPIAlinear, "Einzel3_right_ground")
Einzel3_right_ground.define_2x3_mapping(DXFSolid.rotz)
geom.set_solid(16, Einzel3_right_ground)

Mask = DXFSolid(dxfPIAlinear, "Mask")
Mask.define_2x3_mapping(DXFSolid.rotz)
geom.set_solid(17, Mask)

# // Simulation boundaries

geom.set_boundary(1, Bound(BoundaryTypes.NEUMANN, 0.0))  # xmin
geom.set_boundary(2, Bound(BoundaryTypes.NEUMANN, 0.0))  # xmax
geom.set_boundary(3, Bound(BoundaryTypes.NEUMANN, 0.0))  # ymin
geom.set_boundary(4, Bound(BoundaryTypes.NEUMANN, 0.0))  # ymax
geom.set_boundary(5, Bound(BoundaryTypes.NEUMANN, 0.0))  # zmin

# // User defined electrodes

geom.set_boundary(7, Bound(BoundaryTypes.DIRICHLET, 0.0))  # plasma
geom.set_boundary(8, Bound(BoundaryTypes.DIRICHLET, 0.0))  # gnd
geom.set_boundary(9, Bound(BoundaryTypes.DIRICHLET, -V_Einzel1))  # Einzel1
geom.set_boundary(10, Bound(BoundaryTypes.DIRICHLET, 0.0))  # gnd
geom.set_boundary(11, Bound(BoundaryTypes.DIRICHLET, 0.0))  # gnd
geom.set_boundary(12, Bound(BoundaryTypes.DIRICHLET, -V_Einzel2))  # Einzel2
geom.set_boundary(13, Bound(BoundaryTypes.DIRICHLET, 0.0))  # gnd
geom.set_boundary(14, Bound(BoundaryTypes.DIRICHLET, 0.0))  # gnd
geom.set_boundary(15, Bound(BoundaryTypes.DIRICHLET, -V_Einzel3))  # Einzel3
geom.set_boundary(16, Bound(BoundaryTypes.DIRICHLET, 0.0))  # gnd
geom.set_boundary(17, Bound(BoundaryTypes.DIRICHLET, 0.0))  # Mask

# // Build mesh and surface triangulation
geom.build_mesh()

# // Construct necessary fields
epot = EpotField(geom)

scharge = MeshScalarField(geom)
bfield = MeshVectorField()
efield = EpotEfield(epot)

efldextrpl = [FieldExtrapolE.FIELD_EXTRAPOLATE] * 6


#   field_extrpl_e efldextrpl[6] = {FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
# 				   FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
# 				   FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE};
efield.set_extrapolation(efldextrpl)

# // Construct solver and particle database

solver = EpotBiCGSTABSolver(geom)
pdb = ParticleDataBase3D(geom)

pmirror = [False, False, False, False, False, False]
pdb.set_mirror(pmirror)

# // Bfield suppression in plasma region, will be suppressed if smaller than epot threshold //

psup = PPlasmaBfieldSuppression(epot, 10.0)
pdb.set_bfield_suppression(psup)

rhoe = pdb.get_rhosum()
print(rhoe)
# // Loop for running the vlasov iteration
nIteration = 1
for i in range(nIteration):
    print("Iteration: ", i)
    solver.solve(epot, scharge)
    efield.recalculate()

    # if i == 1:
    #    rhoe = pdb.get_rhosum()
    #    solver.set_pexp_plasma(rhoe, Te, Up)

    pdb.clear()
    # //    double j_p = beamIntensity_mA / (M_PI*rPlasma * rPlasma);  // beam density
    # //    double sizeBeamMax = M_PI * rMaxBeam * rMaxBeam;  // maximum size of beam
    # //    double sizeMeshSquare = mesh * mesh;
    # //    double numberMeshMax = sizeBeamMax / sizeMeshSquare; // number of meshs in max beam
    # //    double trajectoriesBeamP = numberMeshMax *  trajectoriesPerMeshSquare;
    # //    double trajectoriesP = (radiusPlasmaCyl * radiusPlasmaCyl) / (rPlasma * rPlasma) * trajectoriesBeamP;

    # //    cout<< "Proton beam density: " << j_p << "A/m^2" <<endl;
    # //    cout<< "Proton trajectories: " << trajectoriesP << "A/m^2" <<endl;

    pdb.add_cylindrical_beam_with_energy(
        10000,
        1,
        1,
        1,
        30.0e3,
        0,
        0,
        Vec3D(0, 0, -10.0e-3),
        Vec3D(1, 0, 0),
        Vec3D(0, 1, 0),
        68.0e-3,
    )
    pdb.iterate_trajectories(scharge, efield, bfield)

# // Define particle beam


# // Launch interactive plotter
#      if (gtkplotterON)
#      {
#      GTKPlotter plotter( &argc,&argv );
#      plotter.set_geometry( &geom );
#      plotter.set_epot( &epot );
#      plotter.set_particledatabase( &pdb );
#      plotter.set_efield( &efield );
#      plotter.set_bfield( &bfield );
#      plotter.set_scharge( &scharge );
#      plotter.new_geometry_plot_window();
#      plotter.run();
#      }

# }

print("----------------------- All good --------------------------")
