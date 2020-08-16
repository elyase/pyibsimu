#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <limits>
#include <epot_bicgstabsolver.hpp>
#include <multimeshvectorfield.hpp>
#include <dxf_solid.hpp>
#include <stl_solid.hpp>
#include <stlfile.hpp>
#include <mydxffile.hpp>
#include <gtkplotter.hpp>
#include <geomplotter.hpp>
#include <geometry.hpp>
#include <func_solid.hpp>
#include <epot_efield.hpp>
#include <error.hpp>
#include <ibsimu.hpp>
#include <trajectorydiagnostics.hpp>
#include <particledatabase.hpp>
#include <particlediagplotter.hpp>
#include <fielddiagplotter.hpp>
#include <compmath.hpp>

namespace py = pybind11;


void set_message() {
    ibsimu.set_message_threshold(MSG_VERBOSE, 2);
} 

void simu(EpotBiCGSTABSolver solver, EpotField epot, MeshScalarField scharge, EpotEfield efield,MeshVectorField  bfield, ParticleDataBase3D pdb)
{
    using namespace std;
    size_t nIteration = 1;  
    for( size_t i = 0; i < nIteration; i++ ) {
        cout << "Iteration: " << i << endl;
        solver.solve( epot, scharge );
        efield.recalculate();
        pdb.clear();
        pdb.add_cylindrical_beam_with_energy( 10000, 1, 1, 1,   
                    30.0e3, 0, 0,
                    Vec3D(0,0,-10.0e-3),
                    Vec3D(1,0,0),
                    Vec3D(0,1,0),
                    68.e-3);

        pdb.iterate_trajectories( scharge, efield, bfield );
    }
}              

PYBIND11_MODULE(_ibsimu, m)
{

    m.def("set_message", &set_message, "set_message");
    //py::def("set_message_threshold", &ibsimu::set_message_threshold);

    m.def("simu", &simu, "Run script");

    py::class_<Solid>(m, "Solid");
    py::class_<MyDXFEntities>(m, "MyDXFEntities")
        .def("selection_all", &MyDXFEntities::selection_all)
        .def("scale", &MyDXFEntities::scale);

    py::class_<MyDXFEntitySelection>(m, "MyDXFEntitySelection");

    py::class_<MyDXFFile>(m, "MyDXFFile")
        .def(py::init<const std::string>())
        .def("set_warning_level", &MyDXFFile::set_warning_level)
        .def(
            "get_entities", [](MyDXFFile &self) {
                return self.get_entities();
            },
            py::return_value_policy::reference_internal);

    py::class_<DXFSolid, Solid>(m, "DXFSolid")
        .def(py::init<MyDXFFile *, const std::string &>())
        .def("rotz", &DXFSolid::rotz)
        .def("define_2x3_mapping", &DXFSolid::define_2x3_mapping)
        .def("cylindric",
             [](DXFSolid &self) {
                 self.define_2x3_mapping(DXFSolid::rotz);
             });

    py::class_<Vec3D>(m, "Vec3D")
        .def(py::init<double, double, double>())
        .def("norm2", &Vec3D::norm2);

    py::class_<Int3D>(m, "Int3D")
        .def(py::init<int32_t, int32_t, int32_t>());

    py::enum_<geom_mode_e>(m, "GeometryMode")
        .value("MODE_1D", MODE_1D)
        .value("MODE_2D", MODE_2D)
        .value("MODE_CYL", MODE_CYL)
        .value("MODE_3D", MODE_3D)
        .export_values();

    py::enum_<bound_e>(m, "BoundaryTypes")
        .value("DIRICHLET", BOUND_DIRICHLET)
        .value("NEUMANN", BOUND_NEUMANN)
        .export_values();

    py::class_<Bound>(m, "Bound")
        .def(py::init<bound_e, double>());

    py::class_<Geometry>(m, "Geometry")
        .def(py::init<geom_mode_e, Int3D, Vec3D, double>())
        .def("set_solid", &Geometry::set_solid)
        .def("build_mesh", &Geometry::build_mesh)
        .def("set_boundary", &Geometry::set_boundary);

    py::class_<ScalarField>(m, "ScalarField");

    py::class_<MeshScalarField, ScalarField>(m, "MeshScalarField")
        .def(py::init<const Geometry &>());

    py::class_<EpotField, MeshScalarField>(m, "EpotField")
        .def(py::init<const Geometry &>());

    py::class_<VectorField>(m, "VectorField");

    py::class_<MeshVectorField, VectorField>(m, "MeshVectorField")
        .def(py::init<>());

    py::class_<EpotEfield, VectorField>(m, "EpotEfield")
        .def(py::init<const EpotField &>())
        .def("set_extrapolation", [](EpotEfield &self, std::vector<field_extrpl_e> v) {
            // field_extrpl_e extrpl[v.size()];
            // std::copy(v.begin(), v.end(), extrpl);
            return self.set_extrapolation(v.data());
        })
        .def("recalculate", &EpotEfield::recalculate);

    py::enum_<field_extrpl_e>(m, "FieldExtrapolE")
        .value("FIELD_EXTRAPOLATE", FIELD_EXTRAPOLATE)
        .value("FIELD_MIRROR", FIELD_MIRROR)
        .value("FIELD_ANTIMIRROR", FIELD_ANTIMIRROR)
        .value("FIELD_SYMMETRIC_POTENTIAL", FIELD_SYMMETRIC_POTENTIAL)
        .value("FIELD_ZERO", FIELD_ZERO)
        .value("FIELD_NAN", FIELD_NAN)
        .export_values();

    py::class_<EpotMatrixSolver>(m, "EpotMatrixSolver");

    py::class_<EpotBiCGSTABSolver, EpotMatrixSolver>(m, "EpotBiCGSTABSolver")
        .def(py::init<Geometry &>())
        .def("set_pexp_plasma", &EpotBiCGSTABSolver::set_pexp_plasma)
        .def("solve", &EpotBiCGSTABSolver::solve);

    py::class_<CallbackFunctor>(m, "CallbackFunctor");
    py::class_<CallbackFunctorD_V, CallbackFunctor>(m, "CallbackFunctorD_V");

    py::class_<PPlasmaBfieldSuppression, CallbackFunctorD_V>(m, "PPlasmaBfieldSuppression")
        .def(py::init<const MeshScalarField &, double>());

    py::class_<ParticleDataBase3D>(m, "ParticleDataBase3D")
        .def(py::init([](const Geometry &geom) {
            ibsimu.set_message_threshold(MSG_VERBOSE, 2);
            ibsimu.set_thread_count(1);
            return new ParticleDataBase3D(geom);
        }))
        .def("set_mirror", [](ParticleDataBase3D &self, std::vector<bool> v) {
            // auto mirror = std::make_unique<bool[]>(v.size());
            bool mirror[v.size()];
            std::copy(std::begin(v), std::end(v), mirror);
            self.set_mirror(mirror);
        })
        .def("set_bfield_suppression", &ParticleDataBase3D::set_bfield_suppression)        
        .def("get_rhosum", &ParticleDataBase3D::get_rhosum)
        .def("clear", &ParticleDataBase3D::clear)
        .def("iterate_trajectories", &ParticleDataBase3D::iterate_trajectories)
        .def("add_cylindrical_beam_with_energy", &ParticleDataBase3D::add_cylindrical_beam_with_energy);
}