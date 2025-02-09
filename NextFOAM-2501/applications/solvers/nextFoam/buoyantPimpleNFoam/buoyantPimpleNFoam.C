/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    buoyantPimpleNFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent flow of compressible fluids
    for ventilation and heat-transfer, with optional mesh motion
    and mesh topology changes.

    Turbulence is modelled using a run-time selectable compressible RAS or
    LES model.

\*---------------------------------------------------------------------------*/

#define NFOAM

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "turbulentFluidThermoModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "radiationModel.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "bound.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent fluid flow"
        " of compressible fluids, including radiation,"
        " with optional mesh motion and mesh topology changes."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createRDeltaT.H"
    #include "createFields.H"
    #include "createRhoUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "readDyMControls.H"
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        // Store divrhoU from the previous mesh
        // so that it can be mapped and used in correctPhi
        // to ensure the corrected phi has the same divergence
        #include "createDivRhoU.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pvCoupling.loop())
        {
            if (pvCoupling.firstIter() || moveMeshOuterCorrectors)
            {
                // Store momentum to set rhoUf for introduced faces.
                if (rhoUf.valid())
                {
                    rhoU.reset(new volVectorField("rhoU", rho*U));
                }

                // Do any mesh changes
                mesh.update();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & rhoUf();

                        #include "correctPhi.H"
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            if (!pvCoupling.initialState())
            {
                MRF.makeRelative(rhof, phi);
            }

            if (mesh.changing())
            {
                // Make the fluxes relative to the mesh-motion
                fvc::makeRelative(phi, rho, U);
            }

            if (pvCoupling.solveFlow())
            {
                #include "UEqn.H"

                // --- Pressure corrector loop
                while (pvCoupling.correct())
                {
                    #include "pEqn.H"
                }
            }

            #include "createSpeciesDiffusionEnthapyTransport.H"

            if (pvCoupling.solveEnergy())
            {
                #include "EEqn.H"

                if (!thermo.isochoric() || Y.size())
                {
                    #include "updateDensity.H"
                }
                else
                {
                    #include "continuityErrs.H"
                }
            }
            else
            {
                if (!thermo.incompressible() || Y.size())
                {
                    thermo.correctRhoMu();
                    #include "updateDensity.H"
                }
                else
                {
                    #include "continuityErrs.H"
                }
            }

            #include "YEqn.H"

            if (pvCoupling.solveFlow() && pvCoupling.turbCorr())
            {
                turbulence->correct();
            }

            if (pvCoupling.finalIter() || moveMeshOuterCorrectors)
            {
                if (mesh.changing() && correctPhi)
                {
                    //Correct rhoUf if the mesh is moving
                    fvc::correctRhoUf(rhoUf, rho, U, phi);
                }
            }

            MRF.makeAbsolute(rhof, phi);
            fvc::makeAbsolute(phi, rho, U);
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
