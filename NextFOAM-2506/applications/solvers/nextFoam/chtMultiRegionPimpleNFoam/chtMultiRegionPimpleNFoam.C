/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019,2022 OpenCFD Ltd.
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
    chtMultiRegionPimpleNFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent fluid flow and solid heat
    conduction with conjugate heat transfer between solid and fluid regions.

    It handles secondary fluid or solid circuits which can be coupled
    thermally with the main fluid region. i.e radiators, etc.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "turbulentFluidThermoModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "CorrectPhi.H"
#include "multiRegionPimpleControl.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "loopControl.H"
#include "pressureControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "bound.H"
#include "wallFvPatch.H"

#define NFOAM

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent fluid flow and solid heat"
        " conduction with conjugate heat transfer"
        " between solid and fluid regions."
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

    forAll(fluidRegions, i)
    {   
        turbulenceFluid[i]->validate();
    }

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

        while (pvCoupling.loop())
        {
            while(pvCoupling.fluidLoop())
            {
                #include "setRegionFluidFields.H"

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

                // --- Pressure-velocity PIMPLE corrector loop
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

            while(pvCoupling.solidLoop())
            {
                #include "setRegionSolidFields.H"

                if (pvCoupling.solveEnergy())
                {
                    #include "hEqn.H"
                }
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
