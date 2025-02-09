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
    chtMultiRegionSimpleNFoam

Group
    grpHeatTransferSolvers

Description
    Steady-state solver for buoyant, turbulent fluid flow and solid heat
    conduction with conjugate heat transfer between solid and fluid regions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "multiRegionSimpleControl.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "loopControl.H"
#include "wallFvPatch.H"

#define NFOAM

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for buoyant, turbulent fluid flow and solid heat"
        " conduction with conjugate heat transfer"
        " between solid and fluid regions."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    forAll(fluidRegions, i)
    {
        turbulenceFluid[i]->validate();
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (pvCoupling.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while(pvCoupling.fluidLoop())
        {
            #include "setRegionFluidFields.H"

            if (!pvCoupling.initialState())
            {
                MRF.makeRelative(rhof, phi);
            }

            // Pressure-velocity SIMPLE corrector
            if (pvCoupling.solveFlow())
            {
                #include "UEqn.H"
                #include "pEqn.H"
            }

            #include "createSpeciesDiffusionEnthapyTransport.H"

            if (pvCoupling.solveEnergy())
            {
                #include "EEqn.H"

                if (!thermo.isochoric() || Y.size())
                {
                    #include "updateDensity.H"
                }
            }
            else
            {
                if (!thermo.incompressible() || Y.size())
                {
                    thermo.correctRhoMu();
                    #include "updateDensity.H"
                }
            }

            #include "YEqn.H"

            if (pvCoupling.solveFlow())
            {
                turbulence->correct();
            }

            #include "continuityErrs.H"

            MRF.makeAbsolute(rhof, phi);
        }

        while(pvCoupling.solidLoop())
        {
            #include "setRegionSolidFields.H"

            if (pvCoupling.solveEnergy())
            {
                #include "hEqn.H"
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
