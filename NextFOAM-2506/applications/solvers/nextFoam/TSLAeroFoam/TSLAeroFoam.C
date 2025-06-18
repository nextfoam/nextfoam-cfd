/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    TSLAeroFoam

Description
    Density-based compressible steady-state implicit time-marching flow solver.

Author
    JaeHeung Gill, NEXTfoam Co.,Ltd.

\*---------------------------------------------------------------------------*/

#define DBNFOAM

#include "fvCFD.H"
#include "fluidThermo.H"
#include "turbulentFluidThermoModel.H"
#include "bound.H"

#include "thermodynamicConstants.H"
#include "localTimeStep.H"
#include "godunovFlux.H"
#include "spectralRadius.H"
#include "timeMarchingControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#if defined(__clang__)
# pragma clang optimize off
#endif
void setTstd(scalar Tstd) 
{
    const_cast<scalar&>(Foam::constant::thermodynamic::Tstd) = Tstd;
}
#if defined(__clang__)
#pragma clang optimize on
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::lduMatrix::debug = 0;

    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createTimeControls.H"
    #include "createFields.H"

    timeMarchingControl lusgs(mesh, "LU-SGS", false);

    turbulence->validate();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (lusgs.loop())
    {
        #include "readTimeControls.H"
        #include "readFieldBounds.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Solve the approximate Riemann problem for this time step
        // reconstruct numerical fluxes at faces in a proper way
        Godunov->update();
        spectralRadius.update();

        localTimeStep.update
        (
            maxCo,
            spectralRadius.LambdaC()
          + scaleSpectralRadiiV*spectralRadius.LambdaV(),
            adjustLocalTimeStep
        );

        // get access to local time step sizes
        const scalarField rDeltaTau(localTimeStep.rDeltaT());

        const tmp<volScalarField> tgamma(thermo.gamma());
        const volScalarField& gamma(tgamma());

        #include "lusgsSolve.H"

        // Update primitive variables "p, U, T"
        #include "updatePrimitives.H"

        // Update RMS residuals
        #include "updateResiduals.H"

        // Update turbulence after the implicit time integration
        turbulence->correct();

        #include "updateSolverInfo.H"

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
