/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "smoothSolver.H"
#include "profiling.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(smoothSolver, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<smoothSolver>
        addsmoothSolverSymMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<smoothSolver>
        addsmoothSolverAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smoothSolver::smoothSolver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{
    readControls();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::smoothSolver::readControls()
{
    lduMatrix::solver::readControls();
    nSweeps_ = controlDict_.getOrDefault<label>("nSweeps", 1);
}


Foam::solverPerformance Foam::smoothSolver::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    solveScalarField& psi = tpsi.ref();

    // Setup class containing solver performance data
    solverPerformance solverPerf(typeName, fieldName_);

    // If the nSweeps_ is negative do a fixed number of sweeps
    if (nSweeps_ < 0)
    {
        addProfiling(solve, "lduMatrix::smoother.", fieldName_);

        autoPtr<lduMatrix::smoother> smootherPtr = lduMatrix::smoother::New
        (
            fieldName_,
            matrix_,
            interfaceBouCoeffs_,
            interfaceIntCoeffs_,
            interfaces_,
            controlDict_
        );

        smootherPtr->smooth
        (
            psi,
            source,
            cmpt,
            -nSweeps_
        );

        solverPerf.nIterations() -= nSweeps_;
    }
    else
    {
        solveScalar normFactor = 0;
        solveScalarField residual;

        ConstPrecisionAdaptor<solveScalar, scalar> tsource(source);

        {
            solveScalarField Apsi(psi.size());
            solveScalarField temp(psi.size());

            // Calculate A.psi
            matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);

            // Calculate normalisation factor
            /*const lduMesh& mesh = matrix().mesh();

            if (mesh.hasDb() && mesh.thisDb().foundObject<scalarField>("UAp"))
            {
                if 
                (
                    mesh.thisDb().foundObject<scalarField>("Umag")
                 && (
                        solverPerf.fieldName() == "Ux"
                     || solverPerf.fieldName() == "Uy"
                     || solverPerf.fieldName() == "Uz"
                    )
                )
                {   
                    matrix_.Amul
                    (
                        Apsi,
                        mesh.thisDb().lookupObject<scalarField>("Umag"),
                        interfaceBouCoeffs_,
                        interfaces_,
                        cmpt
                    );
                }

                normFactor = max
                (
                    max
                    (
                        gSumMag(Apsi, mesh.comm()),
                        this->normFactor(psi, tsource(), Apsi, temp)
                    ),
                    solverPerformance::small_
                );
            } // by Gill
            else
            {
                normFactor = this->normFactor(psi, tsource(), Apsi, temp);
            }*/
            normFactor = this->normFactor(psi, tsource(), Apsi, temp);

            residual = tsource() - Apsi;

            matrix().setResidualField
            (
                ConstPrecisionAdaptor<scalar, solveScalar>(residual)(),
                fieldName_,
                true
            );

            // Calculate residual magnitude
            solverPerf.initialResidual() =
                gSumMag(residual, matrix().mesh().comm())/normFactor;
            solverPerf.finalResidual() = solverPerf.initialResidual();
        }

        if ((log_ >= 2) || (lduMatrix::debug >= 2))
        {
            Info.masterStream(matrix().mesh().comm())
                << "   Normalisation factor = " << normFactor << endl;
        }


        // Check convergence, solve if not converged
        if
        (
            minIter_ > 0
         || !solverPerf.checkConvergence(tolerance_, relTol_, log_)
        )
        {
            addProfiling(solve, "lduMatrix::smoother.", fieldName_);

            autoPtr<lduMatrix::smoother> smootherPtr = lduMatrix::smoother::New
            (
                fieldName_,
                matrix_,
                interfaceBouCoeffs_,
                interfaceIntCoeffs_,
                interfaces_,
                controlDict_
            );

            // Smoothing loop
            do
            {
                smootherPtr->smooth
                (
                    psi,
                    source,
                    cmpt,
                    nSweeps_
                );

                residual =
                    matrix_.residual
                    (
                        psi,
                        source,
                        interfaceBouCoeffs_,
                        interfaces_,
                        cmpt
                    );

                // Calculate the residual to check convergence
                solverPerf.finalResidual() =
                    gSumMag(residual, matrix().mesh().comm())/normFactor;
            } while
            (
                (
                    (solverPerf.nIterations() += nSweeps_) < maxIter_
                && !solverPerf.checkConvergence(tolerance_, relTol_, log_)
                )
             || solverPerf.nIterations() < minIter_
            );
        }

        matrix().setResidualField
        (
            ConstPrecisionAdaptor<scalar, solveScalar>(residual)(),
            fieldName_,
            false
        );
    }

    return solverPerf;
}


// ************************************************************************* //
