/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

Class
    Foam::multiRegionPimpleControl

Description
    PIMPLE control class to supply convergence information/checks for
    the multi-region PIMPLE loop.

    May also be used to for PISO-based algorithms as PISO controls are a
    sub-set of PIMPLE controls.

\*---------------------------------------------------------------------------*/

#ifndef multiRegionPimpleControl_H
#define multiRegionPimpleControl_H

#include "multiRegionSolutionControl.H"

//- Declare that PIMPLE Control will be used
#define PIMPLE_CONTROL
#define CREATE_MESH createMultiRegionMeshPostProcess.H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class multiRegionPimpleControl Declaration
\*---------------------------------------------------------------------------*/

class multiRegionPimpleControl
:
    public multiRegionSolutionControl
{
    // Private member functions
 
        //- No copy construct
        multiRegionPimpleControl(const multiRegionPimpleControl&) = delete;

        //- No copy assignment
        void operator=(const multiRegionPimpleControl&) = delete;


protected:

    // Protected data

        // Solution controls

            //- Maximum number of PIMPLE correctors
            label nCorrPIMPLE_;

            //- Current PISO corrector
            label corrPISO_;

            //- Converged flag
            bool converged_;

 
    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual void readControl();

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied();

        //- Set the firstIteration flag on the mesh data dictionary
        virtual void setFirstIterFlag
        (
            fvMesh* fvmesh,
            const bool check = true,
            const bool force = false
        );


public:

    // Static Data Members

        //- Run-time type information
        TypeName("multiRegionPimpleControl");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        multiRegionPimpleControl
        (
            const Time& runTime, 
            const PtrList<fvMesh>& fluid,
            const PtrList<fvMesh>& solid,
            bool foamLog = true
        );

        //- Construct from dynamicFvMesh and the name of control sub-dictionary
        multiRegionPimpleControl
        (
            const Time& runTime, 
            const PtrList<dynamicFvMesh>& fluid,
            const PtrList<fvMesh>& solid,
            bool foamLog = true
        );


    //- Destructor
    virtual ~multiRegionPimpleControl() = default;


    // Member Functions

        // Access

            //- Maximum number of PIMPLE correctors
            inline label nCorrPIMPLE() const;

            //- Maximum number of PISO correctors
            inline label nCorrPISO() const;

            //- Current PISO corrector index
            inline label corrPISO() const;

            //- Flag to indicate whether to update density in SIMPLE
            //  rather than PISO mode
            inline bool SIMPLErho() const;


        // Solution control

            //- PIMPLE loop
            virtual bool loop();

            //- Pressure corrector loop control
            inline bool correct();

            //- Return true to store the initial residuals
            inline bool storeInitialResiduals() const;

            //- Return true for first PIMPLE (outer) iteration
            inline bool firstIter() const;

            //- Return true for final PIMPLE (outer) iteration
            inline bool finalIter() const;

            //- Return true for final inner iteration
            inline bool finalInnerIter() const;

            //- Return true to solve for turbulence
            inline bool turbCorr();

            //- Return true to apply ddtCorr
            inline bool ddtCorr() const;

            scalar compressibleCourantNo
            (
                const fvMesh& mesh,
                const Time& runTime,
                const volScalarField& rho,
                const surfaceScalarField& phi
            ) const;

            scalar solidRegionDiffusionNo
            (
                const fvMesh& mesh,
                const Time& runTime,
                const volScalarField& Cprho,
                const volScalarField& kappa
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiRegionPimpleControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
