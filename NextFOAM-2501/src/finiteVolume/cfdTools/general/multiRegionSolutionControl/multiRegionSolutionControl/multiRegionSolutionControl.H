/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::multiRegionSolutionControl

Description
    Base class for multi-region solution control classes

\*---------------------------------------------------------------------------*/

#ifndef multiRegionSolutionControl_H
#define multiRegionSolutionControl_H

#include "fvMesh.H"
#include "Pair.H"

//- Declare that multiRegionSolutionControl will be used
#define MULTIREGION

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dynamicFvMesh;

/*---------------------------------------------------------------------------*\
                 Class multiRegionSolutionControl Declaration
\*---------------------------------------------------------------------------*/

class multiRegionSolutionControl
:
    public objectRegistry
{
public:

    //- Simple convenient storage of field residuals
    struct fieldData
    {
        wordRe name;
        scalar absTol;
        scalar relTol;
        scalar initialResidual;
    };


    // Static Member Functions

        //- Extract maximum residual for the solver performance entry,
        //- provided the corresponding volume field is available on the mesh.
        //
        //  \return initial residual as first member, the final residual
        //      as the second (or last) member
        static Pair<scalar> maxResidual
        (
            const fvMesh* fvmesh,
            const entry& dataDictEntry
        );


protected:

    // Protected data

      //- Reference to the time database
        const Time& runTime_;

        //- Pointers to the multi-region mesh database
        List<const fvMesh*> regions_;

        //- List of residual data per field per mesh region
        List<List<fieldData>> residualControls_;

        //- Pointer to the current region mesh database
        fvMesh* mesh_;

        //- The dictionary name, e.g. SIMPLE, PIMPLE
        const word algorithmName_;

        //- The solution dictionary for the current looping region
        dictionary solutionDict_;

        //- The fvSolution dictionary for global settings
        IOdictionary globalSolutionDict_;


        // Solution controls

            //- Flag to indicate whether to solve for the flow
            bool solveFlow_;

            //- Flag to indicate whether to solve for the energy
            bool solveEnergy_;

            //- Flag to indicate whether to solve for the species
            bool solveSpecies_;

            //- Maximum number of fluid regions
            label nFluidIndex_;

            //- Maximum number of solid regions
            label nSolidIndex_;

            //- Maximum number of non-orthogonal correctors
            label nNonOrthCorr_;

            //- Flag to indicate to solve for momentum
            bool momentumPredictor_;

            //- Flag to indicate to solve using transonic algorithm
            bool transonic_;

            //- Flag to indicate to relax pressure using the
            //  "consistent" approach of SIMPLEC
            bool consistent_;

            //- Flag to indicate that the flow system of equations should not
            //  be evolved
            bool frozenFlow_;


            //- Only for PIMPLE/PISO

                //- Maximum number of PISO correctors
                label nCorrPISO_;

                //- Flag to indicate whether to update density in SIMPLE
                //- rather than PISO mode
                bool SIMPLErho_;

                //- Flag to indicate whether to only solve turbulence on 
                //- final iter
                bool turbOnFinalIterOnly_;

                //- Flag to indicate wheter the final solver is used only
                //- on the final pimple iter
                bool finalOnLastPimpleIterOnly_;

                //- Flag to indicate that ddtCorr should be applied; 
                //- default = yes
                bool ddtCorr_;


        // Evolution

            //- Current fluid region index
            label fluidIndex_;

            //- Current solid region index
            label solidIndex_;

            //- Current corrector loop index
            label corr_;

            //- Current non-orthogonal corrector loop index
            label corrNonOrtho_;


        // Residuals

            bool foamLog_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual void readControl() = 0;

        //- Read controls from each fluid region's fvSolution dictionary
        void readFluidControl(const fvMesh*);

        //- Read controls from each solid region's fvSolution dictionary
        void readSolidControl(const fvMesh*);

        //- Return index of field in residualControl_ if present
        virtual label applyToField
        (
            const List<fieldData>& residualControl,
            const word& fieldName,
            const bool useRegEx = true
        ) const;

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied() = 0;

        //- Store previous iteration fields
        virtual void storePrevIterFields() const;

        //- Store previous iteration field for vol<Type>Fields
        template<class Type>
        void storePrevIter() const;

        //- Initial and final residual of the specified field-name,
        //- provided that the corresponding volume field is available
        //- on the fvMesh.
        //
        //  Populate residuals with initial residual as first member and
        //  the final residual as second (last) member.
        template<class Type>
        static bool maxTypeResidual
        (
            const fvMesh* fvmesh,
            const entry& solverPerfDictEntry,
            Pair<scalar>& residuals
        );

        //- Set the firstIteration flag on the mesh data dictionary
        virtual void setFirstIterFlag
        (
            fvMesh* fvmesh,
            const bool check = true,
            const bool force = false
        );
        

private:

        //- No copy construct
        multiRegionSolutionControl(const multiRegionSolutionControl&) = delete;

        //- No copy assignment
        void operator=(const multiRegionSolutionControl&) = delete;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("multiRegionSolutionControl");


    // Constructors

        //- Construct from mesh
        multiRegionSolutionControl
        (
            const Time& runTime, 
            const word& algorithmName,
            const PtrList<fvMesh>& fluid,
            const PtrList<fvMesh>& solid,
            bool foamLog = true
        );

        //- Construct from dynamicFvMesh
        multiRegionSolutionControl
        (
            const Time& runTime, 
            const word& algorithmName,
            const PtrList<dynamicFvMesh>& fluid,
            const PtrList<fvMesh>& solid,
            bool foamLog = true
        );


    //- Destructor
    virtual ~multiRegionSolutionControl();


    // Member Functions

        // Access

            //- Return the solution dictionary for the current looping region
            virtual const dictionary dict() const;

            //- Return the solution dictionary
            inline const dictionary& dict(const fvMesh*) const;

            //- Current corrector loop index
            inline label corr() const;

            //- Current non-orthogonal corrector index
            inline label corrNonOrtho() const;

            //- Current fluid region index
            inline label fluidIndex() const;

            //- Current solid region index
            inline label solidIndex() const;


        // Solution control

            //- Return true to solve for flow
            inline bool solveFlow() const;

            //- Return true to solve for energy
            inline bool solveEnergy() const;

            //- Return true to solve for species
            inline bool solveSpecies() const;

            //- Maximum number of fluid regions
            inline label nFluidIndex() const;

            //- Maximum number of solid regions
            inline label nSolidIndex() const;

            //- Maximum number of non-orthogonal correctors
            inline label nNonOrthCorr() const;

            //- Helper function to identify final non-orthogonal iteration
            inline bool finalNonOrthogonalIter() const;

            //- Flag to indicate to solve for momentum
            inline bool momentumPredictor() const;

            //- Flag to indicate to solve using transonic algorithm
            inline bool transonic() const;

            //- Flag to indicate to relax pressure using the
            //  "consistent" approach of SIMPLEC
            inline bool consistent() const;

            //- Flag to indicate that the flow system of equations should not
            //  be evolved
            inline bool frozenFlow() const;

            //- Flag to indicate whether the simulation is initialized or not
            //  by checking for the presence of a readable "phi"
            inline bool initialState() const;


        // Evolution

            //- Main control loop
            virtual bool loop() = 0;

            //- loop over fluid regions
            bool fluidLoop();

            //- loop over solid regions
            bool solidLoop();

            //- Non-orthogonal corrector loop
            inline bool correctNonOrthogonal();


        // Residuals printing

            inline bool foamLog() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multiRegionSolutionControlTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiRegionSolutionControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
