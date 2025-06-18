/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::timeMarchingControl

Description
    Pseudo-Time control class to supply convergence information/checks for
    the Pseudo-Time loop.

\*---------------------------------------------------------------------------*/

#ifndef timeMarchingControl_H
#define timeMarchingControl_H

#include "solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class timeMarchingControl Declaration
\*---------------------------------------------------------------------------*/

class timeMarchingControl
:
    public solutionControl
{
    // Private member functions

        //- No copy construct
        timeMarchingControl(const timeMarchingControl&) = delete;

        //- No copy assignment
        void operator=(const timeMarchingControl&) = delete;


protected:

    // Protected data

        // Solution controls

            //- Initialised flag for steady solver
            bool initialised_;

            //- Maximum number of Pseudo-Time correctors
            label nCorrPseudoTime_;

            //- Converged flag
            bool converged_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual bool read();

        //- Read controls from fvSolution dictionary
        virtual bool readPseudoTime();

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied();

        //- Return true if all convergence checks are satisfied
        //- for pseudo-time loop
        virtual bool pseudoTimeCriteriaSatisfied();

        //- Set the firstIteration flag on the mesh data dictionary
        virtual void setFirstIterFlag
        (
            const bool check = true,
            const bool force = false
        );

        //- Store previous iteration fields
        virtual void storePrevIterFields() const;

        //- Store previous iteration field for vol<Type>Fields
        template<class Type>
        void storePrevIter() const;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("timeMarchingControl");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        timeMarchingControl
        (
            fvMesh& mesh,
            const word& dictName,
            const bool pseudoTime = false,
            const bool verbose = true
        );


    //- Destructor
    virtual ~timeMarchingControl() = default;


    // Member Functions

        // Access

            //- Maximum number of Pseudo-Time correctors
            inline label nCorrPseudoTime() const;


        // Solution control

            //- Time loop
            virtual bool loop();

            //- Pseudo-Time loop
            virtual bool pseudoTimeMarching();

            //- Return true to store the initial residuals
            inline bool storeInitialResiduals() const;

            //- Return true for first Pseudo-Time iteration
            inline bool firstIter() const;

            //- Return true for final Pseudo-Time iteration
            inline bool finalIter() const;

            //- Return true for final Pseudo-Time iteration
            inline bool finalInnerIter() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "timeMarchingControlTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeMarchingControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
