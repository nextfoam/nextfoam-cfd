/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::nearWallDist

Description
    Distance calculation for cells with face on a wall.
    Searches pointNeighbours to find closest.

SourceFiles
    nearWallDist.C

\*---------------------------------------------------------------------------*/

#ifndef nearWallDist_H
#define nearWallDist_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                        Class nearWallDist Declaration
\*---------------------------------------------------------------------------*/

class nearWallDist
:
    public volScalarField::Boundary
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Lower limit for wall distance
        dimensionedScalar wallDistMin_; // by Gill


    // Private Member Functions

        //- Do all calculations
        void calculate();

        //- No copy construct
        nearWallDist(const nearWallDist&) = delete;

        //- No copy assignment
        void operator=(const nearWallDist&) = delete;


public:

    // Constructors

        //- Construct from components
        nearWallDist(const fvMesh& mesh);


    //- Destructor
    virtual ~nearWallDist();


    // Member Functions

        const volScalarField::Boundary& y() const
        {
            return *this;
        }

        //- Correct for mesh geom/topo changes
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
