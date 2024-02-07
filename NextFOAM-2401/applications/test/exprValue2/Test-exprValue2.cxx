/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

Description
    Test exprValue uniformity checks

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "labelList.H"
#include "scalarField.H"
#include "vectorField.H"
#include "DynamicList.H"
#include "Random.H"
#include "exprValue.H"
#include "exprValueFieldTag.H"

using namespace Foam;

Ostream& printInfo(const expressions::exprValueFieldTag& tag)
{
    tag.print(Pout);
    return Pout;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    using fieldTag = expressions::exprValueFieldTag;

    Random rnd(123456);

    argList::noCheckProcessorDirectories();

    #include "setRootCase.H"

    {
        scalarField fld1(20);
        scalarField fld2a(20, Zero);
        scalarField fld2b(10, 3.10);
        scalarField fld3;

        forAll(fld1, i)
        {
            fld1[i] = rnd.position<scalar>(0, 20);
        }

        fieldTag tag1(fld1.begin(), fld1.end());
        fieldTag tag2a(fld2a.begin(), fld2a.end());
        fieldTag tag2b(fld2b.begin(), fld2b.end());
        fieldTag tag3(fld3.begin(), fld3.end());

        printInfo(tag1) << nl;
        printInfo(tag2a) << nl;
        printInfo(tag2b) << nl;
        printInfo(tag3) << nl;

        fieldTag result;

        result = fieldTag::combineOp{}(tag1, tag2a);
        printInfo(result) << nl;

        Info<< "combine: ";
        printInfo(tag2a) << " and ";
        printInfo(tag2b) << nl;

        result = fieldTag::combineOp{}(tag2a, tag2b);
        printInfo(result) << nl;
    }

    {
        vectorField fld2a(20, Foam::zero{});
        vectorField fld2b(10, vector::uniform(3.10));

        fieldTag tag2a(fld2a.begin(), fld2a.end());
        fieldTag tag2b(fld2b.begin(), fld2b.end());

        printInfo(tag2a) << nl;
        printInfo(tag2b) << nl;

        fieldTag result;

        Info<< "combine: ";
        printInfo(tag2a) << " and ";
        printInfo(tag2b) << nl;

        result = fieldTag::combineOp{}(tag2a, tag2b);
        printInfo(result) << nl;
    }

    if (UPstream::parRun())
    {
        scalarField fld;

        if (!UPstream::master())
        {
            fld.resize(UPstream::myProcNo(), UPstream::myProcNo());
        }

        fieldTag oldTag(fld.begin(), fld.end());

        printInfo(oldTag) << " input: " << fld << " <- before reduce" << nl;

        fieldTag newTag = returnReduce(oldTag, fieldTag::combineOp());

        printInfo(newTag) << " <- after reduce" << nl;
    }


    Info<< "\nEnd\n" << nl;

    return 0;
}


// ************************************************************************* //
