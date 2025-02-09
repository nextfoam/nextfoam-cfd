/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "basicWindSpeedFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

basicWindSpeedFvPatchVectorField::
basicWindSpeedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    ppMin_((boundBox(p.patch().localPoints())).min()),
    time_(iF.time()),
    patch_(p.patch()),
    flowDir_(nullptr),
    zDir_(nullptr),
    V0_(28.0),
    Kzt_(1.0),
    Roughness_(1),
    Iw_(1.0)
{}


basicWindSpeedFvPatchVectorField::
basicWindSpeedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    ppMin_((boundBox(p.patch().localPoints())).min()),
    time_(iF.time()),
    patch_(p.patch()),
    flowDir_(Function1<vector>::New("flowDir", dict, &time_)),
    zDir_(Function1<vector>::New("zDir", dict, &time_)),
    V0_(dict.getOrDefault("V0", 28.0)),
    Kzt_(dict.getOrDefault("Kzt", 1.0)),
    Roughness_(dict.getOrDefault("Roughness", 1)),
    Iw_(dict.getOrDefault("Iw", 1.0))
{
    fvPatchVectorField::operator==(U(patch()));
}


basicWindSpeedFvPatchVectorField::
basicWindSpeedFvPatchVectorField
(
    const basicWindSpeedFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(pvf, p, iF, mapper),
    ppMin_(pvf.ppMin_),
    time_(pvf.time_),
    patch_(p.patch()),
    flowDir_(pvf.flowDir_.clone()),
    zDir_(pvf.zDir_.clone()),
    V0_(pvf.V0_),
    Kzt_(pvf.Kzt_),
    Roughness_(pvf.Roughness_),
    Iw_(pvf.Iw_)
{}


basicWindSpeedFvPatchVectorField::
basicWindSpeedFvPatchVectorField
(
    const basicWindSpeedFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF),
    ppMin_(pvf.ppMin_),
    time_(pvf.time_),
    patch_(pvf.patch_),
    flowDir_(pvf.flowDir_.clone()),
    zDir_(pvf.zDir_.clone()),
    V0_(pvf.V0_),
    Kzt_(pvf.Kzt_),
    Roughness_(pvf.Roughness_),
    Iw_(pvf.Iw_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vectorField basicWindSpeedFvPatchVectorField::U(const fvPatch& patch) const
{
    const vectorField& pCf = patch.Cf();
    const scalar t = time_.timeOutputValue();
    const scalar groundMin = zDir() & ppMin_;

    scalarField Un
    (
        patch.size(),
        0
    );

    forAll(pCf, i)
    {
        scalar zz = zDir() & pCf[i];
        Un[i] = V0_ * Kzt_ * Iw_ * Kzr(zz);
    }

    const vectorField& flowDirField(flowDir(patch));

    vectorField Uv
    (
        patch.size(),
        Foam::vector(0,0,0)
    );

    forAll(pCf, i)
    {
        Uv[i] = flowDirField[i] * Un[i];
    }

    return Uv;
}

vector basicWindSpeedFvPatchVectorField::zDir() const
{
    const scalar t = time_.timeOutputValue();
    const vector dir(zDir_->value(t));
    const scalar magDir = mag(dir);

    if (magDir < SMALL)
    {
        FatalErrorInFunction
            << "magnitude of " << zDir_->name() << " = " << magDir
            << " vector must be greater than zero"
            << abort(FatalError);
    }

    return dir/magDir;
}

vectorField basicWindSpeedFvPatchVectorField::flowDir(const fvPatch& patch) const
{
    const scalar t = time_.timeOutputValue();
    const vector dir(flowDir_->value(t));
    const scalar magDir = mag(dir);
    vectorField flowDirField
    (   
        patch.size(),
        Foam::vector(0,0,0)
    );

    if (magDir < SMALL)
    {
        Info<< "Magnitude of flowDir is Zero, face-normal flow direction mode starts" << endl;

        const vectorField pnf(patch.nf());

        forAll(pnf, i)
        {
            vector flowDirC = zDir() ^ pnf[i] ^ zDir();
            flowDirField[i] = -flowDirC / mag(flowDirC);
        }
    }
    else
    {
        forAll(flowDirField, i)
        {
            flowDirField[i] = dir/magDir;
        }
    }

    return flowDirField;
}

scalar basicWindSpeedFvPatchVectorField::Kzr(scalar& z) const
{
    switch (Roughness_)
    {
    case 1:
        if (z <= 20.0) {
            return 0.58;
        }
        else if (z > 550) {
            return 0.22 * pow(550, 0.33);
        }
        else {
            return 0.22 * pow(z, 0.33);
        }
    
    case 2:
        if (z <= 15.0) {
            return 0.81;
        }
        else if (z > 450) {
            return 0.45 * pow(450, 0.22);
        }
        else {
            return 0.45 * pow(z, 0.22);
        }
    
    case 3:
        if (z <= 10.0) {
            return 1.0;
        }
        else if (z > 350) {
            return 0.71 * pow(350, 0.15);
        }
        else {
            return 0.71 * pow(z, 0.15);
        }
    
    case 4:
        if (z <= 5.0) {
            return 1.13;
        }
        else if (z > 250) {
            return 0.98 * pow(250, 0.1);
        }
        else {
            return 0.98 * pow(z, 0.1);
        }
    
    default:
        FatalErrorInFunction
            << "Invalid surface roughness scale"
            << abort(FatalError);
    }
}

void basicWindSpeedFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fvPatchVectorField::operator==(U(patch()));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void basicWindSpeedFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void basicWindSpeedFvPatchVectorField::rmap
(
    const fvPatchVectorField& pvf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(pvf, addr);

    const basicWindSpeedFvPatchVectorField& blpvf =
        refCast<const basicWindSpeedFvPatchVectorField>(pvf);

}


void basicWindSpeedFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    
    if (flowDir_)
    {
        flowDir_->writeData(os);
    }
    if (zDir_)
    {
        zDir_->writeData(os);
    }
    os.writeEntry("V0", V0_);
    os.writeEntry("Kzt", Kzt_);
    os.writeEntry("Roughness", Roughness_);
    os.writeEntry("Iw", Iw_);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    basicWindSpeedFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
