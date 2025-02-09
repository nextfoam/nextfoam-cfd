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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "subsonicInletFvPatchField.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Field<Type> subsonicInletFvPatchField<Type>::bcValue
(
    const fvPatchField<vector>& pU,
    const fvPatchField<scalar>& pT,
    const fluidThermo& thermo
) const
{
    const label patchi = this->patch().index();

    const tmp<volScalarField> tgamma = thermo.gamma();
    const scalarField& gamma = tgamma().boundaryField()[patchi];
    const scalarField gammai
    (
        tgamma().boundaryField()[patchi].patchInternalField()
    );

    const scalar molWeight
    (
        readScalar
        (
            thermo.subDict("mixture").subDict("specie").lookup("molWeight")
        )
    );

    const scalar R(RR/molWeight);

    const vectorField n(this->patch().nf());

    const vectorField Ui(pU.patchInternalField());
    const scalarField Ti(pT.patchInternalField());

    scalarField ci(sqrt(gammai*R*Ti));
    scalarField c0sqr(sqr(ci) + 0.5*(gammai - 1.0)*magSqr(Ui));

    scalarField Rm((Ui & n) - 2.0*ci/(gammai - 1.0));

    scalarField ctheta((flowDir_ & n)/mag(flowDir_));

    scalarField cb
    (
        (Rm*(1.0 - gamma)/(2.0 + (gamma - 1.0)*sqr(ctheta)))
        *(1.0 - ctheta*sqrt(((gamma - 1.0)*sqr(ctheta) + 2.0)*c0sqr
                               /(gamma - 1.0)/sqr(Rm) - 0.5*(gamma - 1.0)))
    );

    Field<Type> bcValue(this->patch().size(), pTraits<Type>::zero);

    if(this->internalField().name() == "p")
    {
        scalarField pb(p0_*pow(sqr(cb)/c0sqr, gamma/(gamma - 1.0)));

        bcValue = pb*pTraits<Type>::one;
    }
    else if(this->internalField().name() == "T")
    {
        scalarField Tb(T0_*sqr(cb)/c0sqr);

        bcValue = Tb*pTraits<Type>::one;
    }
    else
    {
        FatalErrorInFunction
            << "This boundary condition should be used for p, U and T."
            << exit(FatalError);
    }

    return bcValue;
}


template<>
Field<vector> subsonicInletFvPatchField<vector>::bcValue
(
    const fvPatchField<vector>& pU,
    const fvPatchField<scalar>& pT,
    const fluidThermo& thermo
) const
{
    const label patchi = this->patch().index();

    const tmp<volScalarField> tgamma = thermo.gamma();
    const scalarField& gamma = tgamma().boundaryField()[patchi];
    const scalarField gammai
    (
        tgamma().boundaryField()[patchi].patchInternalField()
    );

    const scalar molWeight
    (
        readScalar
        (
            thermo.subDict("mixture").subDict("specie").lookup("molWeight")
        )
    );

    const scalar R(RR/molWeight);

    const vectorField n(this->patch().nf());

    const vectorField Ui(pU.patchInternalField());
    const scalarField Ti(pT.patchInternalField());

    scalarField ci(sqrt(gammai*R*Ti));
    scalarField c0sqr(sqr(ci) + 0.5*(gammai - 1.0)*magSqr(Ui));

    scalarField Rm((Ui & n) - 2.0*ci/(gammai - 1.0));

    scalarField ctheta((flowDir_ & n)/mag(flowDir_));

    scalarField cb
    (
        (Rm*(1.0 - gamma)/(2.0 + (gamma - 1.0)*sqr(ctheta)))
        *(1.0 - ctheta*sqrt(((gamma - 1.0)*sqr(ctheta) + 2.0)*c0sqr
                               /(gamma - 1.0)/sqr(Rm) - 0.5*(gamma - 1.0)))
    );

    vectorField Ub
    (
        sqrt(2.0*(c0sqr - sqr(cb))/(gamma - 1.0))*flowDir_/mag(flowDir_)
    );

    Field<vector> bcValue(this->patch().size(), pTraits<vector>::zero);

    if (this->internalField().name() == "U")
    {
        bcValue = Ub;
    }
    else
    {
        FatalErrorInFunction
            << "This boundary condition should be used for p, U and T."
            << exit(FatalError);
    }

    return bcValue;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
subsonicInletFvPatchField<Type>::subsonicInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    p0_(0.0),
    T0_(0.0),
    flowDir_(Zero),
    curTimeIndex_(-1)
{}


template<class Type>
subsonicInletFvPatchField<Type>::subsonicInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    p0_(dict.get<scalar>("p0")),
    T0_(dict.get<scalar>("T0")),
    flowDir_(dict.get<vector>("flowDir")),
    curTimeIndex_(-1)
{
    flowDir_ /= mag(flowDir_);

    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate();
    }
}


template<class Type>
subsonicInletFvPatchField<Type>::subsonicInletFvPatchField
(
    const subsonicInletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    p0_(ptf.p0_),
    T0_(ptf.T0_),
    flowDir_(ptf.flowDir_),
    curTimeIndex_(-1)
{}


template<class Type>
subsonicInletFvPatchField<Type>::subsonicInletFvPatchField
(
    const subsonicInletFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    p0_(ptf.p0_),
    T0_(ptf.T0_),
    flowDir_(ptf.flowDir_),
    curTimeIndex_(-1)
{}


template<class Type>
subsonicInletFvPatchField<Type>::subsonicInletFvPatchField
(
    const subsonicInletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    p0_(ptf.p0_),
    T0_(ptf.T0_),
    flowDir_(ptf.flowDir_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void subsonicInletFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
}


template<class Type>
void subsonicInletFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);
}


template<class Type>
void subsonicInletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Field<Type>& patchField = *this;

        const fvPatch& patch = this->patch();

        static word thermoName(fluidThermo::dictName);

        patchField = this->bcValue
        (
            patch.template lookupPatchField<volVectorField, vector>("U"),
            patch.template lookupPatchField<volScalarField, scalar>("T"),
            this->db().objectRegistry::template lookupObject<fluidThermo>
            (
                thermoName
            )
        );

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void subsonicInletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    os.writeEntry("p0", p0_);
    os.writeEntry("T0", T0_);
    os.writeEntry("flowDir", flowDir_);

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
