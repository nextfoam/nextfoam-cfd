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

\*---------------------------------------------------------------------------*/

#include "exprValueFieldTag.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::expressions::exprValueFieldTag::empty() const noexcept
{
    return
    (
        uniformity_ == Foam::Detail::ListPolicy::uniformity::EMPTY
    );
}


bool Foam::expressions::exprValueFieldTag::is_uniform() const noexcept
{
    return
    (
        uniformity_ == Foam::Detail::ListPolicy::uniformity::UNIFORM
    );
}


bool Foam::expressions::exprValueFieldTag::is_nonuniform() const noexcept
{
    return
    (
        uniformity_ == Foam::Detail::ListPolicy::uniformity::NONUNIFORM
    );
}


bool Foam::expressions::exprValueFieldTag::equal
(
    const exprValueFieldTag& rhs
) const
{
    return (value_ == rhs.value_);
}


void Foam::expressions::exprValueFieldTag::set_nouniform() noexcept
{
    uniformity_ = Foam::Detail::ListPolicy::uniformity::NONUNIFORM;
    value_ = Foam::zero{};
}


void Foam::expressions::exprValueFieldTag::combine
(
    const exprValueFieldTag& b
)
{
    if (b.empty())
    {
        // no-op
        return;
    }

    exprValueFieldTag& a = *this;

    if (a.empty())
    {
        a = b;
    }
    else if (a.is_nonuniform())
    {
        // Already non-uniform/mixed
        // a.uniformity_ |= b.uniformity_;

        a.value_ = Foam::zero{};
    }
    else if (a.is_uniform() && b.is_uniform())
    {
        // Both are uniform, but are they the same value?
        if (!a.equal(b))
        {
            a.set_nouniform();
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

void Foam::expressions::exprValueFieldTag::read(Istream& is)
{
    label uniformTag;

    is >> uniformTag;
    uniformity_ = int(uniformTag);
    value_.read(is);
}


void Foam::expressions::exprValueFieldTag::write(Ostream& os) const
{
    os << label(uniformity_);
    value_.write(os, false);  // No pruning
}


void Foam::expressions::exprValueFieldTag::print(Ostream& os) const
{
    os  << "{ uniform:"
        << label(uniformity_)
        << " type:" << label(value_.typeCode())
        << " value: " << value_ << " }";
}


Foam::Istream& Foam::operator>>
(
    Istream& is,
    expressions::exprValueFieldTag& tag
)
{
    tag.read(is);
    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const expressions::exprValueFieldTag& tag
)
{
    tag.write(os);
    return os;
}


// ************************************************************************* //
