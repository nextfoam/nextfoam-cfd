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

Class
    Foam::expressions::exprValueFieldTag

Description
    A polymorphic single-value container for tracking Field content
    as uniform etc.

SourceFiles
    exprValueFieldTag.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_expressions_exprValueFieldTag_H
#define Foam_expressions_exprValueFieldTag_H

#include "exprValue.H"
#include "List.H"  // For ListPolicy

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

/*---------------------------------------------------------------------------*\
                      Class exprValueFieldTag Declaration
\*---------------------------------------------------------------------------*/

class exprValueFieldTag
{
    // Private data

        //- Uniformity of field
        int uniformity_{};

        //- Representative (uniform) value for the field
        expressions::exprValue value_{};

public:

    // Constructors

        //- Default construct
        exprValueFieldTag() = default;

        //- Construct from a range of values
        template<class Type>
        explicit exprValueFieldTag(const Type* first, const Type* last)
        {
            set(first, last);
        }


    // Member Functions

        //- True if the uniformity is "empty"
        bool empty() const noexcept;

        //- True if the uniformity is "uniform"
        bool is_uniform() const noexcept;

        //- True if the uniformity is "non-uniform"
        bool is_nonuniform() const noexcept;

        //- Test for equality of the values
        bool equal(const exprValueFieldTag& rhs) const;

        //- Set value and uniformity from range of data
        template<class Type>
        void set(const Type* first, const Type* last)
        {
            uniformity_ =
                Foam::Detail::ListPolicy::check_uniformity(first, last);

            if (first != last)  // or uniformity_ != EMPTY
            {
                value_.set<Type>(*first);
            }
            else
            {
                value_.set<Type>(Foam::zero{});
            }
        }

        //- Set uniform type and value
        template<class Type>
        void set_uniform(const Type& val)
        {
            uniformity_ = Foam::Detail::ListPolicy::uniformity::UNIFORM;
            value_.set<Type>(val);
        }

        //- Set as non-uniform
        void set_nouniform() noexcept;


    // Reduction operations

        //- Combine - eg, for global uniformity
        void combine(const exprValueFieldTag& b);

        //- Binary operator to be used by reduce function for detecting
        //- global uniformity
        struct combineOp
        {
            exprValueFieldTag operator()
            (
                const exprValueFieldTag& a,
                const exprValueFieldTag& b
            ) const
            {
                exprValueFieldTag result(a);
                result.combine(b);
                return result;
            }
        };


    // IO Operations

        void read(Istream& is);
        void write(Ostream& os) const;

        //- Print description to Ostream
        void print(Ostream& os) const;
};

} // End namespace expressions


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// IO
Istream& operator>>(Istream&, expressions::exprValueFieldTag&);
Ostream& operator<<(Ostream&, const expressions::exprValueFieldTag&);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
