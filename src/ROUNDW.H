/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::ROUNDWLimiter

Description
    Class with limiter function which returns the limiter for the
    ROUNDW differencing scheme based on r obtained from the LimiterFunc
    class.

ROUNDW scheme is derived from Eq (6.6) 
Deng, Xi. "A unified framework for non-linear reconstruction schemes in a compact stencil. Part 1: Beyond second order." Journal of Computational Physics 481 (2023): 112052.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    ROUNDW.C

\*---------------------------------------------------------------------------*/

#ifndef ROUNDW_H
#define ROUNDW_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ROUNDWLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class ROUNDWLimiter
:
    public LimiterFunc
{

public:

    ROUNDWLimiter(Istream&)
    {}

    scalar limiter
    (
        const scalar cdWeight,
        const scalar faceFlux,
        const typename LimiterFunc::phiType& phiP,
        const typename LimiterFunc::phiType& phiN,
        const typename LimiterFunc::gradPhiType& gradcP,
        const typename LimiterFunc::gradPhiType& gradcN,
        const vector& d
    ) const
    {
        scalar r = LimiterFunc::r
        (
            faceFlux, phiP, phiN, gradcP, gradcN, d
        );

        scalar rr = 0.0;
        
       if(r>0.0){ rr=(0.0000116709 + r*(1.00018 + r*(40.6172 + r*(261.078 + r*(1371.44 + r*(4398.13 + r*(12835.2 + r*(24009.2 + r*(39920.9 + r*(33855.5 + r*(20208. + r*(9247.62 + r*(3150.1 + r*(883.085 +  r*(171.244 +  r*(24.6189 + (1.00035 + 0.0000233419*r)*r))))))))))))))))/(sqr(sqr(13. + r*(2. + r)))*sqr(sqr(0.0769231 + r*(0.153846 + r))));}

        return rr;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
