/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2022 OpenCFD Ltd.
    Copyright (C) 2026 Jiawei Chen, Xiang I. A. Yang, Ping He
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
    Foam::LESModels::DeltaLambda2Delta

Description
    Delta formulation that accounts for the orientation of the local vortex
    axis given by the lambda-2 criterion (Jeong & Hussain 1995). It is the
    DeltaOmegaTilde length scale of Shur et al. (2015) with the vorticity
    direction replaced by the axis n_R:

    \verbatim
        delta = deltaCoeff*max_f(2*mag(n_R ^ (fc_f - cc)))
    \endverbatim

    where fc_f are the face centres and cc is the cell centre of the cell.

    n_R is the unit eigenvector of the real symmetric tensor
    M = S & S + Omega & Omega for its largest eigenvalue, where S and Omega
    are the symmetric and antisymmetric parts of the velocity gradient.
    M is the same for grad(U) and its transpose, so no transposition is
    needed. Since M is real symmetric, n_R exists in every cell and no
    fallback is required.

    Usage
    \verbatim
    hmaxOuter       DeltaLambda2Delta;

    DeltaLambda2DeltaCoeffs
    {
        deltaCoeff      1.035;      // default 1.035
    }
    \endverbatim

    References:
    \verbatim
        Shur, M. L., Spalart, P. R., Strelets, M. K., & Travin, A. K. (2015).
        An enhanced version of DES with rapid transition
        from RANS to LES in separated flows.
        Flow, turbulence and combustion, 95(4), 709-737.
        DOI:10.1007/s10494-015-9618-0

        Jeong, J., & Hussain, F. (1995).
        On the identification of a vortex.
        Journal of Fluid Mechanics, 285, 69-94.
        DOI:10.1017/S0022112095000462
    \endverbatim

SourceFiles
    DeltaLambda2Delta.C

\*---------------------------------------------------------------------------*/

#ifndef LESModels_DeltaLambda2Delta_H
#define LESModels_DeltaLambda2Delta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                    Class DeltaLambda2Delta Declaration
\*---------------------------------------------------------------------------*/

class DeltaLambda2Delta
:
    public LESdelta
{
    // Private Data

        //- Model coefficient (default 1.035, as in DeltaOmegaTildeDelta)
        scalar deltaCoeff_;


    // Private Member Functions

        //- Calculate the delta values
        void calcDelta();

        //- No copy construct
        DeltaLambda2Delta(const DeltaLambda2Delta&) = delete;

        //- No copy assignment
        void operator=(const DeltaLambda2Delta&) = delete;


public:

    //- Runtime type information
    TypeName("DeltaLambda2Delta");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        DeltaLambda2Delta
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );


    //- Destructor
    virtual ~DeltaLambda2Delta() = default;


    // Member Functions

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        //- Correct delta values
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
