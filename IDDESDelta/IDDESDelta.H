/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::IDDESDelta

Description
    IDDESDelta used by the IDDES (improved low Re Spalart-Allmaras DES model)
    The min and max delta are calculated using the face to face distance of
    the cell.

SourceFiles
    IDDESDelta.C

\*---------------------------------------------------------------------------*/

#ifndef LESModels_IDDESDelta_H
#define LESModels_IDDESDelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                         Class IDDESDelta Declaration
\*---------------------------------------------------------------------------*/

class IDDESDelta
:
    public LESdelta
{
    // Private data

        //- Run-time selectable delta for hmax (inner: used in Cw*hmax term)
        //  Defaults to the maxDeltaxyz model if not supplied
        autoPtr<LESdelta> hmaxPtr_;

        //- Run-time selectable delta for hmax (outer: used as upper-bound clip)
        //  Defaults to hmaxPtr_ if not supplied (i.e. same as inner hmax)
        autoPtr<LESdelta> hmaxOuterPtr_;

        scalar Cw_;

        //- Shielding function (f~d) supplied by the IDDES turbulence model,
        //  used to switch the outer upper-bound clip between hmax and
        //  hmaxOuter. Defaults to zero (i.e. behaves as hmaxOuter) until
        //  the owning turbulence model calls updateFd()
        volScalarField fd_;

        //- Threshold on fd_ above which the outer clip is restricted to
        //  hmax (i.e. Delta_max) rather than hmaxOuter
        scalar fdLimit_;

        //- Coefficient for the volumetric lower-bound floor Clim*V^(1/3)
        //  applied to hmaxOuter in the fd_ < fdLimit_ branch, preventing
        //  the outer clip from collapsing on very fine/anisotropic cells.
        //  Defaults to 0.15
        scalar Clim_;

        //- Cached face-to-face geometric distance (expensive to compute,
        //  depends only on the mesh - not on fd_)
        volScalarField faceToFaceMax_;


    // Private Member Functions

        //- No copy construct
        IDDESDelta(const IDDESDelta&) = delete;

        //- No copy assignment
        void operator=(const IDDESDelta&) = delete;

        //- Calculate the (expensive) mesh-geometric quantities.
        //  Only needs to be called on construction and when the mesh changes
        void calcGeometricDelta();

        //- Combine the cached geometric quantities with the current fd_
        //  to (re)calculate delta. Cheap - safe to call every time-step
        void calcDelta();


public:

    //- Runtime type information
    TypeName("IDDESDelta");


    // Constructors

        //- Construct from name, mesh and IOdictionary
        IDDESDelta
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );


    // Destructor

        ~IDDESDelta()
        {}


    // Member Functions

        //- Read the LESdelta dictionary
        void read(const dictionary&);

        //- Return the inner hmax delta field
        const volScalarField& hmax() const
        {
            return *hmaxPtr_;
        }

        //- Return the outer hmax delta field (upper-bound clip)
        const volScalarField& hmaxOuter() const
        {
            return hmaxOuterPtr_ ? *hmaxOuterPtr_ : *hmaxPtr_;
        }

        //- Update the shielding function used to select the outer
        //  upper-bound clip (hmax where fd >= fdLimit_, hmaxOuter
        //  elsewhere), and recalculate delta
        void updateFd(const volScalarField& fd);

        // Correct values
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
