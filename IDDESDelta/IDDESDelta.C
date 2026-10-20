/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "IDDESDelta.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"
#include "maxDeltaxyz.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(IDDESDelta, 0);
    addToRunTimeSelectionTable(LESdelta, IDDESDelta, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::IDDESDelta::calcGeometricDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();

    // Wall-normal vectors
    const volVectorField& n = wallDist::New(mesh).n();

    scalarField& faceToFacenMax = faceToFaceMax_.primitiveFieldRef();

    const cellList& cells = mesh.cells();
    const vectorField& faceCentres = mesh.faceCentres();

    forAll(cells, celli)
    {
        scalar maxDelta = 0.0;
        const labelList& cFaces = cells[celli];
        const vector nci = n[celli];

        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];
            const point& fci = faceCentres[facei];

            forAll(cFaces, cFacej)
            {
                label facej = cFaces[cFacej];
                const point& fcj = faceCentres[facej];
                scalar ndfc = nci & (fcj - fci);

                if (ndfc > maxDelta)
                {
                    maxDelta = ndfc;
                }
            }
        }

        faceToFacenMax[celli] = maxDelta;
    }

    label nD = mesh.nGeometricD();

    if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, LES is not strictly applicable" << nl
            << endl;
    }
    else if (nD != 3)
    {
        FatalErrorInFunction
            << "Case must be either 2D or 3D" << exit(FatalError);
    }
}


void Foam::LESModels::IDDESDelta::calcDelta()
{
    const volScalarField& hmax      = hmaxPtr_();
    const volScalarField& hmaxOuter = hmaxOuterPtr_ ? hmaxOuterPtr_() : hmax;
    const fvMesh& mesh = turbulenceModel_.mesh();

    // Volumetric lower-bound floor Clim*V^(1/3), applied only within the
    // fd < fdLimit_ (outer, LES-favouring) branch to prevent hmaxOuter from
    // collapsing to an unphysically small value on very fine/anisotropic
    // cells. Clim_ defaults to 0.15.
    volScalarField climCbrtV
    (
        IOobject
        (
            "IDDESDelta.climCbrtV",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimLength, Zero)
    );
    climCbrtV.primitiveFieldRef() = Clim_*cbrt(mesh.V());

    // Outer upper-bound clip: use hmax (Delta_max) where the shielding
    // function fd_ is at or above fdLimit_, otherwise fall back to
    // max(hmaxOuter, Clim*V^(1/3))
    const volScalarField hmaxOuterLim
    (
        IOobject
        (
            "hmaxOuterLim",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pos0(fd_ - fdLimit_)*hmax
      + neg(fd_ - fdLimit_)*max(hmaxOuter, climCbrtV)
    );

    delta_.primitiveFieldRef() =
        min
        (
            max
            (
                max
                (
                    Cw_*wallDist::New(mesh).y(),
                    Cw_*hmax           // inner hmax: geometric scale
                ),
                faceToFaceMax_
            ),
            hmaxOuterLim               // outer hmax: upper-bound clip
        );

    // Handle coupled boundaries
    delta_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::IDDESDelta::IDDESDelta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    hmaxPtr_(nullptr),
    hmaxOuterPtr_(nullptr),
    Cw_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "Cw",
            0.15
        )
    ),
    fd_
    (
        IOobject
        (
            "IDDESDelta.fd",
            turbulence.mesh().time().timeName(),
            turbulence.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence.mesh(),
        dimensionedScalar(dimless, Zero)
    ),
    fdLimit_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "fdLimit",
            0.01
        )
    ),
    Clim_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "Clim",
            0.15
        )
    ),
    faceToFaceMax_
    (
        IOobject
        (
            "faceToFaceMax",
            turbulence.mesh().time().timeName(),
            turbulence.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence.mesh(),
        dimensionedScalar(dimLength, Zero)
    )
{
    const dictionary& coeffsDict = dict.optionalSubDict(type() + "Coeffs");

    // Inner hmax (used in Cw*hmax term)
    if (coeffsDict.found("hmax"))
    {
        hmaxPtr_ =
            LESdelta::New
            (
                IOobject::groupName("hmax", turbulence.U().group()),
                turbulence,
                coeffsDict,
                "hmax"
            );
    }
    else
    {
        Info<< "Employing " << maxDeltaxyz::typeName << " for hmax" << endl;

        hmaxPtr_.reset
        (
            new maxDeltaxyz
            (
                IOobject::groupName("hmax", turbulence.U().group()),
                turbulence,
                coeffsDict
            )
        );
    }

    // Outer hmax (used as upper-bound clip); optional, defaults to inner hmax
    if (coeffsDict.found("hmaxOuter"))
    {
        Info<< "Employing user-defined hmaxOuter for upper-bound clip" << endl;

        hmaxOuterPtr_ =
            LESdelta::New
            (
                IOobject::groupName("hmaxOuter", turbulence.U().group()),
                turbulence,
                coeffsDict,
                "hmaxOuter"
            );
    }

    calcGeometricDelta();
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::IDDESDelta::read(const dictionary& dict)
{
    const dictionary& coeffsDict(dict.optionalSubDict(type() + "Coeffs"));

    coeffsDict.readIfPresent<scalar>("Cw", Cw_);
    coeffsDict.readIfPresent<scalar>("fdLimit", fdLimit_);
    coeffsDict.readIfPresent<scalar>("Clim", Clim_);

    calcDelta();
}


void Foam::LESModels::IDDESDelta::updateFd(const volScalarField& fd)
{
    fd_ = fd;
    calcDelta();
}


void Foam::LESModels::IDDESDelta::correct()
{
    if (turbulenceModel_.mesh().changing())
    {
        hmaxPtr_->correct();

        if (hmaxOuterPtr_)
        {
            hmaxOuterPtr_->correct();
        }

        calcGeometricDelta();
        calcDelta();
    }
}


// ************************************************************************* //
