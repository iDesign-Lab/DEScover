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

\*---------------------------------------------------------------------------*/

#include "DeltaLambda2Delta.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(DeltaLambda2Delta, 0);
    addToRunTimeSelectionTable(LESdelta, DeltaLambda2Delta, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::DeltaLambda2Delta::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();

    const volVectorField& U = turbulenceModel_.U();

    // Velocity gradient tensor
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volTensorField& gradU = tgradU.cref();

    const cellList& cells          = mesh.cells();
    const vectorField& cellCentres = mesh.cellCentres();
    const vectorField& faceCentres = mesh.faceCentres();

    forAll(cells, celli)
    {
        const tensor& gradUc = gradU[celli];

        // Symmetric and antisymmetric parts of the velocity gradient.
        // skew(gradU) = -Omega in the OpenFOAM convention, which does not
        // matter here since only Omega & Omega is needed.
        const symmTensor S = symm(gradUc);
        const tensor W = skew(gradUc);

        // M = S & S + Omega & Omega (symm() makes the symmetry exact)
        const symmTensor M = symm((S & S) + (W & W));

        // eigenValues(symmTensor) returns the eigenvalues in ascending
        // order, and the rows of eigenVectors() are the corresponding unit
        // eigenvectors, so the largest eigenvalue is eVals.z() and its
        // eigenvector is the third row eVecs.z()
        const vector eVals = eigenValues(M);
        const tensor eVecs = eigenVectors(M, eVals);

        // Unit vortex axis: eigenvector for the largest eigenvalue of M
        const vector nR = eVecs.z();

        // Delta_R = max over faces of 2*|n_R ^ (fc - cc)|
        const labelList& cFaces = cells[celli];
        const point& cc = cellCentres[celli];

        scalar deltaMaxTmp = 0;
        for (const label facei : cFaces)
        {
            const scalar tmp =
                2.0*mag(nR ^ (faceCentres[facei] - cc));
            if (tmp > deltaMaxTmp)
            {
                deltaMaxTmp = tmp;
            }
        }

        delta_[celli] = deltaCoeff_*deltaMaxTmp;
    }

    const label nD = mesh.nGeometricD();

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

    delta_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::DeltaLambda2Delta::DeltaLambda2Delta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    deltaCoeff_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "deltaCoeff",
            1.035
        )
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::DeltaLambda2Delta::read(const dictionary& dict)
{
    const dictionary& coeffsDict = dict.optionalSubDict(type() + "Coeffs");
    coeffsDict.readIfPresent<scalar>("deltaCoeff", deltaCoeff_);
    calcDelta();
}


void Foam::LESModels::DeltaLambda2Delta::correct()
{
    calcDelta();
}


// ************************************************************************* //
