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

#include "DeltaLiutexTildeDelta.H"
#include "fvcGrad.H"
#include "maxDeltaxyz.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(DeltaLiutexTildeDelta, 0);
    addToRunTimeSelectionTable(LESdelta, DeltaLiutexTildeDelta, dictionary);
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
namespace
{

// Relative threshold on the swirling strength. Cells with
// lambda_ci <= swirlTol*mag(A) are treated as non-vortical. This only
// removes round-off-level swirl (e.g. in a uniform or one-dimensional flow,
// where the discriminant is zero analytically but ~1e-16 numerically).
constexpr scalar swirlTol = 1e-6;

// Unit Liutex (Rortex) vector n_R of the velocity-gradient tensor
// A_ij = du_i/dx_j, sign-aligned with the vorticity.
// Returns zero where no rigid rotation exists.
vector liutexAxis(const tensor& A)
{
    // Invariants of A
    const scalar I1 = tr(A);
    const scalar I2 =
        A.xx()*A.yy() + A.xx()*A.zz() + A.yy()*A.zz()
      - A.xy()*A.yx() - A.xz()*A.zx() - A.yz()*A.zy();
    const scalar I3 = det(A);

    // Characteristic equation lambda^3 - I1 lambda^2 + I2 lambda - I3 = 0.
    // With lambda = y + I1/3:  y^3 + p y + q = 0
    const scalar p = (3.0*I2 - sqr(I1))/3.0;
    const scalar q = (-2.0*pow3(I1) + 9.0*I1*I2 - 27.0*I3)/27.0;

    // One real and two complex-conjugate eigenvalues if disc > 0
    const scalar disc = sqr(q/2.0) + pow3(p/3.0);
    if (disc <= 0)
    {
        return Zero;
    }

    const scalar sqrtDisc = sqrt(disc);
    const scalar S1 = cbrt(-q/2.0 + sqrtDisc);
    const scalar S2 = cbrt(-q/2.0 - sqrtDisc);

    const scalar lambdaR = S1 + S2 + I1/3.0;
    const scalar lambdaCi = mag(0.5*sqrt(3.0)*(S1 - S2));

    if (lambdaCi <= swirlTol*max(mag(A), VSMALL))
    {
        return Zero;
    }

    // Real (right) eigenvector: null vector of B = A - lambdaR*I.
    // B has rank 2, so the null vector is the cross product of two
    // independent rows; take the best-conditioned pair.
    const tensor B(A - lambdaR*tensor::I);
    const vector b0(B.x());
    const vector b1(B.y());
    const vector b2(B.z());

    vector r(b0 ^ b1);
    scalar rMagSqr = magSqr(r);

    const vector c1(b0 ^ b2);
    if (magSqr(c1) > rMagSqr)
    {
        r = c1;
        rMagSqr = magSqr(c1);
    }

    const vector c2(b1 ^ b2);
    if (magSqr(c2) > rMagSqr)
    {
        r = c2;
        rMagSqr = magSqr(c2);
    }

    if (rMagSqr <= VSMALL)
    {
        return Zero;
    }

    r /= sqrt(rMagSqr);

    // Vorticity from the same tensor: omega_i = eps_ijk A_kj
    const vector omega
    (
        A.zy() - A.yz(),
        A.xz() - A.zx(),
        A.yx() - A.xy()
    );

    // Sign correction n_R = sign(omega . r) r
    scalar Wr = omega & r;
    if (Wr < 0)
    {
        r = -r;
        Wr = -Wr;
    }

    // Liutex magnitude R = omega.r - sqrt((omega.r)^2 - 4 lambda_ci^2)
    // (Gao & Liu 2018)
    const scalar R = Wr - sqrt(max(sqr(Wr) - 4.0*sqr(lambdaCi), scalar(0)));

    if (R <= 0)
    {
        return Zero;
    }

    return r;
}

} // End anonymous namespace
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::DeltaLiutexTildeDelta::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();
    const volVectorField& U0 = turbulenceModel_.U();

    // OpenFOAM convention: fvc::grad(U)_ij = du_j/dx_i. The Liutex axis is
    // defined from A_ij = du_i/dx_j, i.e. the transpose.
    tmp<volTensorField> tgradU = fvc::grad(U0);
    const volTensorField& gradU = tgradU.cref();

    vectorField nvecLiutex(mesh.nCells());
    forAll(nvecLiutex, celli)
    {
        nvecLiutex[celli] = liutexAxis(gradU[celli].T());
    }
    tgradU.clear();

    const cellList& cells = mesh.cells();
    const vectorField& cellCentres = mesh.cellCentres();
    const vectorField& faceCentres = mesh.faceCentres();

    forAll(cells, celli)
    {
        const labelList& cFaces = cells[celli];
        const point& cc = cellCentres[celli];
        const vector& nv = nvecLiutex[celli];

        scalar deltaMaxTmp = 0;

        for (const label facei : cFaces)
        {
            const point& fc = faceCentres[facei];
            const scalar tmp = 2.0*mag(nv ^ (fc - cc));

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

    // Handle coupled boundaries
    delta_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::DeltaLiutexTildeDelta::DeltaLiutexTildeDelta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    hmaxPtr_(nullptr),
    deltaCoeff_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "deltaCoeff",
            1.035
        )
    ),
    requireUpdate_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<bool>
        (
            "requireUpdate",
            true
        )
    )
{
    if (dict.optionalSubDict(type() + "Coeffs").found("hmax"))
    {
        // User-defined hmax
        hmaxPtr_ =
            LESdelta::New
            (
                IOobject::groupName("hmax", turbulence.U().group()),
                turbulence,
                dict.optionalSubDict("hmaxCoeffs"),
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
                dict.optionalSubDict("hmaxCoeffs")
            )
        );
    }

    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::DeltaLiutexTildeDelta::read(const dictionary& dict)
{
    const dictionary& coeffsDict = dict.optionalSubDict(type() + "Coeffs");

    coeffsDict.readIfPresent<scalar>("deltaCoeff", deltaCoeff_);
    coeffsDict.readIfPresent<bool>("requireUpdate", requireUpdate_);

    calcDelta();
}


void Foam::LESModels::DeltaLiutexTildeDelta::correct()
{
    if (turbulenceModel_.mesh().changing() && requireUpdate_)
    {
        hmaxPtr_->correct();
    }

    calcDelta();
}


// ************************************************************************* //
