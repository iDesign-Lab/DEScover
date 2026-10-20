# DEScover: vortex-axis-based subgrid length scales for DES

OpenFOAM implementations of the subgrid length scales discovered by **DEScover** (automated **Dis**covery of **DES** closures), an agentic-AI framework in which large-language-model agents propose, review, implement and evaluate turbulence-model closures.

The modelling idea prescribed to the framework is that the subgrid length scale of detached eddy simulation (DES) is the grid width normal to the local vortex axis:

$$
\Delta_{\perp} = C_\delta \max_f \left( 2\,\lvert \mathbf{n}_v \times \mathbf{r}_f \rvert \right), \qquad C_\delta = 1.035,
$$

where $\mathbf{n}_v$ is the unit vector along the local vortex axis and $\mathbf{r}_f$ is the vector from the cell centre to the centre of face $f$. The candidate length scales differ only in how $\mathbf{n}_v$ is obtained from the velocity-gradient tensor. The Liutex-based length scale $\Delta_R$ gave the best results: it alleviates the delayed Kelvin–Helmholtz instability in separating shear layers and preserves the RANS shielding of attached boundary layers, without the empirical sensor of the shear-layer-adapted (SLA) length scale.

> **Paper:** J. Chen, X. I. A. Yang and P. He, *Automated discovery of physically interpretable subgrid length scales for detached eddy simulation via agentic AI*, arXiv:XXXX.XXXXX (2026). Please cite it if you use this code (see [Citation](#citation)).

## Contents

| Directory | OpenFOAM type name | Vortex axis $\mathbf{n}_v$ | Length scale in the paper |
|---|---|---|---|
| `DeltaLiutexTildeDelta/` | `DeltaLiutexTilde` | Liutex (Rortex) vector, sign-aligned with vorticity | $\Delta_R$ (**recommended**) |
| `DeltaSwirlDelta/` | `DeltaSwirlDelta` | real eigenvector of $\nabla\mathbf{u}$ (swirling strength, $\lambda_{ci}$) | $\Delta_{\text{swirl}}$ |
| `DeltaLambda2Delta/` | `DeltaLambda2Delta` | eigenvector of $\mathbf{S}^2+\boldsymbol{\Omega}^2$ for its largest eigenvalue ($\lambda_2$ criterion) | $\Delta_{\lambda_2}$ |
| `IDDESDelta/` | `IDDESDelta` | modified IDDES filter width (see below) | — |

The vorticity-based length scale $\Delta_\omega$ and the SLA length scale $\Delta_{\text{SLA}}$ used for comparison in the paper are available in OpenFOAM v2312 as `DeltaOmegaTilde` and `SLADelta`.

### Fallback in non-vortical cells

| Length scale | Behaviour where no vortex axis exists |
|---|---|
| `DeltaLiutexTilde` | $\mathbf{n}_R=\mathbf{0}$, so $\Delta_R=0$; the outer length scale is then given by the lower bound $C_{\text{lim}}V^{1/3}$ of `IDDESDelta` |
| `DeltaSwirlDelta` | falls back to $C_\delta\max_f(2\lvert\mathbf{r}_f\rvert)$, i.e. $C_\delta\Delta_{\max}$ on hexahedral cells |
| `DeltaLambda2Delta` | no fallback needed ($\mathbf{S}^2+\boldsymbol{\Omega}^2$ is real symmetric) |

### Modified `IDDESDelta`

The IDDES filter width is

$$
\Delta = \min\left[\max\left(C_w d_w,\ C_w \Delta_{\max},\ \Delta_{f2f}\right),\ \Delta_{\text{outer},lim}\right],
\qquad
\Delta_{\text{outer},lim} =
\begin{cases}
\max\left(\Delta_{\text{outer}},\ C_{\text{lim}} V^{1/3}\right), & f_d < \varepsilon,\\
\Delta_{\max}, & f_d \ge \varepsilon,
\end{cases}
$$

where $\Delta_{\text{outer}}$ is selected at run time through `hmaxOuter`. Compared with the standard OpenFOAM `IDDESDelta`, this version

- adds the run-time selectable outer length scale `hmaxOuter` (defaults to `hmax`, i.e. the standard IDDES);
- restricts it to $\Delta_{\max}$ in RANS-shielded regions ($f_d \ge$ `fdLimit`) and bounds it from below by `Clim`$\,V^{1/3}$ elsewhere;
- receives the DDES shielding function $f_d$ through `updateFd(fd)`.

## Requirements

- OpenFOAM **v2312** (www.openfoam.com). Other ESI versions with the `LESdelta` framework may work but have not been tested.

## Installation

<!-- TODO: please check this section against the way the code was built for the paper. -->

1. Copy the four directories into the LES-delta source directory of OpenFOAM, replacing the original `IDDESDelta`:
   ```bash
   cp -r DeltaLiutexTildeDelta DeltaSwirlDelta DeltaLambda2Delta IDDESDelta \
         $FOAM_SRC/TurbulenceModels/turbulenceModels/LES/LESdeltas/
   ```
2. Add the new `.C` files to the `Make/files` of the turbulence-model library, next to the existing `LESdeltas` entries, e.g.
   ```
   $(LESdelta)/DeltaLiutexTildeDelta/DeltaLiutexTildeDelta.C
   $(LESdelta)/DeltaSwirlDelta/DeltaSwirlDelta.C
   $(LESdelta)/DeltaLambda2Delta/DeltaLambda2Delta.C
   ```
3. Recompile the turbulence models:
   ```bash
   cd $FOAM_SRC/TurbulenceModels && ./Allwmake
   ```

## Usage

Select the modified `IDDESDelta` in `constant/turbulenceProperties` and choose the outer length scale with `hmaxOuter`. Example for $\Delta_R$, as used in the paper:

```cpp
simulationType  LES;

LES
{
    LESModel        SpalartAllmarasIDDES;
    turbulence      yes;
    printCoeffs     no;

    delta           IDDESDelta;

    IDDESDeltaCoeffs
    {
        Cw              0.15;               // default 0.15
        fdLimit         0.01;               // epsilon, default 0.01
        Clim            0.15;               // C_lim,   default 0.15

        // hmax         maxDeltaxyz;        // inner length scale in Cw*hmax (default)
        hmaxOuter       DeltaLiutexTilde;   // outer length scale Delta_outer

        DeltaLiutexTildeCoeffs
        {
            deltaCoeff  1.035;              // C_delta, default 1.035
        }
    }
}
```

Use `hmaxOuter DeltaSwirlDelta;` or `hmaxOuter DeltaLambda2Delta;` (with the corresponding `...Coeffs` sub-dictionary) to run the other candidates. Without `hmaxOuter`, the standard IDDES length scale ($\Delta_{\max}$) is recovered.

## Test cases

The length scales were evaluated on four cases of increasing complexity:

1. zero-pressure-gradient flat-plate boundary layer, $Re_L = 10^7$;
2. backward-facing step of Vogel & Eaton (1985), $Re_H = 28\,000$;
3. free shear layer of Delville (1995), $U_1 = 41.54$ m/s, $U_2 = 22.40$ m/s;
4. GLC305 airfoil with the 944 glaze-ice shape, $Re_c = 3.5\times10^6$, $\alpha = 6^\circ$.

## Citation

If you use this code, please cite:

```bibtex
@article{chen2026descover,
  title   = {Automated Discovery of Physically Interpretable Subgrid Length Scales for Detached Eddy Simulation via Agentic {AI}},
  author  = {Chen, Jiawei and Yang, Xiang I. A. and He, Ping},
  journal = {arXiv preprint arXiv:XXXX.XXXXX},
  year    = {2026}
}
```

## License

This code is derived from OpenFOAM and is distributed under the GNU General Public License v3 or later, in line with the licence headers of the source files.

## Acknowledgements

This work was supported by the National Science Foundation under grant numbers 2415347 and 2541867.
