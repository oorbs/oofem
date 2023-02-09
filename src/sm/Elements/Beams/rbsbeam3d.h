/*
 *
 *                      ____ ____  ____  ____ _____
 *                     / __   __ \/ __ \/ __ ) ___/
 *                    / / /  / / / /_/ / __  \__ \
 *                   / /_/  /_/ / _, _/ /_/ /__/ /
 *                   \____ ____/_/ |_/_____/____/
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef rbsbeam3d_h
#define rbsbeam3d_h

#include "sm/Elements/Beams/beambaseelement.h"
#include "sm/CrossSections/fiberedcs.h"
#include "sm/CrossSections/rbsfiberedcs.h"
#include "sm/Materials/winklermodel.h"
#include "dofmanager.h"
#include "vtkxmlexportmodule.h"

#include "beam3d.h"

///@name Input fields for Beam3d
//@{
#define _IFT_RBSBeam3d_Name "rbsbeam3d"
//#define _IFT_Beam3d_dofstocondense "dofstocondense"
//#define _IFT_Beam3d_refnode "refnode"
//#define _IFT_Beam3d_refangle "refangle"
//#define _IFT_Beam3d_zaxis "zaxis"
//#define _IFT_Beam3d_subsoilmat "subsoilmat"
//@}

//#define Beam3d_nSubBeams 10

namespace oofem {

//class FEI3dLineLin;

/**
 * This class extends beam3d to implement RBSM springs formulation in a beam element
 * Beam 3D cannot support nonlinearity since it does not integrate over the element;
 * but RBSM springs don't need integration as they theoretically have zero length.
 * (single integration point)
 * 
 * @author Saeid Mehrpay
 */
 class RBSBeam3d : public Beam3d
{
public:
     //void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
     //void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep) override;

 protected:

public:
    RBSBeam3d(int n, Domain *d);
    virtual ~RBSBeam3d();

    void FiberedCrossSectionInterface_computeStrainVectorInFiber(
        FloatArray &answer, const FloatArray &masterGpStrain,
        GaussPoint *slaveGp, TimeStep *tStep ) override;

    // definition & identification
    const char *giveClassName() const override { return "RBSBeam3d"; }
    const char *giveInputRecordName() const override { return _IFT_RBSBeam3d_Name; }
    void initializeFrom(InputRecord &ir) override;

    ///@todo Introduce interpolator and remove these two:
    integrationDomain giveIntegrationDomain() const override { return _Line; }
    Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }
    //Element_Geometry_Type giveGeometryType() const override { return EGT_Composite; }

    virtual void setNumberOfGaussPoints( int nip ) override;

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext & gc, TimeStep * tStep, UnknownType) override;
#endif

protected:
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &) override;

    double giveKappayCoeff(TimeStep *tStep);
    double giveKappazCoeff(TimeStep *tStep);
    void computeKappaCoeffs(TimeStep *tStep);

    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;

    MaterialMode giveMaterialMode() override { return _3dBeam; }
    int giveNumberOfIPForMassMtrxIntegration() override { return 4; }

    bool hasDofs2Condense() { return ( ghostNodes [ 0 ] || ghostNodes [ 1 ] ); }


    void giveInternalForcesVectorAtPoint(FloatArray &answer, TimeStep *tStep, FloatArray &coords);
    void computeInternalForcesFromBoundaryEdgeLoadVectorAtPoint(
        FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode,
        TimeStep *tStep, FloatArray &pointCoords, double ds, bool global);
    void computeInternalForcesFromBodyLoadVectorAtPoint(
        FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode,
        FloatArray &pointCoords, double ds);
};

} // end namespace oofem
#endif // rbsbeam3d_h
