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

#ifndef rbsmbeam3d_h
#define rbsmbeam3d_h

#include "sm/Elements/Beams/beambaseelement.h"
#include "sm/CrossSections/fiberedcs.h"
#include "sm/CrossSections/rbsfiberedcs.h"
#include "sm/Materials/winklermodel.h"
#include "dofmanager.h"
#include "vtkxmlexportmodule.h"

#include "libeam3d.h"

///@name Input fields for Beam3d
//@{
#define _IFT_RBSMBeam3d_Name "rbsmbeam3d"
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
 * This class is the Mindlin version of the RBSM springs beam class (rbsbeam3d) which
 * instead of extending the beam3d class, extends the libeam3d to implement RBSM springs
 *
 * @author Saeid Mehrpay
 */
class RBSMBeam3d : public LIBeam3d
{
public:
protected:
    double referenceAngle = 0;
    int referenceNode;
    FloatArray zaxis;


public:
    RBSMBeam3d( int n, Domain *d );
    virtual ~RBSMBeam3d();

    void FiberedCrossSectionInterface_computeStrainVectorInFiber(
        FloatArray &answer, const FloatArray &masterGpStrain,
        GaussPoint *slaveGp, TimeStep *tStep ) override;

    int giveLocalCoordinateSystem( FloatMatrix &answer ) override;

    // definition & identification
    const char *giveClassName() const override { return "RBSMBeam3d"; }
    const char *giveInputRecordName() const override { return _IFT_RBSMBeam3d_Name; }
    void initialize();
    void initializeFrom( InputRecord &ir ) override;

    virtual void setNumberOfGaussPoints( int nip );


protected:
    void giveInternalForcesVectorAtPoint( FloatArray &answer, TimeStep *tStep, FloatArray &coords );
    void computeInternalForcesFromBoundaryEdgeLoadVectorAtPoint(
        FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode,
        TimeStep *tStep, FloatArray &pointCoords, double ds, bool global );
    void computeInternalForcesFromBodyLoadVectorAtPoint(
        FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode,
        FloatArray &pointCoords, double ds );
};

} // end namespace oofem
#endif // rbsbeam3d_h
