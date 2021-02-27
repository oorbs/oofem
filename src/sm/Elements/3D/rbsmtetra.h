/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
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

#ifndef rbsmtetra_h
#define rbsmtetra_h

#include "sm/Elements/structural3delement.h"
#include "sm/ErrorEstimators/directerrorindicatorrc.h"
#include "sm/ErrorEstimators/huertaerrorestimator.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "spatiallocalizer.h"
#include "sm/ErrorEstimators/zzerrorestimator.h"
#include "mmashapefunctprojection.h"

#include "oofemtxtinputrecord.h"

#define _IFT_RBSMTetra_Name "rbsmtetra"

namespace oofem {
class FEI3dTetLin;

/**
 * This class implements a tetrahedral rigid body - spring model.
 * Each element has one master node with 6 degree of freedom.
 * @author Saeid Mehrpay
 */
class RBSMTetra : public Structural3DElement
                /*,
                  public ZZNodalRecoveryModelInterface,
                  public NodalAveragingRecoveryModelInterface,
                  public SPRNodalRecoveryModelInterface,
                  public SpatialLocalizerInterface,
                  public ZZErrorEstimatorInterface,
                  public HuertaErrorEstimatorInterface
                */
{
protected:
    static FEI3dTetLin interpolation;

    static std::vector< IntArray > cornerCoords;
    /// number of vertices
    int numberOfCornerNodes;
    /// global ID of rigid body cell central node
    int centerDofmanager;
    /// mesh nodes defining corners of the rigid body
    IntArray cornerNodes;
    /// facets indices of rigid body
    IntArray facetArray;
    /*
    // cloned corner nodes for debug purpose
    IntArray clonedNodes;
    */


public:
    RBSMTetra(int n, Domain * d);
    virtual ~RBSMTetra() { }

    FEInterpolation *giveInterpolation() const override;

    void initializeFrom(InputRecord &ir) override;
    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    int giveNumberOfIPForMassMtrxIntegration() override { return 4; }
    Interface *giveInterface(InterfaceType it) override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_RBSMTetra_Name; }
    const char *giveClassName() const override { return "RBSMTetra"; }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawSpecial(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

    double giveRelativeSelfComputationalCost() override { return 2.15; }

/*
    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep) override;
*/

/*
    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override;
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override;
*/

    // HuertaErrorEstimatorInterface
/*
    void HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                          IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                          HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                          int &localNodeId, int &localElemId, int &localBcId,
                                                          IntArray &controlNode, IntArray &controlDof,
                                                          HuertaErrorEstimator :: AnalysisMode aMode) override;
    void HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
*/
    /// Make central and cloned nodes for rigid body
    void makeDofmanagers( InputRecord &ir );
    std::vector<FloatArray> coordsFromIr( InputRecord &ir );
    int makeDofmanager( InputRecord &dummyIr );
    int makeDofmanager( InputRecord &dummyIr, int id );
    int nextDofmanagerID();
    void rbsmDummyIr( InputRecord &irIn, std::vector<OOFEMTXTInputRecord> &irOut, int master );
    void setCornerNodesFromIr( InputRecord &ir );
    /**
     * Local renumbering support. For some tasks (parallel load balancing, for example) it is necessary to
     * renumber the entities.
     * Since RBSM Tetra class cannot access dofManLabelMap we have to use this override as a work around.
     * This should be amended once dofManLabelMap is a property of Domain class.
     * @param f is a functor that decides the renumbering will not be used in this override.
     */
    void updateLocalNumbering( EntityRenumberingFunctor &f ) override;
};
} // end namespace oofem
#endif // rbsmtetra_h
