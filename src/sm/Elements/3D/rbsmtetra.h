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
// .: properties :.
public:

protected:
    static FEI3dTetLin interpolation;

    /// map geometry (mesh) node number to cloned nodes local number
    static std::map<int, std::set<int>> clonesOfGeoNode;
    /// map geometry (mesh) node number to RBSM element local number
    static std::map<int, std::set<int>> cellElementsOfGeoNode;
    /// map set of "facet nodes set" to "elements" local number
    static std::map<std::vector<int>, std::set<int>> cellElementsOfFacets;
    /// number of vertices
    int numberOfCornerNodes; // make static constant
    /// global ID of rigid body cell central node
    int centerDofmanager;
    /// global number of mesh nodes defining corners of the rigid body
    IntArray geoNodes;
    /// global number of the springs beam elements
    std::vector<IntArray> springsBeams;
    /// facets indices of rigid body
    IntArray facetArray;


// .: methods :.
public:
    RBSMTetra(int n, Domain * d);
    virtual ~RBSMTetra() { }

    FEInterpolation *giveInterpolation() const override;

    /// @returns index numbers of clone nodes associated with the geometry node
    static std::set<int> giveClonesOfGeoNode(int geoNode);
    /// @returns index numbers of elements associated with the geometry node
    static std::set<int> giveCellElementsOfGeoNode(int geoNode);
    /// @returns index numbers of cell node of elements associated with the geometry node
    //static std::set<int> giveCellNodesOfGeoNode( int geoNode );

    /// @returns index number of rigid body cell central node
    int giveCellDofmanagerNumber() { return centerDofmanager; }

    void initializeFrom(InputRecord &ir) override;

    void setCrossSection(int csIndx) override;

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    int giveNumberOfIPForMassMtrxIntegration() override { return 4; }
    Interface *giveInterface(InterfaceType it) override;
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_RBSMTetra_Name; }
    const char *giveClassName() const override { return "RBSMTetra"; }

/*
#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawSpecial(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif
*/

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
    /**
     * Makes a DOF manager and returns assigned local number
     * @returns local number which is based on an arithmetic
     * progression starting from 1 and is assigned automatically.
     */
    int makeDofmanager( InputRecord &dummyIr );
    /**
     * Makes a DOF manager and returns assigned local number
     * @param globalNumber represents the number inside an input file.
     * @returns local number which is based on an arithmetic
     * progression starting from 1 and is assigned automatically.
     */
    int makeDofmanager( InputRecord &dummyIr, int globalNumber );
    /**
     * @returns next available global number for a new DOF manager,
     * global number represents the number given inside an input file.
     */
    int nextDofmanagerGlobalNumber();
    /**
     * @returns next available global number for a new element,
     * global number represents the number given inside an input file.
     */
    int nextElementGlobalNumber();

    /**
     * Makes springs beam element and returns assigned local number
     * @returns local number which is based on an arithmetic
     * progression starting from 1 and is assigned automatically.
     */
    int makeSpringsBeam( int globalNumber, int dmanA, int dmanB );


protected:
    void rbsmDummyIr( InputRecord &irIn, std::vector<OOFEMTXTInputRecord> &irOut, int master );
    /// set the geometry nodes numbers from input file
    void setGeoNodesFromIr( InputRecord &ir );
    /// update the geometry nodes to clone nodes map
    ///@param dofManArrayIsGlobal if true, means clone numbers should be converted to local
    void updateClonesOfGeoNode(bool dofManArrayIsGlobal = true);
    /// update the geometry nodes to center nodes map
    void updateCellElementsOfGeoNode();
    /// update the facet to element map
    std::map<int, IntArray> updateCellElementsOfFacets();
};
} // end namespace oofem
#endif // rbsmtetra_h
