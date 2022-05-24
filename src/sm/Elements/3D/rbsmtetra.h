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
#include "sm/Materials/DiscreteModels/rbsitz.h"
#include "oofemtxtinputrecord.h"

#define _IFT_RBSMTetra_Name "rbsmtetra"
#define _IFT_RBSMTetra_itz "itz"
#define MINDLIN

namespace oofem {
class FEI3dTetLin;
class RBSMBeam3d;
class RBSBeam3d;

/**
 * This class implements a tetrahedral rigid body - spring model.
 * Each element has one master node with 6 degree of freedom.
 * @author Saeid Mehrpay
 */
class RBSMTetra : public Structural3DElement, public ItzInterface
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
    /// facets' stresses used by rigid body todo: make private
    std::vector<FloatArray> tempFacetsStressVector;

protected:
    static FEI3dTetLin interpolation;

    // Added to this class to avoid domain modification todo: actually belong to domain
    /// missing domain class properties as a workaround to avoid modifying domain class
    static int domain_nDofman, domain_nElements, domain_buffer;
    /// largest global number assigned to an existing DOF manager
    static int domain_maxDofGlNum;

    /// last reported preprocessing progress
    static int prepProgress;

    /// map geometry (mesh) node number to cloned nodes local number
    static std::map<int, std::set<int>> clonesOfGeoNode;
    /// map geometry (mesh) node number to RBSM element local number
    static std::map<int, std::set<int>> cellElementsOfGeoNode;
    /// map set of "facet nodes set" to "elements" local number
    static std::map<std::vector<int>, std::set<int>> mapFacetElement;
    /// number of vertices
    int numberOfCornerNodes; // make static constant
    /// global ID of rigid body cell central node
    int centerDofmanager;
    // Interfacial material
    //int itzMaterial = 0;
    /// global number of mesh nodes defining corners of the rigid body
    IntArray geoNodes;
    /// map facet number to array of existing sisters
    std::map<int, IntArray> existingSisters;
    /// global number of the springs beam elements
    std::vector<IntArray> springsBeams;
    /// facets' vertex indices of rigid body
    std::vector<IntArray> facetArray;


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

    /// @returns index number of geometry nodes
    IntArray giveGeoNodes() { return geoNodes; }

    /**
     * calculates distance of the fiber zone centroid from facet centroid
     * @param nFacet the number of targeted facet of the element
     * @return distance coordinates for all fiber zones in a vector
     */
    std::vector<FloatArray> giveFiberZonesOffsetsOfFacet_3TriaDissect( int nFacet );

    /**
     * calculates confining stress of facet
     * @param nFacet the number of targeted facet of the element
     * @return confining stress vector
     */
    FloatArray giveConfiningStress( int nFacet, TimeStep *tStep );

    void initializeFrom(InputRecord &ir) override;

    void postInitialize() override;

    void updateLocalNumbering(EntityRenumberingFunctor &f) override;

    void setCrossSection(int csIndx) override;


    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;

    double computeLength() override { return 0.; };

    int giveNumberOfIPForMassMtrxIntegration() override { return 4; }

    void giveInternalForcesVector(FloatArray &answer, TimeStep *, int useUpdatedGpRecord = 0) override;


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

    /**
     * calculates the area of facet
     * @param nFacet the target facet number
     */
    double giveAreaOfFacet(int nFacet);

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

    /**
     * Make central and cloned nodes for rigid body
     */
    void makeDofmanagers( InputRecord &ir );

    /**
     * Obtains the coordinates of rigid body cell
     * @param ir input record
     * @return array containing coordinates of the central node [0]
     * and corner nodes [1-4] of the rigid body cell
     */
    std::vector<FloatArray> coordsFromIr( InputRecord &ir );

    /**
     * Makes a DOF manager and returns assigned local number
     * @returns local number which is based on an arithmetic
     * progression starting from 1 and is assigned automatically.
     */
    int makeDofmanager( InputRecord &dummyIr,
                        const bool domainDofListResize = true, const int resizeExtraRoom = 0 );
    /**
     * Makes a DOF manager and returns assigned local number
     * @param globalNumber represents the number inside an input file.
     * @param Number represents the number of DOF manager.
     * @param domainDofListResize increase size of DOF managers list by 1 (default true).
     * @param resizeExtraRoom extra increment of DOF list (default 0).
     * @returns local number which is based on an arithmetic
     * progression starting from 1 and is assigned automatically.
     */
    int makeDofmanager( InputRecord &dummyIr, const int globalNumber, const int number,
                        const bool domainDofListResize = true, const int resizeExtraRoom = 0 );
    /**
     * @param nDofmanPlus increase number of DOF managers by (default 1).
     * @param skipNumber skip available global DOF managers by (default 0).
     * @param number stores next available DOF number plus skip number.
     * @returns next available global number plus skip number for a new DOF manager,
     * global number represents the number given inside an input file.
     */
    int nextDofmanagerGlobalNumber( int &number, int nDofmanPlus = 1, int skipNumber = 0 );
    /**
     * finds maximum global number used for a new Dof manager
     * @returns maximum global number of all existing DOFs,
     * global number represents the number given inside an input file.
     */
    int findMaxDofmanagerGlobalNumber();
    /**
     * checks global number can be used for a new DOF manager,
     * all DOF managers must be initialized
     * @returns true for a valid global number for a new DOF manager,
     * global number represents the number given inside an input file.
     */
    bool isValidDofmanagerGlobalNumber( int globalNum );
    /**
     * @returns next available global number for a new element,
     * global number represents the number given inside an input file.
     */
    int nextElementGlobalNumber( int baseNumber = 1 );

    /**
     * Makes springs beam element and returns assigned local number
     * @returns local number which is based on an arithmetic
     * progression starting from 1 and is assigned automatically.
     */
    int makeSpringsBeam( int globalNumber, int dmanA, int dmanB );

    /**
     * Makes springs beam element cross-section and returns assigned number
     * @param globalNumber
     * @param nFacet the number of facet whose cross-section will be made
     * @return assigned number to the new cross-section
     */
    int makeSpringsBeamCrossSection_3TriaDissect( int nFacet );
    int makeSpringsBeamCrossSection_CircularDist( int nFacet, int numberOfFibers, int materialNumber );
    int findSpringsBeamMaterial( int nFacet );

    void giveCharacteristicMatrix( FloatMatrix &answer, CharType type, TimeStep *tStep ) override;
    void giveCharacteristicVector( FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep ) override;

    //bool isActivated( TimeStep *tStep ) override;

protected:
    /**
     * @param irIn input record to obtain node coordinates.
     * @param irOut generated dummy input record for making RBSM nodes including rigid arm nodes.
     * @param masterGlNum global number of master node to tie rigid arms too.
     * global number represents the number given inside an input file.
     */
    void rbsmDummyIr( InputRecord &irIn, std::vector<OOFEMTXTInputRecord> &irOut, int masterGlNum );
    /// set the geometry nodes numbers from input file
    void setGeoNodesFromIr( InputRecord &ir );
    /// update the geometry nodes to clone nodes map
    ///@param dofManArrayIsGlobal if true, means clone numbers should be converted to local
    void updateClonesOfGeoNode(bool dofManArrayIsGlobal = true);
    /// update the geometry nodes to center nodes map
    void updateCellElementsOfGeoNode();
    /// update the facet to element map
    ///@return map of existing sister element for each facet
    std::map<int, IntArray> updateCellElementsOfFacets();
#ifdef MINDLIN
    RBSMBeam3d *giveSpringsBeamOfFacet( int nFacet );
#else
    RBSBeam3d *giveSpringsBeamOfFacet( int nFacet );
#endif
};

/**
 * The element interface required by RBSMTetra.
 */
//class RBSMTetraInterface : public Interface
//{
//public:
//    RBSMTetraInterface() { }
//
//    /**
//     * Computes full 3d strain vector in element fiber. This function is necessary
//     * if layered cross section is specified.
//     * @param answer Full fiber strain vector.
//     * @param masterGpStrain Strain vector at master gauss point.
//     * @param slaveGp Slave integration point representing particular fiber.
//     * @param tStep Time step.
//     */
//    virtual void RBSMTetraInterface_computeConfinedStressVector(FloatArray &answer, const FloatArray &masterGpStrain,
//        GaussPoint *slaveGp, TimeStep *tStep) = 0;
//};
} // end namespace oofem
#endif // rbsmtetra_h
