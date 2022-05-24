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

#ifndef rbsconcrete1_h
#define rbsconcrete1_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/isolinearelasticmaterial.h"

///@name Input fields for RBSConcrete1
//@{
#define _IFT_RBSConcrete1_Name "rbsconcrete1"
#define _IFT_RBSConcrete1_fc "fc"
#define _IFT_RBSConcrete1_tangentmodulus "et" //Et
//#define _IFT_RBSConcrete1_hardeningmoduli "h"
//@}

namespace oofem {
class Domain;

/**
 * This class implements a nonlinear concrete-like material for RBSM element
 * @author Saeid Mehrpay
 */
class RBSConcrete1 : public StructuralMaterial
{
protected:
    /// Elastic modulus.
    double E = 0.;
    /// tangent modulus.
    double Et = 0.;
    /// plastic modulus.
    double H = 0.;

    /// Poisson's ratio
    double nu = 0.;

    /// Peak (uniaxial) compressive yield stress.
    double fc = 0.;
    /// Peak tensile stress.
    double ft = 0.;
    /// Peak shear stress.
    double fs = 0.;

    //  *** M U L T I L I N E A R  S H E A R ***
    /// number of multilinear stages
    static const int maxNKShear = 4, maxNKTensile = 2;
    /// number of final crack mouth opening stages
    static const int maxNCmdS = 2, maxNCmdT = 1;
    // these can be constants
    /// linear stress coefficient before f'c (between 0. and 1.)
    double linearStressRatio = 0.;
    /// shear coefficient
    double shearCoef = 0.;
    // macroscopic compressive strain corresponding to peak
    double criticalStrain = 0.;
    /// shear spring strain corresponding to peak
    double criticalShearStrain = 0.;
    /// tensile spring strain corresponding to peak
    double criticalTensileStrain = 0.;

    /// todo: should be material input properties.
    /// shear displacement cracks in mm //***
    double tensile_cmd1 = 0.3, shear_cmd1 = 0.3, shear_cmd2 = 1.;
    /// critical strain and shear crack mouth openings [mm] stage 1 & 2
    FloatArray shearCmdKeyPoints;
    /// critical strain and tensile crack mouth openings [mm] stage 1
    FloatArray tensileCmdKeyPoints;

    /// elastic shear modulus
    double G;
    /// shear stresses (sigma) & hardening (Gt)
    FloatArray fs_k, G_k;
    //FloatArrayF<maxNK> sigma_k, G_k;
    /// strains corresponding to strain-strain key-points
    FloatArray eps_k, epsP_k;
    /// plastic modulus (H)
    FloatArray H_k;
    //FloatArrayF<maxNK+1> eps_k, epsP_k;
    //FloatArrayF<maxNK> H_k;


    IsotropicLinearElasticMaterial D;

public:
    RBSConcrete1(int n, Domain * d);

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &ir) override;
    const char *giveInputRecordName() const override { return _IFT_RBSConcrete1_Name; }
    const char *giveClassName() const override { return "RBSConcrete1"; }
    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return true; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    FloatArrayF<6> giveRealStressVector_3d( const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep ) const override;
    FloatArrayF<6> giveRealStressVector_3dBeamSubSoil( const FloatArrayF<6> &reducedStrain, GaussPoint *gp, TimeStep *tStep ) const override;
    /// Overrides the default implementation to preserve the stresses induced due to confinement
    FloatArray giveRealStressVector_StressControl(const FloatArray &reducedE, const IntArray &strainControl, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    FloatArrayF<6> giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const override;
};


class RBSConcrete1Status : public StructuralMaterialStatus
{
protected:
    /// Temporary plastic strain (the given iteration)
    FloatArrayF<6> tempPlasticStrain;

    ///  Last equilibriated plastic strain (end of last time step)
    FloatArrayF<6> plasticStrain;

    FloatArrayF<6> tempDevTrialStress;

    double tempK = 0.;
    double k = 0.;

    double tempKs1 = 0.;
    double ks1 = 0.;

    double tempKs2 = 0.;
    double ks2 = 0.;

    int tempNormalState = 0;
    int normalState     = 0;

    int tempShearState1 = 0;
    int shearState1     = 0;

    int tempShearState2 = 0;
    int shearState2     = 0;

    // Crack Mouth Displacement to Strain Conversion factor
    double cmdToStrain;

public:
    RBSConcrete1Status(GaussPoint * g);

    const FloatArrayF<6> &givePlasticStrain() const { return plasticStrain; }

    void letTempPlasticStrainBe(const FloatArrayF<6> &values) { tempPlasticStrain = values; }

    double giveK() const { return this->k; }

    double giveKs1() const { return this->ks1; }

    double giveKs2() const { return this->ks2; }

    double giveNormalState() const { return this->normalState; }

    double giveShearState1() const { return this->shearState1; }

    double giveShearState2() const { return this->shearState2; }

    void letTempKBe(double value) { tempK = value; }

    void letTempKs1Be(double value) { tempKs1 = value; }

    void letTempKs2Be(double value) { tempKs2 = value; }

    void letTempNormalStateBe(double value) { tempNormalState = value; }

    void letTempShearState1(int value) { tempShearState1 = value; }

    void letTempShearState2(int value) { tempShearState2 = value; }

    void letTempDevTrialStressBe(const FloatArrayF<6> &values) { tempDevTrialStress = values; }
    const FloatArrayF<6> &giveTempDevTrialStress() const { return tempDevTrialStress; }

    const char *giveClassName() const override { return "RBSConcrete1Status"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    // semi optional methods
    //void printOutputAt(FILE *file, TimeStep *tStep) override;
    //void saveContext(DataStream &stream, ContextMode mode) override;
    //void restoreContext(DataStream &stream, ContextMode mode) override;

    /**
     * Updates multi-linear parameters for the cracked stage where slip displacement
     * is used instead of strain
     * @param nKS Number of multilinear stages for tensile behavior
     * @param nCmdS number of shear slip crack displacement key points
     * @param eps_K array of strain key points
     * @param epsP_K array of plastic strain key points
     * @param G_k array of tangent shear moduli
     * @param H_k array of plastic moduli
     * @param fs_k array of stress key points
     * @param shearCmdKeyPoints array of shear crack slip key displacements
     */
    void updateMaterialCrackSlipParameters(
        int nKS, int nCmdS,
        FloatArray &eps_k, FloatArray &epsP_k, FloatArray &G_k, FloatArray &H_k,
        FloatArray fs_k, FloatArray shearCmdKeyPoints );

    /**
     * Updates multilinear parameters for the cracked stage where opening
     * displacement is used instead of strain
     * @param nKT Number of multilinear stages for shear behavior
     * @param nCmdT number of tensile opening displacements key points
     * @param H plastic modulus
     * @param E Young modulus
     * @param ft tensile strength
     * @param tensCritStrain tensile strain at peak stress
     * @param tensileCmdU ultimate tensile crack mouth opening
     */
    void updateMaterialCrackOpeningParameters(
        int nKT, int nCmdT,
        double &H,
        double E, double ft, double tensCritStrain, double tensileCmdU );
};

} // end namespace oofem
#endif // rbsconcrete1_h
