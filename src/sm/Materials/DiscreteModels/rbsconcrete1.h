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
#define _IFT_RBSConcrete1_yieldstress "fc"
#define _IFT_RBSConcrete1_tangentmodulus "et" //Et
//#define _IFT_RBSConcrete1_hardeningmoduli "h"
//@}

namespace oofem {
class Domain;

/**
 * This class implements a isotropic plastic linear material (J2 plasticity condition is used).
 * @author Jim Brozoulis
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

    /// Initial (uniaxial) yield stress.
    double sig0 = 0.;

    IsotropicLinearElasticMaterial D;

public:
    RBSConcrete1(int n, Domain * d);

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &ir) override;
    const char *giveInputRecordName() const override { return _IFT_RBSConcrete1_Name; }
    const char *giveClassName() const override { return "RBSConcrete1"; }
    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return true; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override;

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

    int tempNormalState = 0;
    int normalState     = 0;

public:
    RBSConcrete1Status(GaussPoint * g);

    const FloatArrayF<6> &givePlasticStrain() const { return plasticStrain; }

    void letTempPlasticStrainBe(const FloatArrayF<6> &values) { tempPlasticStrain = values; }

    double giveK() const { return this->k; }

    double giveNormalState() const { return this->normalState; }

    void letTempKBe(double value) { tempK = value; }

    void letTempNormalStateBe(double value) { tempNormalState = value; }

    void letTempDevTrialStressBe(const FloatArrayF<6> &values) { tempDevTrialStress = values; }
    const FloatArrayF<6> &giveTempDevTrialStress() const { return tempDevTrialStress; }

    const char *giveClassName() const override { return "RBSConcrete1Status"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    // semi optional methods
    //void printOutputAt(FILE *file, TimeStep *tStep) override;
    //void saveContext(DataStream &stream, ContextMode mode) override;
    //void restoreContext(DataStream &stream, ContextMode mode) override;
};

} // end namespace oofem
#endif // rbsconcrete1_h
