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

#include "rbsconcrete1.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "sm/Elements/structuralelement.h"
#include "mathfem.h"

namespace oofem {
REGISTER_Material(RBSConcrete1);

RBSConcrete1::RBSConcrete1(int n, Domain *d) : StructuralMaterial(n, d), D(n, d)
{}


void
RBSConcrete1::initializeFrom(InputRecord &ir)
{
    StructuralMaterial::initializeFrom(ir);

    D.initializeFrom(ir);
    IR_GIVE_FIELD(ir, this->sig0, _IFT_RBSConcrete1_yieldstress);
    IR_GIVE_FIELD(ir, this->Et, _IFT_RBSConcrete1_tangentmodulus);

    this->E = this->D.giveYoungsModulus();
    this->H = E * Et / ( E - Et );
}


void RBSConcrete1::giveInputRecord(DynamicInputRecord &ir)
{
    StructuralMaterial::giveInputRecord(ir);
    D.giveInputRecord(ir);
    ir.setField(this->sig0, _IFT_RBSConcrete1_yieldstress);
    ir.setField(this->Et, _IFT_RBSConcrete1_tangentmodulus);
    //ir.setField(this->H, _IFT_RBSConcrete1_hardeningmoduli);
}


MaterialStatus *RBSConcrete1::CreateStatus(GaussPoint *gp) const
{
    return new RBSConcrete1Status(gp);
}


FloatArrayF<6>
RBSConcrete1::giveRealStressVector_3d( const FloatArrayF<6> &totalStrain, GaussPoint *gp, TimeStep *tStep ) const
{
    auto status = static_cast<RBSConcrete1Status *>( this->giveStatus( gp ) );

    // subtract stress thermal expansion
    auto thermalStrain = this->computeStressIndependentStrainVector_3d( gp, tStep, VM_Total );
    auto strain = totalStrain - thermalStrain;

    auto trialElastStrain = strain - status->givePlasticStrain();

    const auto &elasticStiffness = D.giveTangent();
    auto trialStress = dot(elasticStiffness, trialElastStrain);

    //auto [devTrialStress, meanTrialStress] = computeDeviatoricVolumetricSplit(); // c++17
    //auto tmp = computeDeviatoricVolumetricSplit(trialStress);
    //auto devTrialStress = tmp.first;
    //auto meanTrialStress = tmp.second;

    //double J2 = this->computeSecondStressInvariant(devTrialStress);

    double trialNormalStress = trialStress.at( 1 );

    // evaluate the yield surface
    double k = status->giveK();
    double sigma_y = this->sig0 + H * k;
    double tr_f = trialNormalStress - sigma_y;

    FloatArrayF<6> stress;
    if ( tr_f <= 0.0 ) { // elastic
        stress = trialStress;

        status->letTempPlasticStrainBe(status->givePlasticStrain());
    } else { // plastic loading
        double E = D.giveYoungsModulus();
        double dPlStrain = tr_f / ( E + H ); // plastic multiplier
        // radial return
        auto corNormalStress = trialNormalStress - E * dPlStrain;
        stress = trialStress;
        stress.at( 1 ) = corNormalStress;
        k += dPlStrain;

        auto plasticStrain = status->givePlasticStrain();
        plasticStrain.at(1) += dPlStrain;
        status->letTempPlasticStrainBe(plasticStrain);
    }

    // Store the temporary values for the given iteration
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(stress);
    status->letTempKBe(k);
    //status->letTempDevTrialStressBe(devTrialStress);
    return stress;
}


FloatMatrixF<6,6>
RBSConcrete1::give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< RBSConcrete1Status * >( this->giveStatus(gp) );
    double trialNormalStress    = status->giveStressVector().at( 1 );

//    const auto &devTrialStress = status->giveTempDevTrialStress();
//    double J2 = this->computeSecondStressInvariant(devTrialStress);
//    double effectiveTrialStress = sqrt(3 * J2);


    // evaluate the yield surface
    double k = status->giveK();
    double sigma_y = this->sig0 + H * k;
    double tr_f = trialNormalStress - sigma_y;
    
    auto elasticStiffness = D.giveTangent();

    if ( tr_f < 0.0 ) { // elastic
        return elasticStiffness;
    } else { // plastic loading
        elasticStiffness.at(1, 1) = this->Et; // set Et ***
        return elasticStiffness;
    }
}


FloatArrayF<6>
RBSConcrete1::giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const
{
    return D.giveAlpha();
}


int
RBSConcrete1::giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    RBSConcrete1Status *status = static_cast< RBSConcrete1Status * >( this->giveStatus(gp) );
    if ( type == IST_PlasticStrainTensor ) {
        answer = status->givePlasticStrain();
        return 1;
    } else {
        return StructuralMaterial::giveIPValue(answer, gp, type, tStep);
    }
}


//=============================================================================


RBSConcrete1Status::RBSConcrete1Status(GaussPoint * g) :
    StructuralMaterialStatus(g)
{
    strainVector.resize(6);
    stressVector.resize(6);
    tempStressVector = stressVector;
    tempStrainVector = strainVector;
}

void
RBSConcrete1Status::initTempStatus()
{
    //StructuralMaterialStatus::initTempStatus();

    // reset temp vars.
    tempStressVector = stressVector;
    tempStrainVector = strainVector;
    tempPVector      = PVector;
    tempFVector      = FVector;


    tempPlasticStrain = plasticStrain;
    tempK = k;
    tempDevTrialStress = zeros<6>();
}


void 
RBSConcrete1Status::updateYourself(TimeStep *tStep)
{
    // Copy the temporary values to the convered ones. 
    // This method is called after equilibrium has been reached and should 
    // set all...
    StructuralMaterialStatus::updateYourself(tStep);

    plasticStrain = tempPlasticStrain;
    k = tempK;
    // deviatoric trial stress is not really a state variable and was used not to repeat some code...
}

} // end namespace oofem
