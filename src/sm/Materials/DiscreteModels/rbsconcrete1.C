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

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#include "rbsconcrete1.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "sm/Elements/structuralelement.h"
#include "mathfem.h"
#pragma clang diagnostic pop

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

    this->shearCoef = 2.;
    this->G       = D.giveShearModulus();

    sigma_k.resize(maxNK);
    G_k.resize( maxNK );
    eps_k.resize( maxNK + 1 );
    epsP_k.resize( maxNK + 1 );
    H_k.resize( maxNK + 1 );

    // shear yield stress
    this->sigma_k = {
        this->sig0 / shearCoef * .333333,
        this->sig0 / shearCoef * .666667,
        this->sig0 / shearCoef,
        0.
    };

    ///TODO: this one should be inside the status since it depends on GP
    // nonlinear shear (hardening/softening) moduli
    this->G_k = { G, G / 2., G / 4., -G / 10. };
    // calculate other multilinear specifications
    eps_k( 0 )  = sigma_k( 0 ) / G_k( 0 );
    epsP_k( 0 ) = 0;
    for ( int i = 1; i < maxNK; ++i ) {
        eps_k( i )  = eps_k( i - 1 ) + ( sigma_k( i ) - sigma_k( i - 1 ) ) / G_k( i );
        epsP_k( i ) = eps_k( i  ) - sigma_k( i ) / G_k( 0 );
        if ( eps_k( i ) != eps_k( i ) || epsP_k( i ) != epsP_k( i ) ){
            OOFEM_ERROR("Unable to calculate concrete yield stress due to divide to zero error.")
        }
        H_k( i - 1 ) = ( sigma_k( i ) - sigma_k( i - 1 ) ) / ( epsP_k( i ) - epsP_k( i - 1 ) );
    }
    eps_k( maxNK )  = 1.e+16; // INFINITY
    epsP_k( maxNK ) = 1.e+16; // INFINITY
    H_k( maxNK ) = 0.;
    //
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

    // subtract plastic strains
    auto plasticStrain = status->givePlasticStrain();
    FloatArrayF<6> trialElasticStrain;
    for ( int i = 1; i <= 6; ++i ) {
        trialElasticStrain.at( i ) = strain.at( i ) - sgn( strain.at( i ) ) * plasticStrain.at( i );
    }

    const auto &elasticStiffness = D.giveTangent();
    auto trialStress = dot(elasticStiffness, trialElasticStrain );

    // Trial stresses
    double trialNormalStress = trialStress.at( 1 );
    double trialShearStress1 = trialStress.at( 5 );
    double trialShearStress2 = trialStress.at( 6 );

    //auto [devTrialStress, meanTrialStress] = computeDeviatoricVolumetricSplit(); // c++17
    //double J2 = this->computeSecondStressInvariant(devTrialStress);

    FloatArrayF<6> stress;
    stress = trialStress;

    /// Spring nonlinear curve state
    int nKn0        = status->giveNormalState();
    int nKs1        = status->giveShearState1();
    int nKs2        = status->giveShearState2();
    double k = status->giveK();
    double ks1 = status->giveKs1();
    double ks2 = status->giveKs2();

    // 1. Normal stress correction

    // evaluate the yield surface
    double sigma_y = this->sig0 + H * k;    ///fixme!
    double tr_f = trialNormalStress - sigma_y;

    if ( nKn0 ) {
        // damaged normal spring
        stress.at( 1 ) = 0;
    } else if ( tr_f <= 0.0 ) { // elastic
        //status->letTempPlasticStrainBe( status->givePlasticStrain() );
    } else { // plastic loading
        // plastic strain inc.
        double dPlStrain = tr_f / ( E + H );
        // radial return
        auto corNormalStress = trialNormalStress - E * dPlStrain;
        if ( corNormalStress < 0. ) {
            corNormalStress = 0.;
            nKn0            = 1;
        }
        stress.at( 1 ) = corNormalStress;
        k += dPlStrain;

        //auto plasticStrain = status->givePlasticStrain();
        plasticStrain.at( 1 ) += dPlStrain;
        //status->letTempPlasticStrainBe( plasticStrain );
    }


    // 2. Shear stress correction

    if ( status->giveNormalState() ) {
        // damaged normal spring
        stress.at( 5 ) = 0;
        stress.at( 6 ) = 0;
    } else {
        // Shear spring 1:
        {
            double sigma_ys1, tr_fs1;
            sigma_ys1 = sigma_k( nKs1 ) + H_k( nKs1 ) * ( ks1 - epsP_k( nKs1 ) );   ///fixme!
            tr_fs1    = fabs( trialShearStress1 ) - sigma_ys1;
            //
            if ( nKn0 ) { // if normal spring has failed, there is no shear stress
                stress.at( 5 ) = 0;
            } else if ( tr_fs1 <= 0.0 ) { // elastic
                // status->letTempPlasticStrainBe( status->givePlasticStrain() );
            } else { // plastic loading
                double dPlStrain = tr_fs1 / ( G + H_k( nKs1 ) ); // plastic multiplier
                // radial return
                auto corShearStress = trialShearStress1 - sgn( trialShearStress1 ) * G * dPlStrain;
                if ( G_k( nKs1 ) < 0 && corShearStress * trialShearStress1 < 0. ) {
                    corShearStress = 0.;
                    // shear failure cause normal failure
                    // normalState++;
                }
                ks1 += dPlStrain;
                plasticStrain.at( 5 ) += dPlStrain;

                // this can be a loop but it is deliberately avoided
                if ( ks1 > epsP_k( nKs1 + 1 ) ) {
                    nKs1++;
                    // Evaluate new yield stress from extra plastic strain
                    sigma_ys1 = sigma_k( nKs1 ) + H_k( nKs1 ) * ( ks1 - epsP_k( nKs1 ) );
                    tr_fs1    = fabs( corShearStress ) - sigma_ys1;
                    // Calculate complimentary part of plastic strain
                    dPlStrain      = tr_fs1 / ( G + H_k( nKs1 ) );
                    corShearStress = corShearStress - sgn( corShearStress ) * G * dPlStrain;
                    if ( sgn( trialShearStress1 ) * sgn( corShearStress ) < 0 ) {
                        OOFEM_WARNING( "Invalid shear stress for spring S1" );
                    }

                    ks1 += dPlStrain;
                    plasticStrain.at( 5 ) += dPlStrain;
                }

                stress.at( 5 ) = corShearStress;
            }
        }
        // Shear 2
        {
            double sigma_ys2, tr_fs2;
            sigma_ys2 = sigma_k( nKs2 ) + H_k( nKs2 ) * ( ks2 - epsP_k( nKs2 ) );
            tr_fs2    = fabs( trialShearStress2 ) - sigma_ys2;
            //
            if ( nKn0 ) { // if normal spring has failed, there is no shear stress
                stress.at( 6 ) = 0;
            } else if ( tr_fs2 <= 0.0 ) { // elastic
                // status->letTempPlasticStrainBe( status->givePlasticStrain() );
            } else { // plastic loading
                double dPlStrain = tr_fs2 / ( G + H_k( nKs2 ) ); // plastic multiplier
                // radial return
                auto corShearStress = trialShearStress2 - sgn( trialShearStress2 ) * G * dPlStrain;
                if ( G_k( nKs2 ) < 0 && corShearStress * trialShearStress2 < 0. ) {
                    corShearStress = 0.;
                    // shear failure cause normal failure
                    // normalState++;
                }
                ks2 += dPlStrain;
                plasticStrain.at( 6 ) += dPlStrain;

                // this can be a loop but it is deliberately avoided
                if ( ks2 > epsP_k( nKs2 + 1 ) ) {
                    nKs2++;
                    // Evaluate new yield stress from extra plastic strain
                    sigma_ys2 = sigma_k( nKs2 ) + H_k( nKs2 ) * ( ks2 - epsP_k( nKs2 ) );
                    tr_fs2    = fabs( corShearStress ) - sigma_ys2;
                    // Calculate complimentary part of plastic strain
                    dPlStrain      = tr_fs2 / ( G + H_k( nKs2 ) );
                    corShearStress = corShearStress - sgn( corShearStress ) * G * dPlStrain;
                    if ( sgn( trialShearStress2 ) * sgn( corShearStress ) < 0 ) {
                        OOFEM_ERROR( "Invalid shear stress for spring S2" );
                    }

                    ks2 += dPlStrain;
                    plasticStrain.at( 6 ) += dPlStrain;
                }

                stress.at( 6 ) = corShearStress;
            }
        }
    }

    // 3. Store the temporary values for the given iteration
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(stress);
    status->letTempPlasticStrainBe( plasticStrain );
    status->letTempKBe(k);
    status->letTempKs1Be(ks1);
    status->letTempShearState1(nKs1);
    status->letTempKs2Be(ks2);
    status->letTempShearState2(nKs2);
    //status->letTempDevTrialStressBe(devTrialStress);
    status->letTempNormalStateBe( nKn0 );
    return stress;
}


FloatMatrixF<6,6>
RBSConcrete1::give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< RBSConcrete1Status * >( this->giveStatus(gp) );
    double trialNormalStress    = status->giveStressVector().at( 1 );
    double trialShearStress1    = status->giveStressVector().at( 5 );
    double trialShearStress2    = status->giveStressVector().at( 6 );

    auto elasticStiffness = D.giveTangent();

    double k = status->giveK();
    double sigma_y = this->sig0 + H * k;

    // 1. Normal spring

    // evaluate the yield surface
    double tr_f = trialNormalStress - sigma_y;

    if(status->giveNormalState()) {
        // prevent instability by not using 0.
        elasticStiffness.at( 1, 1 ) = 1.; // almost zero
    } else if ( tr_f < 0.0 ) {
        // elastic loading
    } else {
        // plastic loading: Et
        // 1. use Et > 0 or small positive value
        /*elasticStiffness.at(1, 1 ) = max(this->Et,
            min( 1., elasticStiffness.at( 1, 1 ) / 10000. ));*/
        // 2. use 1.
        /*elasticStiffness.at(1, 1 ) = 1.;*/
        // 3. use Et
        /*elasticStiffness.at(1, 1 ) = this->Et;*/
        // 4. use Et > 0 or Ee
        elasticStiffness.at(1, 1 ) = this->Et > 0. ? this->Et : this->E;
    }

    // 2. Shear springs

    if ( status->giveNormalState() ) {
        // damaged normal spring
        elasticStiffness.at( 5, 5 ) = 1.;   // almost zero
        elasticStiffness.at( 6, 6 ) = 1.;   // almost zero
    } else {
        int nKs1        = status->giveShearState1();
        int nKs2        = status->giveShearState2();
        double ks1 = status->giveKs1();
        double ks2 = status->giveKs2();
        //double sigma_ys1 = this->sig0 + Hs1 * ks1;
        double sigma_ys1 = sigma_k( nKs1 ) + H_k( nKs1 ) * ( ks1 - epsP_k( nKs1 ) );       ///FIXME!
        double tr_fs1 = fabs( trialShearStress1 ) - sigma_ys1;
        double sigma_ys2 = sigma_k( nKs2 ) + H_k( nKs2 ) * ( ks2 - epsP_k( nKs2 ) );       ///FIXME!
        double tr_fs2 = fabs( trialShearStress2 ) - sigma_ys2;

        // Shear 1
        if ( tr_fs1 <= 0.0 ) { // elastic loading
            // continue
        } else { // plastic loading;
            // use Gt for hardening (faster convergence) and elastic G for softening (avoid divergence)
            elasticStiffness.at( 5, 5 ) = this->G_k(nKs1) > 0. ? this->G_k(nKs1) : G;
        }

        // Shear 2
        if ( tr_fs2 <= 0.0 ) { // elastic loading
            // continue
        } else { // plastic loading
            elasticStiffness.at( 6, 6 ) = this->G_k(nKs2) > 0. ? this->G_k(nKs2) : G;
        }
    }

    return elasticStiffness;
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
    tempNormalState = normalState;
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
    k   = tempK;
    ks1 = tempKs1;
    ks2 = tempKs2;
    normalState   = tempNormalState;
    shearState1   = tempShearState1;
    shearState2   = tempShearState2;
    // deviatoric trial stress is not really a state variable and was used not to repeat some code...
}

} // end namespace oofem
