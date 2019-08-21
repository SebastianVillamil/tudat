/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/altimeterCrossoverPartial.h"

namespace tudat
{

namespace observation_partials
{

////! Update the scaling object to the current times and states
//void AltimeterCrossoverScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
//                                 const std::vector< double >& times,
//                                 const observation_models::LinkEndType fixedLinkEnd,
//                                 const Eigen::VectorXd currentObservation )
//{
//    // Compute Euclidean distance vector
//    Eigen::Vector3d rangeVector = linkEndStates[ 1 ].segment( 0, 3 ) - linkEndStates[ 0 ].segment( 0, 3 );
//    Eigen::Matrix< double, 1, 3 > rangeVectorNormalized = rangeVector.transpose( ) / rangeVector.norm( );

//    // Compute scaling for receiver reference
//    if( fixedLinkEnd == observation_models::receiver )
//    {
//        referenceLightTimeCorrectionScaling_ = 1.0 / ( 1.0 - rangeVectorNormalized.transpose( ).dot( linkEndStates[ 0 ].segment( 3, 3 ) ) /
//                physical_constants::SPEED_OF_LIGHT );
//        referenceScalingFactor_ =  rangeVectorNormalized * referenceLightTimeCorrectionScaling_;
//    }

//    // Compute scaling for transmitter reference
//    else if( fixedLinkEnd == observation_models::transmitter )
//    {
//        referenceLightTimeCorrectionScaling_ =
//                1.0 / ( 1.0 - rangeVectorNormalized.transpose( ).dot( linkEndStates[ 1 ].segment( 3, 3 ) ) /
//                physical_constants::SPEED_OF_LIGHT );
//        referenceScalingFactor_ =  rangeVectorNormalized * referenceLightTimeCorrectionScaling_;
//    }

//    currentLinkEndType_ = fixedLinkEnd;
//}

//! Function to calculate the observation partial(s) at required time and state
AltimeterCrossoverPartial::AltimeterCrossoverPartialReturnType AltimeterCrossoverPartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const Eigen::Vector1d& currentObservation )
{
    AltimeterCrossoverPartialReturnType returnPartial;

    Eigen::Matrix< double, 1, 3 > firstArcPartialWrtCurrentPosition;
    Eigen::Matrix< double, 1, 3 > secondArcPartialWrtCurrentPosition;

    // Definitions for dh/dr1:
    double dR2dt2;
    Eigen::Matrix< double, 1, 3 > dt2dr1;
    Eigen::Matrix< double, 1, 3 > dR1dr1;
    double dR1dt1;
    Eigen::Matrix< double, 1, 3 > dt1dr1;

    // Definitions for dh/dr2:
    Eigen::Matrix< double, 1, 3 > dR2dr2;
//    double dR2dt2;
    Eigen::Matrix< double, 1, 3 > dt2dr2;
//    double dR1dt1;
    Eigen::Matrix< double, 1, 3 > dt1dr2;

    Eigen::MatrixXd observationPartialWrtCurrentState = Eigen::MatrixXd::Zero( 1, 6 );

    // Iterate over all link ends
    for( positionPartialIterator_ = positionPartialList_.begin( ); positionPartialIterator_ != positionPartialList_.end( );
         positionPartialIterator_++ )
    {
        // Compute dhdr1
        if( positionPartialIterator_->first == observation_models::first_arc_body )
        {
            currentStateArc1_  = states[ 0 ];
            currentTimeArc1_ = times[ 0 ];
            currentStateArc2_  = states[ 1 ];
            currentTimeArc2_ = times[ 1 ];

            Eigen::Matrix< double, 3, Eigen::Dynamic > currentInertialPositionPartialWrtParameter =
                    positionPartialIterator_->second->calculatePartialOfPosition(
                                          currentStateArc1_ , currentTimeArc1_ );

            // Obtain body-fixed (bf) states
            Eigen::Vector6d bf_StateArc1 = transformStateToTargetFrame(
                        currentStateArc1_, currentTimeArc1_, centralBodyRotationModel_ );
            Eigen::Vector6d bf_StateArc2 = transformStateToTargetFrame(
                        currentStateArc2_, currentTimeArc2_, centralBodyRotationModel_ );

            // Definition of dR2dt2
            Eigen::Vector3d r2UnitVector = currentStateArc2_.segment( 0, 3 ) / currentStateArc2_.segment( 0, 3 ).norm();
            dR2dt2 = currentStateArc2_.segment( 3, 3 ).dot( r2UnitVector );

            // Definition of dt2dr1
            Eigen::Vector3d bf_r2UnitVector = bf_StateArc2.segment( 0, 3 ) / bf_StateArc2.segment( 0, 3 ).norm();
            Eigen::Vector3d bf_v1UnitVector = bf_StateArc1.segment( 3, 3 ) / bf_StateArc1.segment( 3, 3 ).norm();
            Eigen::Vector3d bf_v2UnitVector = bf_StateArc2.segment( 3, 3 ) / bf_StateArc2.segment( 3, 3 ).norm();
            Eigen::Matrix3d A_dt2dr1;
            A_dt2dr1 << bf_r2UnitVector.transpose(), bf_v1UnitVector.transpose(), bf_v2UnitVector.transpose();
            Eigen::Vector3d v2Vertical = (bf_StateArc2.segment( 3, 3 ).dot(bf_r2UnitVector))*bf_r2UnitVector;
            Eigen::Vector3d v2InPlane = bf_StateArc2.segment( 3, 3 ) - v2Vertical;
            Eigen::Vector3d b_dt2dr1( 0, 0, ( 1 / v2InPlane.norm( ) ) );
            Eigen::Vector6d bf_dt2dr1;
            bf_dt2dr1 << ( A_dt2dr1.inverse( ) * b_dt2dr1 ), 0, 0, 0;
            dt2dr1 << transformStateToGlobalFrame(
                          bf_dt2dr1, currentTimeArc1_, centralBodyRotationModel_ ).segment( 0, 3 ).transpose();

            // Definition of dR1dr1
            double r1_norm = currentStateArc1_.segment( 0, 3 ).norm( );
            dR1dr1 << ( currentStateArc1_.segment( 0, 3 ) / r1_norm ).transpose();

            // Definition of dR1dt1
            Eigen::Vector3d r1UnitVector = currentStateArc1_.segment( 0, 3 ) / currentStateArc1_.segment( 0, 3 ).norm();
            dR1dt1 = currentStateArc1_.segment( 3, 3 ).dot( r1UnitVector );

            // Definition of dt1dr1
            Eigen::Vector3d bf_r1UnitVector = bf_StateArc1.segment( 0, 3 ) / bf_StateArc1.segment( 0, 3 ).norm();
            Eigen::Matrix3d A_dt1dr1;
            A_dt1dr1 << bf_r1UnitVector.transpose(), bf_v2UnitVector.transpose(), bf_v1UnitVector.transpose();
            Eigen::Vector3d v1Vertical = (bf_StateArc1.segment( 3, 3 ).dot(bf_r1UnitVector))*bf_r1UnitVector;
            Eigen::Vector3d v1InPlane = bf_StateArc1.segment( 3, 3 ) - v1Vertical;
            Eigen::Vector3d b_dt1dr1( 0, 0, - ( 1 / v1InPlane.norm( ) ) );
            Eigen::Vector6d bf_dt1dr1;
            bf_dt1dr1 << ( A_dt1dr1.inverse( ) * b_dt1dr1 ), 0, 0, 0;
            dt1dr1 << transformStateToGlobalFrame(
                          bf_dt1dr1, currentTimeArc1_, centralBodyRotationModel_ ).segment( 0, 3 ).transpose();

            firstArcPartialWrtCurrentPosition << ( dR2dt2*dt2dr1 - dR1dr1 - dR1dt1*dt1dr1 );
//            firstArcPartialWrtCurrentPosition << ( - dR1dr1 );

            observationPartialWrtCurrentState.block( 0, 0, 1, 3 ) = firstArcPartialWrtCurrentPosition;

            returnPartial.push_back(
                        std::make_pair( observationPartialWrtCurrentState, currentTimeArc1_ ) );

//            std::cout << "\n (dR2dt2): \n" << dR2dt2 << std::endl;
//            std::cout << "\n (dt2dr1): \n" << dt2dr1 << std::endl;
//            std::cout << "\n dR2dt2*dt2dr1: \n" << dR2dt2*dt2dr1 << std::endl;

//            std::cout << "\n dR1dr1: \n" << dR1dr1 << std::endl;

//            std::cout << "\n (dR1dt1): \n" << dR1dt1 << std::endl;
//            std::cout << "\n (dt1dr1): \n" << dt1dr1 << std::endl;
//            std::cout << "\n dR1dt1*dt1dr1: \n" << dR1dt1*dt1dr1 << std::endl;

//            std::cout << "\n firstArcPartialWrtCurrentPosition: \n" << firstArcPartialWrtCurrentPosition << std::endl;
//            std::cout << "\n observationPartialWrtCurrentState: \n" << observationPartialWrtCurrentState << std::endl;
        }
        // Compute dhdr2
        else if( positionPartialIterator_->first == observation_models::second_arc_body )
        {
            currentStateArc1_  = states[ 0 ];
            currentTimeArc1_ = times[ 0 ];
            currentStateArc2_  = states[ 1 ];
            currentTimeArc2_ = times[ 1 ];

            Eigen::Matrix< double, 3, Eigen::Dynamic > currentInertialPositionPartialWrtParameter =
                    positionPartialIterator_->second->calculatePartialOfPosition(
                                          currentStateArc2_ , currentTimeArc2_ );

            // Obtain body-fixed (bf) states
            Eigen::Vector6d bf_StateArc1 = transformStateToTargetFrame(
                        currentStateArc1_, currentTimeArc1_, centralBodyRotationModel_ );
            Eigen::Vector6d bf_StateArc2 = transformStateToTargetFrame(
                        currentStateArc2_, currentTimeArc2_, centralBodyRotationModel_ );

            // Definition of dR2dr2
            double r2_norm = currentStateArc2_.segment( 0, 3 ).norm( );
            dR2dr2 << ( currentStateArc2_.segment( 0, 3 ) / r2_norm ).transpose();

            // Definition of dR2dt2
            Eigen::Vector3d r2UnitVector = currentStateArc2_.segment( 0, 3 ) / currentStateArc2_.segment( 0, 3 ).norm();
            dR2dt2 = currentStateArc2_.segment( 3, 3 ).dot( r2UnitVector );

            // Definition of dt2dr2
            Eigen::Vector3d bf_r2UnitVector = bf_StateArc2.segment( 0, 3 ) / bf_StateArc2.segment( 0, 3 ).norm();
            Eigen::Vector3d bf_v1UnitVector = bf_StateArc1.segment( 3, 3 ) / bf_StateArc1.segment( 3, 3 ).norm();
            Eigen::Vector3d bf_v2UnitVector = bf_StateArc2.segment( 3, 3 ) / bf_StateArc2.segment( 3, 3 ).norm();
            Eigen::Matrix3d A_dt2dr2;
            A_dt2dr2 << bf_r2UnitVector.transpose(), bf_v1UnitVector.transpose(), bf_v2UnitVector.transpose();
            Eigen::Vector3d v2Vertical = (bf_StateArc2.segment( 3, 3 ).dot(bf_r2UnitVector))*bf_r2UnitVector;
            Eigen::Vector3d v2InPlane = bf_StateArc2.segment( 3, 3 ) - v2Vertical;
            Eigen::Vector3d b_dt2dr2( 0, 0, - ( 1 / v2InPlane.norm( ) ) );
            Eigen::Vector6d bf_dt2dr2;
            bf_dt2dr2 << ( A_dt2dr2.inverse( ) * b_dt2dr2 ), 0, 0, 0;
            dt2dr2 << transformStateToGlobalFrame(
                          bf_dt2dr2, currentTimeArc2_, centralBodyRotationModel_ ).segment( 0, 3 ).transpose();

            // Definition of dR1dt1
            Eigen::Vector3d r1UnitVector = currentStateArc1_.segment( 0, 3 ) / currentStateArc1_.segment( 0, 3 ).norm();
            dR1dt1 = currentStateArc1_.segment( 3, 3 ).dot( r1UnitVector );

            // Definition of dt1dr2
            Eigen::Vector3d bf_r1UnitVector = bf_StateArc1.segment( 0, 3 ) / bf_StateArc1.segment( 0, 3 ).norm();
            Eigen::Matrix3d A_dt1dr2;
            A_dt1dr2 << bf_r1UnitVector.transpose(), bf_v2UnitVector.transpose(), bf_v1UnitVector.transpose();
            Eigen::Vector3d v1Vertical = (bf_StateArc1.segment( 3, 3 ).dot(bf_r1UnitVector))*bf_r1UnitVector;
            Eigen::Vector3d v1InPlane = bf_StateArc1.segment( 3, 3 ) - v1Vertical;
            Eigen::Vector3d b_dt1dr2( 0, 0, ( 1 / v1InPlane.norm( ) ) );
            Eigen::Vector6d bf_dt1dr2;
            bf_dt1dr2 << ( A_dt1dr2.inverse( ) * b_dt1dr2 ), 0, 0, 0;
            dt1dr2 << transformStateToGlobalFrame(
                          bf_dt1dr2, currentTimeArc2_, centralBodyRotationModel_ ).segment( 0, 3 ).transpose();

            secondArcPartialWrtCurrentPosition << ( dR2dr2 + dR2dt2*dt2dr2 - dR1dt1*dt1dr2 );
//            secondArcPartialWrtCurrentPosition << ( dR2dr2 );

            observationPartialWrtCurrentState.block( 0, 0, 1, 3 ) = secondArcPartialWrtCurrentPosition;

            returnPartial.push_back(
                        std::make_pair( observationPartialWrtCurrentState, currentTimeArc2_ ) );

//            std::cout << "\n dR2dr2: \n" << dR2dr2 << std::endl;

//            std::cout << "\n (dR2dt2): \n" << dR2dt2 << std::endl;
//            std::cout << "\n (dt2dr2): \n" << dt2dr2 << std::endl;
//            std::cout << "\n dR2dt2*dt2dr2: \n" << dR2dt2*dt2dr2 << std::endl;

//            std::cout << "\n (dR1dt1): \n" << dR1dt1 << std::endl;
//            std::cout << "\n (dt1dr2): \n" << dt1dr2 << std::endl;
//            std::cout << "\n dR1dt1*dt1dr2: \n" << dR1dt1*dt1dr2 << std::endl;

//            std::cout << "\n secondArcPartialWrtCurrentPosition: \n" << secondArcPartialWrtCurrentPosition << std::endl;
//            std::cout << "\n observationPartialWrtCurrentState: \n" << observationPartialWrtCurrentState << std::endl;
        }
    }
    return returnPartial;
}

} // namespace observation_partials

} // namespace tudat
