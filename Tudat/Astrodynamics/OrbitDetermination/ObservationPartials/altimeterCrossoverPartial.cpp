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

    Eigen::Matrix< double, 1, 3 > firstArcPartialWrtCurrentPosition;  // 1, 3 as dh/dr1
    Eigen::Matrix< double, 1, 3 > secondArcPartialWrtCurrentPosition;  // 1, 3 as dh/dr2

    double dR2dt2;
    Eigen::Matrix< double, 1, 3 > dt2dr1;  // 1, 3 due to d/dr1
    Eigen::Matrix< double, 1, 3 > dR1dr1;  // 1, 3 due to .transpose( )
    double dR1dt1;
    Eigen::Matrix< double, 1, 3 > dt1dr1;  // 1, 3 due to d/dr1

    Eigen::Matrix< double, 1, 3 > dR2dr2;  // 1, 3 due to .transpose( )
//    Eigen::Matrix< double, 3, 1 > dR2dt2;
    Eigen::Matrix< double, 1, 3 > dt2dr2;  // 1, 3 due to d/dr2
//    Eigen::Matrix< double, 3, 1 > dR1dt1;
    Eigen::Matrix< double, 1, 3 > dt1dr2;  // 1, 3 due to d/dr2

    Eigen::MatrixXd observationPartialWrtCurrentState = Eigen::MatrixXd::Zero( 1, 6 );

    // Iterate over all link ends
    for( positionPartialIterator_ = positionPartialList_.begin( ); positionPartialIterator_ != positionPartialList_.end( );
         positionPartialIterator_++ )
    {
        if( positionPartialIterator_->first == observation_models::first_arc_body )
        {
            currentState_  = states[ 0 ];
            currentTime_ = times[ 0 ];

            Eigen::Matrix< double, 3, Eigen::Dynamic > currentInertialPositionPartialWrtParameter =
                    positionPartialIterator_->second->calculatePartialOfPosition(
                                          currentState_ , currentTime_ );

            // def of dR2dt2
            Eigen::Vector3d r2UnitVector = states[ 1 ].segment( 0, 3 ) / states[ 1 ].segment( 0, 3 ).norm();
            dR2dt2 = states[ 1 ].segment( 0, 3 ).dot( r2UnitVector );

            // def of dt2dr1
            Eigen::Vector3d r1UnitVector = states[ 0 ].segment( 0, 3 ) / states[ 0 ].segment( 0, 3 ).norm();
            Eigen::Vector3d v1UnitVector = states[ 0 ].segment( 3, 3 ) / states[ 0 ].segment( 3, 3 ).norm();
            Eigen::Vector3d v2UnitVector = states[ 1 ].segment( 3, 3 ) / states[ 1 ].segment( 3, 3 ).norm();
            Eigen::Matrix3d A_dt2dr1;
            A_dt2dr1 << r1UnitVector.transpose(), v1UnitVector.transpose(), v2UnitVector.transpose();
            Eigen::Vector3d v2InPlane = states[ 1 ].segment( 3, 3 ) -
                    states[ 1 ].segment( 3, 3 ).cwiseProduct(r2UnitVector);
            Eigen::Vector3d b_dt2dr1( 0, 0, ( 1 / v2InPlane.norm( ) ) );
            dt2dr1 << ( A_dt2dr1.inverse( ) * b_dt2dr1 );

            // def of dR1dr1
            double r1_norm = states[ 0 ].segment( 0, 3 ).norm( );
            dR1dr1 << ( (1/r1_norm)*states[ 0 ].segment( 0, 3 ) ).transpose();

            // def of dR1dt1
            dR1dt1 = states[ 0 ].segment( 0, 3 ).dot( r1UnitVector );

            // def of dt1dr1
            Eigen::Matrix3d A_dt1dr1;
            A_dt1dr1 << r1UnitVector.transpose(), v2UnitVector.transpose(), v1UnitVector.transpose();
            Eigen::Vector3d v1InPlane = states[ 0 ].segment( 3, 3 ) -
                    states[ 0 ].segment( 3, 3 ).cwiseProduct( r1UnitVector );
            Eigen::Vector3d b_dt1dr1( 0, 0, - ( 1 / v1InPlane.norm( ) ) );
            dt1dr1 << ( A_dt1dr1.inverse( ) * b_dt1dr1 );

            firstArcPartialWrtCurrentPosition << ( dR2dt2*dt2dr1 - dR1dr1 - dR1dt1*dt1dr1 );
//            firstArcPartialWrtCurrentPosition << ( - dR1dr1 );
            observationPartialWrtCurrentState.block( 0, 0, 1, 3 ) = firstArcPartialWrtCurrentPosition;

            // Beware the MINUS, since xoverObs = |pos(t2)|-|pos(t1)|
//            firstArcPartialWrtCurrentPosition << - ( (1/rho)*currentState_.segment( 0, 3 ) );

//            observationPartialWrtCurrentState.block( 0, 0, 1, 3 ) =
//                    firstArcPartialWrtCurrentPosition.transpose() * 1.0; //  * -1.0 as ad-hoc fix

            returnPartial.push_back(
                        std::make_pair( observationPartialWrtCurrentState, currentTime_ ) );

/*            if( currentTime_ == 1045408937.2387744 || currentTime_ == 1045406388.7029006 )
            {
                std::cout << std::endl << "-- PARTIAL DERIVATIVE TEST (from altimeterCrossoverPartial.cpp) --"
                          << std::endl;
                std::cout << std::setprecision(17) << "for t1: " << currentTime_ << std::endl;
                std::cout << "\n The state vector s(t1) is: \n" <<
                             currentState_.transpose() << std::endl;

                std::cout << "\n The rossover observable is: " << currentObservation << std::endl;

                std::cout << "\n The pos.norm is: " << rho << std::endl;

                std::cout << "\n currentInertialPositionPartialWrtParameter: " << std::endl <<
                             currentInertialPositionPartialWrtParameter << std::endl;

                std::cout << "\n firstArcPartialWrtCurrentPosition: " << std::endl <<
                             firstArcPartialWrtCurrentPosition << std::endl;

                std::cout << "\n observationPartialWrtCurrentState: " << std::endl <<
                             observationPartialWrtCurrentState << std::endl << std::endl;
            } // */
        }
        else if( positionPartialIterator_->first == observation_models::second_arc_body )
        {
            currentState_  = states[ 1 ];
            currentTime_ = times[ 1 ];

            Eigen::Matrix< double, 3, Eigen::Dynamic > currentInertialPositionPartialWrtParameter =
                    positionPartialIterator_->second->calculatePartialOfPosition(
                                          currentState_ , currentTime_ );

            // def of dR2dr2
            double r2_norm = states[ 1 ].segment( 0, 3 ).norm( );
            dR2dr2 << ( (1/r2_norm)*states[ 1 ].segment( 0, 3 ) ).transpose();

            // def of dR2dt2
            Eigen::Vector3d r2UnitVector = states[ 1 ].segment( 0, 3 ) / states[ 1 ].segment( 0, 3 ).norm();
            dR2dt2 = states[ 1 ].segment( 0, 3 ).dot( r2UnitVector );

            // def of dt2dr2
            Eigen::Vector3d v1UnitVector = states[ 0 ].segment( 3, 3 ) / states[ 0 ].segment( 3, 3 ).norm();
            Eigen::Vector3d v2UnitVector = states[ 1 ].segment( 3, 3 ) / states[ 1 ].segment( 3, 3 ).norm();
            Eigen::Matrix3d A_dt2dr2;
            A_dt2dr2 << r2UnitVector.transpose(), v1UnitVector.transpose(), v2UnitVector.transpose();
            Eigen::Vector3d v2InPlane = states[ 1 ].segment( 3, 3 ) -
                    states[ 1 ].segment( 3, 3 ).cwiseProduct( r2UnitVector );
            Eigen::Vector3d b_dt2dr2( 0, 0, - ( 1 / v2InPlane.norm( ) ) );
            dt2dr2 << ( A_dt2dr2.inverse( ) * b_dt2dr2 );

            // def of dR1dt1
            Eigen::Vector3d r1UnitVector = states[ 1 ].segment( 0, 3 ) / states[ 1 ].segment( 0, 3 ).norm();
            dR1dt1 = states[ 0 ].segment( 0, 3 ).dot( r1UnitVector );

            // def of dt1dr2
            Eigen::Matrix3d A_dt1dr2;
            A_dt1dr2 << r1UnitVector.transpose(), v2UnitVector.transpose(), v1UnitVector.transpose();
            Eigen::Vector3d v1InPlane = states[ 0 ].segment( 3, 3 ) -
                    states[ 0 ].segment( 3, 3 ).cwiseProduct( r1UnitVector );
            Eigen::Vector3d b_dt1dr2( 0, 0, ( 1 / v1InPlane.norm( ) ) );
            dt1dr2 << ( A_dt1dr2.inverse( ) * b_dt1dr2 );

            secondArcPartialWrtCurrentPosition << ( dR2dr2 + dR2dt2*dt2dr2 + dR1dt1*dt1dr2 );
//            secondArcPartialWrtCurrentPosition << ( dR1dt1*dt1dr2 );
            observationPartialWrtCurrentState.block( 0, 0, 1, 3 ) = secondArcPartialWrtCurrentPosition;

//            secondArcPartialWrtCurrentPosition << ( (1/rho)*currentState_.segment( 0, 3 ) );

//            observationPartialWrtCurrentState.block( 0, 0, 1, 3 ) =
//                    secondArcPartialWrtCurrentPosition.transpose() * 1.0; //  * -1.0 as ad-hoc fix

            returnPartial.push_back(
                        std::make_pair( observationPartialWrtCurrentState, currentTime_ ) );

/*            if( times[ 0 ] == 1045408937.2387744 || times[ 0 ] == 1045406388.7029006 )
            {
                std::cout << std::endl << "-- PARTIAL DERIVATIVE TEST (from altimeterCrossoverPartial.cpp) --"
                          << std::endl;
                std::cout << std::setprecision(17) << "for t2: " << currentTime_ << std::endl;
                std::cout << "\n The state vector s(t2) is: \n" <<
                             currentState_.transpose() << std::endl;

                std::cout << "\n The rossover observable is: " << currentObservation << std::endl;

                std::cout << "\n The pos.norm is: " << rho << std::endl;

                std::cout << "\n currentInertialPositionPartialWrtParameter: " << std::endl <<
                             currentInertialPositionPartialWrtParameter << std::endl;

                std::cout << "\n secondArcPartialWrtCurrentPosition: " << std::endl <<
                             secondArcPartialWrtCurrentPosition << std::endl;

                std::cout << "\n observationPartialWrtCurrentState: " << std::endl <<
                             observationPartialWrtCurrentState << std::endl << std::endl;
            } // */
        }
    }
    return returnPartial;
}

} // namespace observation_partials

} // namespace tudat
