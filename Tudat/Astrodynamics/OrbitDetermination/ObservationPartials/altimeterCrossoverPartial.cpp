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

    Eigen::Matrix< double, 3, 1 > firstArcPartialWrtCurrentPosition;
    Eigen::Matrix< double, 3, 1 > secondArcPartialWrtCurrentPosition;
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
            double rho = currentState_.segment( 0, 3 ).norm( );

            // Beware the MINUS, since xoverObs = |pos(t2)|-|pos(t1)|
            firstArcPartialWrtCurrentPosition << - ( currentInertialPositionPartialWrtParameter *
                                                   (1/rho)*currentState_.segment( 0, 3 ) );

            observationPartialWrtCurrentState.block( 0, 0, 1, 3 ) =
                    firstArcPartialWrtCurrentPosition.transpose();

            returnPartial.push_back(
                        std::make_pair( observationPartialWrtCurrentState, currentTime_ ) );

/*            if( currentTime_ == 1045408937.2387744 || currentTime_ == 1045406388.7029006 )
            {
                std::cout << std::endl << "-- PARTIAL DERIVATIVE TEST (from altimeterCrossoverPartial.cpp) --"
                          << std::endl;
                std::cout << std::setprecision(17) << "for t1: " << currentTime_ << std::endl;
                std::cout << "The state vector s(t1) is: \n" <<
                             currentState_.transpose() << std::endl;
                std::cout << "The rossover observable is: " << std::endl << currentObservation << std::endl;
                std::cout << "The pos.norm is: " << std::endl << rho << std::endl;
                std::cout << "currentInertialPositionPartialWrtParameter: " << std::endl <<
                             currentInertialPositionPartialWrtParameter << std::endl;
                std::cout << "firstArcPartialWrtCurrentPosition: " << std::endl <<
                             firstArcPartialWrtCurrentPosition << std::endl;
                std::cout << "observationPartialWrtCurrentState: " << std::endl <<
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
            double rho = currentState_.segment( 0, 3 ).norm( );
            secondArcPartialWrtCurrentPosition << ( currentInertialPositionPartialWrtParameter *
                                                    (1/rho)*currentState_.segment( 0, 3 ) );

            observationPartialWrtCurrentState.block( 0, 0, 1, 3 ) =
                    secondArcPartialWrtCurrentPosition.transpose();

            returnPartial.push_back(
                        std::make_pair( observationPartialWrtCurrentState, currentTime_ ) );
        }
    }
    return returnPartial;
}

} // namespace observation_partials

} // namespace tudat
