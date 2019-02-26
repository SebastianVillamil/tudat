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
    std::cout << "You've made it to altimeterCrossoverPartial.cpp!" << std::endl;
    // std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > >
    AltimeterCrossoverPartialReturnType returnPartial;

    // Iterate over all link ends
    for( positionPartialIterator_ = positionPartialList_.begin( ); positionPartialIterator_ != positionPartialList_.end( );
         positionPartialIterator_++ )
    {
        // The current partial relates to the state at arc 1.
        if( positionPartialIterator_->first == observation_models::first_arc_body )
        {
            std::cout << "first loop!" << std::endl;
            currentState_  = states[ 0 ];
            currentTime_ = times[ 0 ];
        }
        // The current partial relates to the state at arc 2.
        else if( positionPartialIterator_->first == observation_models::second_arc_body )
        {
            std::cout << "second loop!" << std::endl;
            currentState_  = states[ 1 ];
            currentTime_ = times[ 1 ];
        }

        Eigen::Matrix< double, 3, Eigen::Dynamic > currentInertialPositionPartialWrtParameter =
                positionPartialIterator_->second->calculatePartialOfPosition(
                                      currentState_ , currentTime_ );

        // IMPLEMENT HERE: PARTIAL DERIVATIVE OF YOUR CROSSOVER OBSERVATION W.R.T. THE CURRENT STATE
        // originally: Eigen::Matrix< double, 1, 3 >
        double rho = currentState_.segment( 0, 3 ).norm( );
        Eigen::Matrix< double, 3, 1 > observationPartialWrtCurrentPosition = (1/rho)*currentState_.segment( 0, 3 );
//        observationPartialWrtCurrentPosition << 1, 1, 1;

        // Set partial output
        returnPartial.push_back(
                    std::make_pair(
                        currentInertialPositionPartialWrtParameter * observationPartialWrtCurrentPosition,
                        currentTime_ ) );
    }

    return returnPartial;
}

} // namespace observation_partials

} // namespace tudat
