/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ALTIMETRYCROSSOVEROBSERVATIONMODEL_H
#define TUDAT_ALTIMETRYCROSSOVEROBSERVATIONMODEL_H

#include <map>

#include <functional>
#include <boost/make_shared.hpp>

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType = double,
          typename TimeType = double >
class AltimetryCrossoverObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:    
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;

    AltimetryCrossoverObservationModel(
            const std::function< StateType( const TimeType ) > firstArcBodyStateFunction,
            const std::function< StateType( const TimeType ) > secondArcBodyStateFunction,
            const std::string& centralBody,
            const std::map< double, double >& crossoverTimes,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( altimetry_crossover, observationBiasCalculator ),
        firstArcBodyStateFunction_( firstArcBodyStateFunction ), secondArcBodyStateFunction_( secondArcBodyStateFunction ),
        centralBody_( centralBody ), crossoverTimes_( crossoverTimes ){ }

    //! Destructor
    ~AltimetryCrossoverObservationModel( ){ }

    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime )
    {
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > linkEndStates;

        return computeIdealObservationsWithLinkEndData(
                    time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates );
    }

    //! Function to compute altimetry crossover observables without any corrections.
    /*!
     *  Function to compute one-way range  observable without any corrections, i.e. the true physical range as computed
     *  from the defined link ends. Note that this observable does include light-time
     *  corrections, which represent physically true corrections. It does not include e.g. system-dependent measurement
     *  errors, such as biases or clock errors.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which observable is to be evaluated.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Ideal one-way range observable.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
                    const TimeType time,
                    const LinkEndType linkEndAssociatedWithTime,
                    std::vector< double >& linkEndTimes,
                    std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
    {
        linkEndTimes.clear( );
        linkEndStates.clear( );

        ObservationScalarType crossoverAltimetryObservation = TUDAT_NAN;
        // Add if checks for the existency of t1 & t2
        TimeType firstArcTime = TUDAT_NAN;
        TimeType secondArcTime = TUDAT_NAN;
        StateType firstArcState;
        StateType secondArcState;

        switch ( linkEndAssociatedWithTime )
        {
        case first_arc_body:
            firstArcTime = time;
            secondArcTime = crossoverTimes_[ time ];
            firstArcState = firstArcBodyStateFunction_( firstArcTime );
            secondArcState = secondArcBodyStateFunction_( secondArcTime );

            crossoverAltimetryObservation = ( secondArcState.segment( 0, 3 ).norm( ) -
                                              firstArcState.segment( 0, 3 ).norm( ) );
            break;

//        case second_arc_body:
//            std::string ArcErrorMessage = "Error, linkEndAssociatedWithTime: " +
//                    std::to_string( linkEndAssociatedWithTime ) + "not yet implemented.";
//            throw std::runtime_error( ArcErrorMessage );
//            break;

        default:
            std::string errorMessage = "Error, cannot have link end type: " +
                    std::to_string( linkEndAssociatedWithTime ) + "for altimetry crossover";
            throw std::runtime_error( errorMessage );
        }

        linkEndTimes.push_back( static_cast< double >( firstArcTime ) );
        linkEndTimes.push_back( static_cast< double >( secondArcTime ) );

        linkEndStates.push_back( firstArcState.template cast< double >(  ) );
        linkEndStates.push_back( secondArcState.template cast< double >(  ) );

        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << crossoverAltimetryObservation ).finished( );
    }

    //! Function to compute one-way range observable without any corrections.
    /*!
     *  Function to compute one-way range  observable without any corrections, i.e. the true physical range as computed
     *  from the defined link ends. Note that this observable does include light-time
     *  corrections, which represent physically true corrections. It does not include e.g. system-dependent measurement
     *  errors, such as biases or clock errors.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which observable is to be evaluated.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Ideal one-way range observable.
     */
//    std::vector< double > computeCrossoverData( const TimeType time )
    Eigen::VectorXd computeCrossoverData( const TimeType time )
    {
        TimeType firstArcTime = time, secondArcTime = crossoverTimes_[ time ];
        StateType firstArcState = firstArcBodyStateFunction_( firstArcTime );
        StateType secondArcState = secondArcBodyStateFunction_( secondArcTime );

        ObservationScalarType crossoverAltimetryObservation;
        crossoverAltimetryObservation = ( secondArcState.segment( 0, 3 ).norm( ) -
                                          firstArcState.segment( 0, 3 ).norm( ));

        double PosDotVel1 = (firstArcState.segment( 0, 3 ).transpose()*(firstArcState.segment( 3, 3 )));
        double RadialVelt1 = PosDotVel1/(firstArcState.segment( 0, 3 ).norm( ));

        double PosDotVel2 = (secondArcState.segment( 0, 3 ).transpose()*(secondArcState.segment( 3, 3 )));
        double RadialVelt2 = PosDotVel2/(secondArcState.segment( 0, 3 ).norm( ));

        Eigen::VectorXd XoverDataVector( 8 );
        XoverDataVector[ 0 ] = secondArcTime;
        XoverDataVector[ 1 ] = ( firstArcState.segment( 0, 3 ).norm( ) );
        XoverDataVector[ 2 ] = ( secondArcState.segment( 0, 3 ).norm( ) );
        XoverDataVector[ 3 ] = crossoverAltimetryObservation;
        XoverDataVector[ 4 ] = RadialVelt1;
        XoverDataVector[ 5 ] = RadialVelt2;
        XoverDataVector[ 6 ] = ( firstArcState.segment( 3, 3 ).norm( ) );
        XoverDataVector[ 7 ] = ( secondArcState.segment( 3, 3 ).norm( ) );

/*        if( time == 1045408937.2387744 || time == 1045406388.7029006 )
        {
            std::cout << std::endl << "-- PARTIAL DERIVATIVE TEST (from altimeterCrossoverPartial.cpp) --"
                      << std::endl;
            std::cout << std::setprecision(17) << "s(t1): \n" << firstArcState << std::endl;
            std::cout << "firstArcState.segment( 3, 3 ).norm(): \n" <<
                         firstArcState.segment( 3, 3 ).norm( ) << std::endl;
        } // */
        return ( XoverDataVector );
    }

    Eigen::Vector3d computeDt2Dr2( const TimeType time )
    {
        TimeType firstArcTime = time, secondArcTime = crossoverTimes_[ time ];
        StateType firstArcState = firstArcBodyStateFunction_( firstArcTime );
        StateType secondArcState = secondArcBodyStateFunction_( secondArcTime );

        Eigen::Vector3d XoverPartialsVector;

        Eigen::Vector3d r2UnitVector = secondArcState.segment( 0, 3 )/secondArcState.segment( 0, 3 ).norm();
        Eigen::Vector3d v1UnitVector = firstArcState.segment( 3, 3 )/firstArcState.segment( 3, 3 ).norm();
        Eigen::Vector3d v2UnitVector = secondArcState.segment( 3, 3 )/secondArcState.segment( 3, 3 ).norm();
        Eigen::Matrix3d A;
        A << r2UnitVector.transpose(), v1UnitVector.transpose(), v2UnitVector.transpose();
//        std::cout << std::setprecision(17) << "r2Vector: \n" << secondArcState.segment( 0, 3 ) << std::endl;
//        std::cout << "r2UnitVector: \n" << r2UnitVector << std::endl;
//        std::cout << "v1UnitVector: \n" << v1UnitVector << std::endl;
//        std::cout << "v2UnitVector: \n" << v2UnitVector << std::endl;
//        std::cout << "A: \n" << A << std::endl;

//        double r2_norm = secondArcState.segment( 0, 3 ).norm( );
//        Eigen::Matrix< double, 1, 3 > dR2dr2;
//        dR2dr2 << ( (1/r2_norm) * secondArcState.segment( 0, 3 ) ).transpose();
//        std::cout << "dR2dr2: \n" << dR2dr2 << std::endl;

        Eigen::Vector3d v2UnitVectorHor = secondArcState.segment( 3, 3 ) -
                secondArcState.segment( 3, 3 ).cwiseProduct(r2UnitVector);
        Eigen::Vector3d b( 0, 0, -( 1 / v2UnitVectorHor.norm( )) );
//        std::cout << std::setprecision(17) << "A: \n" << A << std::endl;
//        std::cout << "A determinant: \n" << A.determinant() << std::endl;
//        std::cout << "A inverse: \n" << A.inverse() << std::endl;
//        std::cout << "A inverse times b \n" << A.inverse()*b << std::endl;

        return ( XoverPartialsVector );
    }

private:

    std::function< StateType( const TimeType ) > firstArcBodyStateFunction_;

    std::function< StateType( const TimeType ) > secondArcBodyStateFunction_;

    std::string centralBody_;

    std::map< double, double > crossoverTimes_;

};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_ALTIMETRYCROSSOVEROBSERVATIONMODEL_H
