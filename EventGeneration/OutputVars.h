/*
    OutputVars.h

    Variables used to specify the desired output from a simulation to file
    This can be electron kinematics or output voltages
*/

#ifndef OUTPUT_VARS_H
#define OUTPUT_VARS_H

namespace rad
{
    enum OutputVar
    {
        kPos,
        kVel,
        kAcc,
        kBField,
        kParticleTime,

        // E and B fields at the antenna
        kAntEField,
        kAntBField,
        // Induced voltage at an antenna
        kAntVoltage,

        // Sample time at some hypothetical digitizer
        kSampleTime,
        // In-phase and quadrature components of some voltage
        kVI,
        kVQ
    };
}

#endif