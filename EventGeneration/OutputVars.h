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
        kXPos,
        kYPos,
        kZPos,
        kXVel,
        kYVel,
        kZVel,
        kXAcc,
        kYAcc,
        kZAcc,
        kBField
    };
}

#endif