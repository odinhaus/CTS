using System;
using System.Collections.Generic;
using System.Linq;


namespace Experiment
{
    public class PID
    {
        public PID(double propGain, double intGain, double dervGain, Func<double> setPoint, Func<double> presVal, Action<double> newVal)
        {
            ProportionalGain = propGain;
            IntegralGain = intGain;
            DerivativeGain = dervGain;
            SetPoint = setPoint;
            PresentValue = presVal;
            NewValue = newVal;
        }

        public double ProportionalGain { get; set; }
        public double IntegralGain { get; set; }
        public double DerivativeGain { get; set; }
        public double Error { get; private set; }
        public double Integral { get; private set; }
        public double Derivative { get; private set; }

        public Func<double> SetPoint { get; private set; }
        public Func<double> PresentValue { get; private set; }
        public Action<double> NewValue { get; private set; }

        public double Update(double deltaT)
        {
            double error = SetPoint() - PresentValue();
            Integral += error*deltaT;
            Derivative = (error - Error) / deltaT;
            double val = ProportionalGain*error + IntegralGain*Integral + DerivativeGain*Derivative;
            NewValue(val);
            Error = error;
            return val;
        }
    }
}
