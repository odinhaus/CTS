using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Experiment
{
    public class Synapse
    {
        public Synapse(Neuron signaler, Neuron receiver)
        {
            Signaler = signaler;
            Receiver = receiver;
        }
        public Neuron Signaler { get; private set; }
        public Neuron Receiver { get; private set; }
        public double Weight { get; set; }
    }
}
