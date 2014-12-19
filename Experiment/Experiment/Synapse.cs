using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Experiment
{
    [Serializable]
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

        public override string ToString()
        {
            return string.Format("{0} <-[{2}]-> {1}", Signaler.Name, Receiver.Name, Weight.ToString("#0.0####"));
        }
    }
}
