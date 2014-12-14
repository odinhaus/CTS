using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Experiment
{
    public class XOrNetwork : Network<bool>
    {
        public XOrNetwork(params Neuron[] inputs) :
            base((idx, signal) => signal[idx] ? 1 : 0, inputs)
        {
        }
    }
}