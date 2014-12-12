using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Experiment
{
    public class Input : Neuron
    {
        public Input(string name = null, IEnumerable<Neuron> outputs = null)
            : base(name, null, null, outputs)
        {
            this.Activation = () => this.Value + this.Bias;
        }
    }
}
