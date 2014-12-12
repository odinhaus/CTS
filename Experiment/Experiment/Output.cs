using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Experiment
{
    public class Output : Neuron
    {
        public Output(string name = null, IEnumerable<Neuron> inputs = null) 
            : base(name, null, inputs, null)
        {
            this.Activation = () => this.Value + this.Bias;
        }
    }
}
