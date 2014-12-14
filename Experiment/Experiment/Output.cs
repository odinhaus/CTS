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
            : base(name, null, null, null, inputs, null)
        {
            this.Activation = () =>
            {
                LastActivation = this.Value + this.Bias;
                return LastActivation;
            };
            this.Sum = () => this.Target;
            this.Delta = () =>
            {
                LastDelta = Sum() - this.LastActivation;
                return LastDelta;
            };
        }

        public double Target { get; set; }
    }
}
