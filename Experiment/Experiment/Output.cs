using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Experiment
{
    [Serializable]
    public class Output : Neuron
    {
        public Output(string name = null, IEnumerable<Neuron> inputs = null) 
            : base(name, null, null, null, inputs, null)
        {
            this._Activation = () =>
            {
                return this.Value + this.Bias;
            };
            this._Sum = () => this.Target;
            // use this if you output simple linear activations
            this._Delta = () =>
            {
                return Sum() - this.LastActivation;
            };
            // use this if you output sigmoid activation values
            //this.Delta = () =>
            //{
            //    this.LastDelta = (this.Target - this.LastActivation) * this.LastActivation * (1 - this.LastActivation);
            //    return this.LastDelta;
            //};
        }

        public double Target { get; set; }
    }
}
