using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Experiment
{
    [Serializable]
    public class Input : Neuron
    {
        public Input(string name = null, IEnumerable<Neuron> outputs = null)
            : base(name, null, null, null, outputs)
        {
            this._Activation = () =>
            {
                return this.Value + this.Bias;
            };
            this._Delta = () =>
            {
                return Sum() - this.LastActivation;
            };
        }
    }
}
