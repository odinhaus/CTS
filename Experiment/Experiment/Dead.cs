using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Experiment
{
    [Serializable]
    public class Dead : Neuron
    {
        public Dead(string name = null)
            : base(name, null, null, null, null)
        {
            this._Activation = () =>
            {
                return 0.0;
            };
            this._Delta = () =>
            {
                return 0.0;
            };
            this._Sum = () =>
            {
                return 0.0;
            };
            this.Axons.Clear();
            this.Dendrites.Clear();
        }
    }
}
