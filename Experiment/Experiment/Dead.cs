using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Experiment
{
    public class Dead : Neuron
    {
        public Dead(string name = null)
            : base(name, null, null, null, null)
        {
            this.Activation = () =>
            {
                return 0.0;
            };
            this.Delta = () =>
            {
                return 0.0;
            };
            this.Sum = () =>
            {
                return 0.0;
            };
            this.Axons.Clear();
            this.Dendrites.Clear();
        }
    }
}
