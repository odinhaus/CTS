using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Experiment
{
    public class Neuron
    {
        private static int _index = 0;
        public Neuron(string name = null, Func<double> activationFunction = null, IEnumerable<Neuron> inputs = null, IEnumerable<Neuron> outputs = null)
        {
            Axons = new List<Synapse>();
            Dendrites = new List<Synapse>();

            if (activationFunction == null)
            {
                activationFunction = () =>
                {
                    return 1.0/(1.0 + Math.Exp(-(this.Value + this.Bias))); // sigmoid activation function
                }; 
            }

            Activation = activationFunction;
            if (inputs != null)
            {
                foreach (var n in inputs)
                {
                    Dendrites.Add(new Synapse(n, this));
                }
            }
            if (outputs != null)
            {
                foreach (var n in outputs)
                {
                    Dendrites.Add(new Synapse(this, n));
                }
            }

            if (name == null)
                name = (_index++).ToString();
            this.Name = name;
        }

        public string Name { get; private set; }
        public IList<Synapse> Axons { get; set; }
        public IList<Synapse> Dendrites { get; set; }
        public Func<double> Activation { get; set; }

        public double Value { get; set; }
        public double Bias { get; set; }
    }
}
