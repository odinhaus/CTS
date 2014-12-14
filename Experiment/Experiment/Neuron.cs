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
        protected string _id = Guid.NewGuid().ToString();

        public Neuron(string name = null, 
            Func<double> activationFunction = null, 
            Func<double> deltaFunction = null, 
            Func<double> sumFunction = null, 
            IEnumerable<Neuron> inputs = null, 
            IEnumerable<Neuron> outputs = null)
        {
            Axons = new List<Synapse>();
            Dendrites = new List<Synapse>();

            if (activationFunction == null)
            {
                activationFunction = () =>
                {
                    LastActivation = 1.0/(1.0 + Math.Exp(-(this.Value + this.Bias))); // sigmoid activation function
                    return LastActivation;
                }; 
            }

            if (sumFunction == null)
            {
                sumFunction = () =>
                {
                    double sum = 0;
                    foreach (var n in this.Axons)
                    {
                        sum += n.Receiver.Delta() * n.Weight;
                    }
                    return sum;
                };
            }

            if (deltaFunction == null)
            {
                deltaFunction = () =>
                {
                    this.LastDelta=  Sum() * this.LastActivation * (1 - this.LastActivation);
                    return this.LastDelta;
                };
            }

            Activation = activationFunction;
            Delta = deltaFunction;
            Sum = sumFunction;

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
        public Func<double> Delta { get; set; }
        public Func<double> Sum { get; set; }

        public double LastActivation { get; protected set; }
        public double LastDelta { get; protected set; }
        public double Value { get; set; }
        public double Bias { get; set; }

        public override string ToString()
        {
            return string.Format("{0}, {1}, {2}, {3}", Name, Bias, Value, Axons.Count());
        }
    }
}
