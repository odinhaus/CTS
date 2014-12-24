using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Experiment
{
    [Serializable]
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
                    return 1.0/(1.0 + Math.Exp(-(this.Value + this.Bias))); // sigmoid activation function
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
                    return Sum() * this.LastActivation * (1 - this.LastActivation);
                };
            }

            _Activation = activationFunction;
            _Delta = deltaFunction;
            _Sum = sumFunction;

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
        public Func<double> Activation { get; private set; }
        public Func<double> Delta { get; private set; }
        public Func<double> Sum { get; private set; }

        bool isActivating = false;
        protected Func<double> _Activation
        {
            set
            {
                Activation = () =>
                {
                    if (!isActivating)
                    {
                        isActivating = true;
                        LastActivation = value();
                        isActivating = false;
                    }
                    return LastActivation;
                };
            }
        }
        bool isDeltaing = false;
        protected Func<double> _Delta
        {
            set
            {
                Delta = () =>
                {
                    if (!isDeltaing)
                    {
                        isDeltaing = true;
                        LastDelta = value();
                        isDeltaing = false;
                    }
                    return LastDelta;
                };
            }
        }
        bool isSumming = false;
        protected Func<double> _Sum
        {
            set
            {
                Sum = () =>
                {
                    if (!isSumming)
                    {
                        isSumming = true;
                        LastSum = value();
                        isSumming = false;
                    }
                    return LastSum;
                };
            }
        }

        public double LastActivation { get; private set; }
        public double LastDelta { get; private set; }
        public double LastSum { get; private set; }
        public double Value { get; set; }
        public double Bias { get; set; }

        public override string ToString()
        {
            return string.Format("{0}, {1}, {2}, {3}", Name, Bias, Value, Axons.Count());
        }

        public int RecursionCount { get; set; }

        public void Reset()
        {
            this.Value = 0.0;
            this.LastActivation = 0.0;
            this.LastDelta = 0.0;
            this.RecursionCount = 0;
        }
    }
}
