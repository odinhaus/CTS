using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace Experiment
{
    public class Network<T>
    {
        /// <summary>
        /// Create a network from a string which describes the neurons, connections and weights
        /// </summary>
        /// <param name="networkMap">of the form L:I:B [W] L:I:B ... where L = neuron level, starting with 1, 
        /// I = neuron level index, starting with 1, 
        /// B = value bias for the neuron and 
        /// W = connection weight</param>
        /// <param name="maxDepth"></param>
        /// <param name="inputTransform"></param>
        /// <returns></returns>
        public static Network<T> Create(string networkMap, int maxDepth, Func<int, T, double> inputTransform)
        {
            Dictionary<string, Neuron> neurons = new Dictionary<string, Neuron>();

            string pat = @"(?=((?<connection>(?<source>\d:\d(:(?<sourceBias>[-\d\.]*))*)\s*\[(?<weight>[-\d\.]*)\]\s*(?<dest>\d:\d(:(?<destBias>[-\d\.]*))*))))";
            Regex r = new Regex(pat);
            Match m = r.Match(networkMap);
            while (m.Success)
            {
                string source = m.Groups["source"].Value;
                string dest = m.Groups["dest"].Value;
                string weight = m.Groups["weight"].Value;
                double sourceBias = 0, destBias = 0;
                if (!string.IsNullOrEmpty(m.Groups["sourceBias"].Value))
                {
                    sourceBias = double.Parse(m.Groups["sourceBias"].Value);
                    source = source.Replace(":" + m.Groups["sourceBias"].Value, "");
                }
                if (!string.IsNullOrEmpty(m.Groups["destBias"].Value))
                {
                    destBias = double.Parse(m.Groups["destBias"].Value);
                    dest = dest.Replace(":" + m.Groups["destBias"].Value, "");
                }

                Neuron sourceN = null, destN = null;

                if (!neurons.TryGetValue(source, out sourceN))
                {
                    if (source.StartsWith("1"))
                    {
                        sourceN = new Input(source);
                    }
                    else
                    {
                        sourceN = new Neuron(source);
                    }
                    sourceN.Bias = sourceBias;
                    neurons.Add(source, sourceN);
                }

                if (!neurons.TryGetValue(dest, out destN))
                {
                    if (dest.StartsWith(maxDepth + ":"))
                    {
                        destN = new Output(dest);
                    }
                    else
                    {
                        destN = new Neuron(dest);
                    }

                    destN.Bias = destBias;
                    neurons.Add(dest, destN);
                }

                var s = new Synapse(sourceN, destN) {Weight = double.Parse(weight)};
                sourceN.Axons.Add(s);
                destN.Dendrites.Add(s);

                m = m.NextMatch();
            }

            return new Network<T>(inputTransform, neurons.Values.Where(n => n is Input).ToArray()) { _maxDepth = maxDepth };
        }

        private Neuron[] _inputs;
        private Neuron[] _outputs;
        private Func<int, T, double> _inputTransform;
        private int _maxDepth;
        public Network(Func<int, T, double> inputTransform, params Neuron[] inputs)
        {
            _inputs = inputs;
            var outputs = new List<Neuron>();
            Traverse((n, level) => true, (n, level) =>  n is Output, outputs.Add);
            _outputs = outputs.ToArray();
            _inputTransform = inputTransform;
        }

        public double[] Update(T inputSignal)
        {
            SetInputs(inputSignal);
            Propagate();
            return ReturnOutputs();
        }

        private double[] ReturnOutputs()
        {
            double[] outputs = new double[_outputs.Length];
            for (int i = 0; i < _outputs.Length; i++)
            {
                Output n = _outputs[i] as Output;
                n.Value = 0;
                foreach (var input in n.Dendrites)
                {
                    input.Receiver.Value += input.Signaler.Activation() * input.Weight;
                }

                outputs[i] = _outputs[i].Activation();
            }
            return outputs;
        }

        private void Propagate()
        {
            Traverse((n, level) => level < _maxDepth, (n, level) => !(n is Input), (n) =>
            {
                n.Value = 0;
                foreach (var input in n.Dendrites)
                {
                    input.Receiver.Value += input.Signaler.Activation() * input.Weight;
                }
            }, true);
        }

        private void Traverse(Func<Neuron, int, bool> traversalSelector, Func<Neuron, int, bool> actionSelector, Action<Neuron> action, bool allowMultipleVisits = false)
        {
            _visited.Clear();
            int level = 1;
            foreach (var n in _inputs)
            {
                if (traversalSelector(n, level))
                    TraverseRecurse(n, traversalSelector, actionSelector, action, allowMultipleVisits, ref level);
            }
        }

        private HashSet<Neuron> _visited = new HashSet<Neuron>();
        private void TraverseRecurse(Neuron n, Func<Neuron, int, bool> traversalSelector, Func<Neuron, int, bool> actionSelector, Action<Neuron> action, bool allowMultipleVisits, ref int level)
        {
            if (allowMultipleVisits || !_visited.Contains(n))
            {
                if (actionSelector(n, level))
                {
                    action(n);
                }
                _visited.Add(n);
            }

            foreach (var nn in n.Axons)
            {
                level++;
                if (traversalSelector(nn.Receiver, level))
                {
                    TraverseRecurse(nn.Receiver, traversalSelector, actionSelector, action, allowMultipleVisits, ref level);
                }
                level--;
            }
        }

        protected void SetInputs(T inputSignal)
        {
            _inputIndex = 0;
            Traverse((n, level) => n is Input, (n, level) => true, (n) => SetInputs(n, inputSignal), false);
        }

        private int _inputIndex = 0;
        protected virtual void SetInputs(Neuron n, T inputSignal)
        {
            n.Value = _inputTransform(_inputIndex, inputSignal);
            _inputIndex++;
        }

        public virtual void Train(T[] inputSignals, double[] outputSignals)
        {
        }

        public Neuron[] Inputs { get { return _inputs; } }
        public double[] Outputs { get; set; }
    }
}
